import re
import os
import shutil
import gzip
import time
import uuid
import subprocess
import json
import argparse
import multiprocessing
import uuid

import numpy as np
import pandas as pd

import zmq

try:
    import pyarrow.feather
except ImportError:
    pass


from .configuration import Configuration
from . import bin
from . import Rpkg
from .serdef import serialize
from .generators import *


def start_broker():
    dealer_path = os.path.join(os.path.dirname(bin.__file__), "dealer")

    proc = subprocess.Popen([dealer_path], stdout=subprocess.PIPE)

    lines = []
    while True:
        lines.append(proc.stdout.readline().decode("utf-8"))
        if lines[-1] == "<>\n":
            break

    in_addr = lines[1].rstrip().split("=")[1]
    out_addr = lines[2].rstrip().split("=")[1]

    return proc, {"in": in_addr, "out": out_addr}


def start_worker(worker_id, out_addr, monitor_addr, worker_script):
    worker_base = os.path.join(os.path.dirname(Rpkg.__file__), "R")

    cmd = [
        "Rscript", "--vanilla",
        "-e", "source('{}')".format(os.path.join(worker_base, "serdef.R")),
        "-e", "source('{}')".format(os.path.join(worker_base, "worker.R")),
        "-e", "source('{}')".format(worker_script),
        str(worker_id), out_addr, monitor_addr
    ]

    proc = subprocess.Popen(cmd)

    return proc, {}


def send_data(socket, metadata, data):
    socket.send_multipart([
        json.dumps(metadata).encode("utf-8"), serialize(data)
    ])


def monitor(state):
    ctx = zmq.Context()
    monitor_sock = ctx.socket(zmq.REP)

    port = monitor_sock.bind_to_random_port("tcp://0.0.0.0")
    state["monitor_addr"] = f"tcp://0.0.0.0:{port}"
    state["ready"] = 0
    state["done_pushing"] = 0

    print(f"Monitor ready ({port})")

    while True:
        message = monitor_sock.recv()

        if message.decode("utf-8").startswith("ready"):
            state["ready"] += 1
            monitor_sock.send(b"OK")

        elif message.decode("utf-8").startswith("is_done?"):
            monitor_sock.send(bytes(state["done_pushing"]))

        else:
            print("Python got weird message: ", message)
            monitor_sock.send(b"WEIRD")


def collect_and_teardown_analysis(output_prefix):
    """Teardown a sandboxed analysis by collecting worker results."""
    # Pattern to define worker artifacts.
    # The format is:
    #  results_worker_ID_output_filename.ext
    pat = (r"^results_worker_(?P<worker_id>\d)_"
            "(?P<target>\w+)"
            "\.(?P<ext>[a-zA-Z]+)$")

    pat = re.compile(pat)

    # Infer the analysis id from current directory.
    cwd = os.getcwd()
    analysis_id = os.path.split(cwd)[-1]

    # There may be multiple file to collect, we detect it.
    worker_files = [fn for fn in os.listdir() if re.match(pat, fn)]
    artifacts = [re.match(pat, fn).groupdict() for fn in worker_files]

    output_files = {}

    def close_files():
        for f in output_files.values():
            f.close()

    try:
        for o in artifacts:
            target = o["target"]
            fn = os.path.join("..", "{}_{}.{}.gz".format(
                output_prefix, o["target"], o["ext"])
            )

            if target not in output_files:
                print(f"Opening output file: {fn}")
                output_files[target] = gzip.open(fn, "wb")

    except Exception as e:
        close_files()
        raise e

    # Detect and write headers.
    header_pat = re.compile(r"^header_(?P<target>\w+)\.[a-zA-Z]+$")
    for fn in os.listdir():
        match = re.match(header_pat, fn)
        if match:
            # We found a header and we need to write it.
            target = match.groupdict()["target"]

            if target not in output_files:
                print("WARNING ignoring header for unknown target: '{}'"
                      "".format(target))
                continue

            with open(fn, "rb") as f:
                output_files[target].write(f.read())

    # Aggregate output files appropriately.
    for filename, artifact in zip(worker_files, artifacts):
        out = output_files[artifact["target"]]
        with open(filename, "rb") as f:
            out.write(f.read())

    close_files()


def setup_analysis(configuration):
    """Setup a directory and configuration files for the analysis."""
    # Create an analysis id and change to the directory for the current
    # analysis.
    # This is to create some kind of sandboxing where workers can write
    # output files during analyses.
    analysis_id = str(uuid.uuid4())[:5]
    os.mkdir(analysis_id)

    filename = os.path.join(analysis_id, "analysis_configuration.json")
    with open(filename, "w") as f:
        f.write(configuration.to_json())

    return analysis_id


def run_phewas(configuration, n_workers, data_generators, output_prefix,
               worker_script):

    context = zmq.Context()

    # Setup the zmq system.
    processes = []
    workers = []
    try:
        # Start the monitor.
        monitor_manager = multiprocessing.Manager()
        monitor_state = monitor_manager.dict()
        monitor_t = multiprocessing.Process(
            target=monitor,
            args=(monitor_state, ),
            daemon=True
        )
        processes.append(monitor_t)

        print("Py: Starting monitor process")
        monitor_t.start()
        print("Py: ok")

        while True:
            # Wait for the monitor to be ready.
            if "monitor_addr" in monitor_state:
                break

        # Start the broker.
        broker, broker_info = start_broker()
        processes.append(broker)
        print(broker_info)

        # Start workers.
        for i in range(n_workers):
            worker, _ = start_worker(
                i, broker_info["out"], monitor_state["monitor_addr"],
                worker_script
            )
            workers.append(worker)
            processes.append(worker)

        # Wait for all the workers to have signaled ready.
        while True:
            if monitor_state["ready"] == n_workers:
                print("Py: All workers are ready.")
                break

        # Connect to the broker and send data for the workers.
        socket = context.socket(zmq.PUSH)
        socket.connect(broker_info["in"])

        for generator in data_generators:
            try:
                for metadata, data in generator:
                    send_data(socket, metadata, data)
            except MissingDataSourceError:
                print(f"Skipping phenotypes generated by '{generator}' "
                       "(missing data source).")

        monitor_state["done_pushing"] = True
        print("Py: Setting done pushing")

        # Poll workers until they are done.
        while True:
            worker_polls = [p.poll() for p in workers]
            if None in worker_polls:
                # Some workers are not done, wait a bit.
                time.sleep(5)

            else:
                if all([ret == 0 for ret in worker_polls]):
                    print("Py: All workers completed successfully.")
                else:
                    print("Py: Some workers exited with an error.")

                break

    finally:
        print("Py: Killing remaining processes.")
        for proc in processes:
            if type(proc) is subprocess.Popen:
                proc.kill()

            elif type(proc) is multiprocessing.Process:
                proc.terminate()

            else:
                raise ValueError()


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cpus",
        default=multiprocessing.cpu_count(),
        type=int
    )

    parser.add_argument(
        "--configuration",
        help="Path to the script that defines the python configuration object.",
        required=True
    )

    parser.add_argument(
        "--analysis-type", "-a",
        choices=["ICD10_3CHAR", "ICD10_BLOCK", "ICD10_RAW",
                 "CONTINUOUS_VARIABLE", "SELF_REPORTED", "CV_ENDPOINTS",
                 "PHECODES"],
        action="append"
    )

    parser.add_argument(
        "--output", "-o",
        help="Output prefix",
        default="results"
    )

    # Sex-based filtering
    parser.add_argument(
        "--sex-path",
        help="Path to the sex data in the databundle. For example "
             "'worker_data.male' where worker_data is the name of the data "
             "source (this information is used for subsetting).",
        default=None
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--male-only", action="store_true")
    group.add_argument("--female-only", action="store_true")

    return parser.parse_args()


def detect_sample_sex(configuration, sex_path):
    # Look for sex_column in the covariables.
    sex = None
    cache = configuration.get_cache()

    ds_name, sex_col = sex_path.split(".")

    if ds_name not in cache:
        raise ValueError("Could not find data source '{}' which is supposed "
                         "to hold sex information.".format(ds_name))

    if sex_col not in cache[ds_name]:
        raise ValueError("Could not find sex column '{}' in data source '{}'"
                         "".format(sex_col, ds_name))

    sex = cache[ds_name][["sample_id", sex_col]]

    values = sex[sex_col].dropna().unique()

    if set(values) == {0, 1}:
        # Sex is encoded as a boolean vector. Try to infer meaning from
        # column name.
        if sex_col == "male":
            # Do Nothing
            pass

        elif sex_col == "female":
            sex["male"] = 1 - sex["female"]

        else:
            raise ValueError(
                "Can't infer sex meaning from column '{}'. Name the sex "
                "column 'male', 'female' or use male = 1, female = 2, "
                "missing = 0."
                "".format(sex_col)
            )

    elif set(values) == {0, 1, 2}:
        print("Assuming 1=male, 2=female, 0=missing for column '{}'."
              "".format(sex_col))

        sex["male"] = np.nan
        sex.loc[sex[sex_col] == 1, "male"] = 1
        sex.loc[sex[sex_col] == 2, "male"] = 0

    # We don't keep individuals with unknown sex.
    return sex.set_index("sample_id", verify_integrity=True)["male"]


def main():
    args = parse_args()

    # Parse the configuration.
    configuration = Configuration.from_file(args.configuration)

    if args.sex_path:
        # We will infer the sex of samples. Useful for subsetting or for
        # sex-based exclusions.
        configuration.set_sample_sex(
            detect_sample_sex(configuration, args.sex_path)
        )

        configuration.sex_path = args.sex_path

    if args.male_only or args.female_only:
        configuration.sex_stratified = True
        configuration._analyzed_sex = "MALE" if args.male_only else "FEMALE"

        if not configuration.sample_sex_known():
            raise RuntimeError("Requested sex filtering, but no --sex-path "
                               "provided. Or sex incorrectly detected.")

        if args.male_only:
            configuration.subset = configuration.get_males().tolist()

        elif args.female_only:
            configuration.subset = configuration.get_females().tolist()

        n = len(configuration.subset)
        sex = "male" if args.male_only else "female"

        print(f"Keeping only {sex} individuals (n={n})")

    # We set the number of workers based on the number of CPUs available.
    # Basically, we keep one for the data generator.
    # If this causes heavy loads, we could consider leaving another one for
    # the dealer (currently not accounted for).
    n_workers = max(1, args.cpus - 1)

    # We select the appropriate data generator for the analysis type.
    # If none is provided, we take all of them.
    #
    # ANALYSIS_TO_GENERATOR is defined in the generators module and maps
    # analysis types to tuples of (variable_type, data_generator)
    # variable_type is either "linear" or "binary"
    data_generators = []
    if args.analysis_type is None:
        data_generators.extend(ANALYSIS_TO_GENERATOR.values())

    else:
        generators = [
            ANALYSIS_TO_GENERATOR[analysis_type]
            for analysis_type in args.analysis_type
        ]
        data_generators.extend(generators)

    # We run different analyses for linear and binary phenotypes.
    phewas_kwargs = {
        "n_workers": n_workers,
    }

    # We select and initialize the data generators for both the linear
    # and binary traits.
    if configuration.should_do_linear():
        print("Running analyses for continuous phenotypes.")

        phewas_kwargs["data_generators"] = [
            f(configuration) for t, f in data_generators if t == "linear"
        ]
        phewas_kwargs["worker_script"] = configuration.linear_conf.worker_script
        phewas_kwargs["output_prefix"] = "{}_continuous".format(args.output)

        _run_phewas_in_sandbox(configuration, **phewas_kwargs)

    else:
        print("Skipping analyses for continuous phenotypes.")

    if configuration.should_do_binary():
        print("Running analyses for binary phenotypes.")

        phewas_kwargs["data_generators"] = [
            f(configuration) for t, f in data_generators if t == "binary"
        ]
        phewas_kwargs["worker_script"] = configuration.binary_conf.worker_script
        phewas_kwargs["output_prefix"] = "{}_binary".format(args.output)

        _run_phewas_in_sandbox(configuration, **phewas_kwargs)

    else:
        print("Skipping analyses for binary phenotypes.")


def _run_phewas_in_sandbox(configuration, **kwargs):
    # If the configuration asks for an analysis but the CLI specified analysis
    # types don't include any generators, then this list will be empty and
    # we can skip the analysis.
    if not kwargs["data_generators"]:
        return

    analysis_id = setup_analysis(configuration)
    os.chdir(analysis_id)

    try:
        # Run the PheWAS using the created temporary directory as the working
        # directory.
        run_phewas(configuration, **kwargs)

        # Collect results and remove temporary files.
        collect_and_teardown_analysis(kwargs["output_prefix"])
    finally:
        os.chdir("..")
        shutil.rmtree(analysis_id)
