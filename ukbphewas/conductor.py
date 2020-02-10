import os
import shutil
import time
import itertools
import uuid
import subprocess
import json
import argparse
import multiprocessing
import uuid

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather
import zmq


from .configuration import Configuration
from . import bin
from . import Rpkg
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
        "-e", "source('{}')".format(os.path.join(worker_base, "worker.R")),
        "-e", "source('{}')".format(worker_script),
        str(worker_id), out_addr, monitor_addr
    ]

    proc = subprocess.Popen(cmd)

    return proc, {}


def send_data(socket, metadata, data):
    table = pa.Table.from_pandas(data, preserve_index=False)
    sink = pa.BufferOutputStream()
    stream = pa.RecordBatchStreamWriter(sink, table.schema)
    stream.write_table(table)
    buf = sink.getvalue()

    socket.send_multipart([
        json.dumps(metadata).encode("utf-8"),
        buf.to_pybytes()
    ])


def monitor(ctx, state):
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


def collect_results(output_filename):
    # Infer the analysis id from current directory.
    cwd = os.getcwd()
    analysis_id = os.path.split(cwd)[-1]

    if os.path.isfile("header.csv"):
        with open("header.csv", "r") as f:
            header = f.read()

    else:
        header = None

    with open(os.path.join("..", output_filename), "w") as out:
        if header:
            out.write(header)

        for filename in os.listdir():
            # We only concatenate files that start with results_worker as a
            # convention.
            if not filename.startswith("results_worker"):
                continue

            with open(filename, "r") as f:
                out.write(f.read())


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


def run_phewas(configuration, n_workers, data_generators, output_filename,
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
            args=(context, monitor_state),
            daemon=True
        )
        processes.append(monitor_t)
        monitor_t.start()

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

        for metadata, data in itertools.chain(*data_generators):
            send_data(socket, metadata, data)

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

        # Collect results.
        collect_results(output_filename)


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
        "--sex-column",
        help="Label of the column with the sex variable "
             "(used for subsetting).",
        default=None
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--male-only", action="store_true")
    group.add_argument("--female-only", action="store_true")

    return parser.parse_args()


def create_sex_subset(configuration, sex_column, keep_male):
    # Look for sex_column in the covariables.
    sex = None

    for filename in configuration.covars_filenames:
        if filename.endswith(".csv") or filename.endswith(".csv.gz"):
            df = pd.read_csv(filename)

        elif filename.endswith(".feather"):
            df = pyarrow.feather.read_feather(filename)

        else:
            raise ValueError("Unknown format for the covariables.")

        if sex_column in df.columns:
            sex = df[["sample_id", sex_column]]
            break

    if sex is None:
        raise ValueError(
            f"Could not find column '{sex_column}' in covariables."
        )

    values = sex[sex_column].unique()
    values = set(values[~np.isnan(values)])

    if set(values) == {0, 1}:
        # Sex is encoded as a boolean vector. Try to infer meaning from
        # column name.
        if sex_column == "male":
            # Do Nothing
            pass

        elif sex_column == "female":
            sex["male"] = 1 - sex["female"]

        else:
            raise ValueError(
                "Can't infer sex meaning from column '{}'. Name the sex "
                "column 'male', 'female' or use male = 1, female = 2, "
                "missing = 0."
                "".format(sex_column)
            )

    elif set(values) == {0, 1, 2}:
        print("Assuming 1=male, 2=female, 0=missing for column '{}'."
              "".format(sex_column))

        sex["male"] = np.nan
        sex.loc[sex[sex_column] == 1, "male"] = 1
        sex.loc[sex[sex_column] == 2, "male"] = 0

    # We don't keep individuals with unknown sex.
    sex = sex.dropna(subset=["male"])

    # Now we return the samples to keep.
    if keep_male:
        return sex.loc[sex.male == 1, "sample_id"].unique().tolist()

    else:
        return sex.loc[sex.male == 0, "sample_id"].unique().tolist()


def main():
    args = parse_args()

    # Parse the configuration.
    configuration = Configuration.from_file(args.configuration)

    if args.male_only or args.female_only:
        if args.sex_column is None:
            raise RuntimeError("Requested sex filtering, but not --sex-column "
                               "provided.")

        configuration.subset = create_sex_subset(
            configuration,
            args.sex_column,
            keep_male=args.male_only
        )

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
        phewas_kwargs["output_filename"] = "{}_continuous.csv".format(
            args.output
        )

        _run_phewas_in_sandbox(configuration, **phewas_kwargs)

    else:
        print("Skipping analyses for continuous phenotypes.")

    if configuration.should_do_binary():
        print("Running analyses for binary phenotypes.")

        phewas_kwargs["data_generators"] = [
            f(configuration) for t, f in data_generators if t == "binary"
        ]
        phewas_kwargs["worker_script"] = configuration.binary_conf.worker_script
        phewas_kwargs["output_filename"] = "{}_binary.csv".format(args.output)

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
        run_phewas(configuration, **kwargs)
    finally:
        os.chdir("..")
        shutil.rmtree(analysis_id)
