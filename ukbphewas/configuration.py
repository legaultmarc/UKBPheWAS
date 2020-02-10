# These should set the worker script dynamically.
# For now, only linear LRT is implemented.
import os
import json

from . import Rpkg
from .db import load_or_create_data_cache


class _BinaryConfiguration(object):
    pass


class SkipBinary(_BinaryConfiguration):
    pass


class DoBinary(_BinaryConfiguration):
    def __init__(
        self,
        include_death_records=True,
        include_secondary_hospit=True,
        min_num_cases=50,
        skip_icd10_raw=False,
    ):
        self.include_death_records = include_death_records
        self.include_secondary_hospit = include_secondary_hospit
        self.min_num_cases = min_num_cases
        self.skip_icd10_raw = skip_icd10_raw


class DoBinaryLRT(DoBinary):
    def __init__(self, augmented_variables, *args, **kwargs):
        self.worker_script = os.path.abspath(os.path.join(
            os.path.dirname(Rpkg.__file__),
            "R", "logistic_deviance_diff_test_worker.R"
        ))

        self.augmented_variables = augmented_variables
        super().__init__(*args, **kwargs)


class _LinearConfiguration(object):
    pass


class SkipLinear(_LinearConfiguration):
    pass


class DoLinear(_LinearConfiguration):
    def __init__(self):
        raise NotImplementedError()


class DoFTest(_LinearConfiguration):
    def __init__(self, augmented_variables):
        self.worker_script = os.path.abspath(os.path.join(
            os.path.dirname(Rpkg.__file__),
            "R", "linear_f_test_worker.R"
        ))

        self.augmented_variables = augmented_variables


class Configuration(object):
    def __init__(
        self,
        covars_filenames,
        model_rhs,
        db_password=None,
        limit=None,
        continuous_variables_path="/data/projects/uk_biobank/data/reports/SGR-2094/",
        raw_data_cache=None,
        subset=None,

        linear_conf=SkipLinear(),
        binary_conf=SkipBinary(),

    ):

        # We accept many covar filenames which will be (outer) joined if needed
        if isinstance(covars_filenames, str):
            covars_filenames = [covars_filenames]

        self.covars_filenames = [
            os.path.abspath(fn) for fn in covars_filenames
        ]

        self.raw_data_cache = os.path.abspath(raw_data_cache)
        self._cache = None

        self.model_rhs = model_rhs
        self.db_password = db_password
        self.limit = limit
        self.continuous_variables_path = os.path.abspath(
            continuous_variables_path
        )
        self.subset = subset

        self.linear_conf = linear_conf
        self.binary_conf = binary_conf

    def should_do_linear(self):
        return not isinstance(self.linear_conf, SkipLinear)

    def should_do_binary(self):
        return not isinstance(self.binary_conf, SkipBinary)

    @classmethod
    def from_file(cls, filename):
        """Initialize a Configuration from a file."""
        global_vars = {}
        with open(filename) as f:
            conf_script = f.read()

        exec(conf_script, global_vars)

        if "configuration" not in global_vars:
            raise ValueError(
                "The configuration script should define a global variable "
                "named 'configuration'."
            )

        configuration = global_vars["configuration"]
        if not isinstance(configuration, cls):
            raise ValueError(
                "The 'configuration' variable in the configuration script "
                "should be an instance of {} (got {})."
                "".format(cls, type(configuration))
            )

        return configuration

    def to_json(self):
        try:
            d = self.__dict__.copy()

            if d.get("_cache") is not None:
                del d["_cache"]

            d["linear_conf"] = self.linear_conf.__dict__
            d["linear_conf"]["type"] = self.linear_conf.__class__.__name__

            d["binary_conf"] = self.binary_conf.__dict__
            d["binary_conf"]["type"] = self.binary_conf.__class__.__name__

            return json.dumps(d)

        except Exception as e:
            print(self.__dict__)
            print(e)
            raise e

    def get_cache(self):
        if self._cache is not None:
            return self._cache

        self._cache = load_or_create_data_cache(self)
        return self._cache
