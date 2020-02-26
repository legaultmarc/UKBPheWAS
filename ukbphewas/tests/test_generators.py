import pandas as pd
import numpy as np

from ..generators import *
from ..configuration import Configuration


DummyConfiguration = lambda: Configuration(
    covars_filenames = [],
    model_rhs = None,
    db_password = None,
    limit = None,
    raw_data_cache = ""
)


class DummyLogistic(object):
    def __init__(self):
        self.min_num_cases = 1


def test_icd10_3chars_cancer_exclusion():
    conf = DummyConfiguration()
    conf.binary_conf = DummyLogistic()
    conf._cache = {}

    conf._cache["diseases"] = pd.DataFrame({
        "eid": ["s1", "s2", "s3"],
        "diag_icd10": ["C502", "C64", "I21"]
    })

    conf._cache["cancer_excl_from_controls"] = np.array(["s1", "s2"])

    expected = {
        "C50": pd.DataFrame({
            "eid": ["s1", "s2"],
            "y": [1, np.nan]
        }),
        "C64": pd.DataFrame({
            "eid": ["s1", "s2"],
            "y": [np.nan, 1]
        }),
        "I21": pd.DataFrame({
            "eid": ["s3"],
            "y": [1]
        }),
    }

    for meta, data in data_generator_icd10_3chars(conf):
        assert meta["variable_id"] in expected
        expct = expected[meta["variable_id"]]

        data = data.sort_values("eid").reset_index(drop=True)

        pd.testing.assert_frame_equal(expct, data)


def test_icd10_block_data_generator():
    conf = DummyConfiguration()
    conf.binary_conf = DummyLogistic()

    conf._cache = {"diseases": pd.DataFrame({
        "eid": ["s1", "s2", "s3", "s4", "s5"],
        "diag_icd10": ["I20", "I25.9", "I513", "A031", "I23.4"]
    })}

    expected = {
        "A00-A09": pd.DataFrame({
            "eid": ["s4"],
            "y": [1]
        }),
        "I20-I25": pd.DataFrame({
            "eid": ["s1", "s2", "s5"],
            "y": [1, 1, 1]
        }),
        "I30-I52": pd.DataFrame({
            "eid": ["s3"],
            "y": [1]
        }),
    }

    for meta, data in data_generator_icd10_block(conf):
        assert meta["variable_id"] in expected
        expct = expected[meta["variable_id"]]

        pd.testing.assert_frame_equal(expct, data)


def test_phecodes_sex_exclusion():
    conf = DummyConfiguration()

    # 187.1 is male only (C60, C63, ...)
    # 174.11 is female only (C50, ...)

    sex = pd.Series(
        [1, 1, np.nan, 0, 0, 1],
        index=[f"s{i+1}" for i in range(6)]
    )
    conf.set_sample_sex(sex)

    conf._cache = {}
    conf._cache["diseases"] = pd.DataFrame({
        "eid": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "diag_icd10": ["C60", "C63", "I21", "C50", "I25", "I200"],
    })

    expected = {
        187.1: pd.DataFrame({
            "eid": ["s1", "s2", "s4", "s5"],
            "y": [1, 1, np.nan, np.nan]
        }),
        # 411.2 excludes Phecodes 410-414.99 (so s5 and s6)
        411.2: pd.DataFrame({
            "eid": ["s3", "s5", "s6"],
            "y": [1, np.nan, np.nan]
        }),
        174.11: pd.DataFrame({
            "eid": ["s1", "s2", "s4", "s6"],
            "y": [np.nan, np.nan, 1, np.nan]
        }),
        # 411.8 excludes Phecodes 410-414.99 (so s3 and s6)
        411.8: pd.DataFrame({
            "eid": ["s3", "s5", "s6"],
            "y": [np.nan, 1, np.nan]
        }),
        # 411.1 excludes Phecodes 410-414.99 (so s3 and s5)
        411.1: pd.DataFrame({
            "eid": ["s3", "s5", "s6"],
            "y": [np.nan, np.nan, 1]
        })
    }

    for i, (meta, data) in enumerate(data_generator_phecodes(conf)):
        assert meta["variable_id"] in expected

        expct = expected[meta["variable_id"]]
        data = data.sort_values("eid").reset_index(drop=True)

        pd.testing.assert_frame_equal(expct, data)

    assert i == len(expected) - 1


def test_phecodes_data_generator():
    conf = DummyConfiguration()

    # We'll simulate:
    # - Some cases of 411.2
    # + I23.4 is not in that definition, but it's in the exclusion range.
    # - A an exclusion from controls

    conf._cache = {"diseases": pd.DataFrame({
        "eid": ["s1", "s2", "s3", "s4", "s5"],
        "diag_icd10": ["I21", "I23.8", "I513", "A031", "I23.4"]
    })}

    expected_status = {
        411.2: [1, 1, 1, 0, np.nan],
        8.5: [0, 0, 0, 1, 0],
        414: [np.nan, np.nan, np.nan, 1]
    }

    for meta, data in data_generator_phecodes(conf):
        expected = expected_status[meta["variable_id"]]

        # Data generators never yield controls so we can drop them.
        expected = np.array([i for i in expected if i != 0])

        data = data.sort_values("eid")

        np.testing.assert_array_equal(data.y.values, expected)


def test_continuous_data_generator():
    conf = DummyConfiguration()
    conf._cache = {"continuous": pd.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s1", "s2"],
        "variable": ["v1", "v1", "v1", "v2", "v2"],
        "value": list(range(5))
    })}

    expected = [
        {
            "variable_id": "v1",
            "data": pd.DataFrame({
                "sample_id": ["s1", "s2", "s3"],
                "v1": [0.0, 1.0, 2.0]
            })
        },
        {
            "variable_id": "v2",
            "data": pd.DataFrame({
                "sample_id": ["s1", "s2"],
                "v2": [3.0, 4.0]
            })
        },
    ]

    generator = data_generator_continuous_variables(conf, _normalize=False)
    for metadata, data in generator:
        expct = expected.pop(0)

        assert "CONTINUOUS_VARIABLE" == metadata["analysis_type"]
        assert expct["variable_id"] == metadata["variable_id"]
        pd.testing.assert_frame_equal(expct["data"], data)


def test_exclusion_range_parser():
    expected = [
        ((1, 9.99), ),
        ((330, 337.99), (341, 349.99)),
    ]

    assert phecode_parse_exclusion_range("001-009.99") == expected.pop(0)
    assert (
        phecode_parse_exclusion_range("330-337.99, 341-349.99") ==
        expected.pop(0)
    )


def test_phecode_in_exclusion_range():
    ranges = phecode_parse_exclusion_range("330-337.99, 341-349.99")

    assert phecode_in_exclusion_range("001", ranges) is False
    assert phecode_in_exclusion_range("329", ranges) is False
    assert phecode_in_exclusion_range("330", ranges) is True
    assert phecode_in_exclusion_range("336", ranges) is True
    assert phecode_in_exclusion_range("337.99", ranges) is True
    assert phecode_in_exclusion_range("338", ranges) is False
    assert phecode_in_exclusion_range("341", ranges) is True
    assert phecode_in_exclusion_range("345", ranges) is True
    assert phecode_in_exclusion_range("349.99", ranges) is True
    assert phecode_in_exclusion_range("350", ranges) is False
