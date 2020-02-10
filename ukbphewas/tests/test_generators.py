import pandas as pd

from ..generators import *
from ..configuration import Configuration


DummyConfiguration = lambda: Configuration(
    covars_filenames = [],
    model_rhs = None,
    db_password = None,
    limit = None,
    raw_data_cache = ""
)


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
