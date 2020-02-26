import functools
import os
from typing import Tuple, Sequence

import pyarrow.feather
import pandas as pd
import numpy as np

from . import data


DATA_ROOT = os.path.abspath(os.path.dirname(data.__file__))


ANALYSIS_TO_GENERATOR = {}


class analysis_type(object):
    """Decorator to automatically populate the ANALYSIS_TO_GENERATOR."""
    def __init__(self, analysis_type):
        self.analysis_type = analysis_type

        if self.analysis_type == "CONTINUOUS_VARIABLE":
            self.variable_type = "linear"

        else:
            self.variable_type = "binary"

    def __call__(self, data_gen):
        global ANALYSIS_TO_GENERATOR
        ANALYSIS_TO_GENERATOR[self.analysis_type] = (
            self.variable_type, data_gen
        )

        return data_gen


@analysis_type("CONTINUOUS_VARIABLE")
def data_generator_continuous_variables(configuration, _normalize=True):
    raw_data = configuration.get_cache()["continuous"]

    # Pivot into wide format.
    raw_data = raw_data.pivot_table(
        index="sample_id",
        values="value",
        columns="variable"
    )

    raw_data = raw_data.reset_index()

    # Use a regular column instead of pandas Index.
    raw_data.columns = raw_data.columns.to_list()

    assert raw_data.index.duplicated().sum() == 0

    # Apply subset.
    if configuration.subset:
        raw_data = raw_data.loc[
            raw_data.sample_id.isin(configuration.subset), :
        ]

    for i, col in enumerate(raw_data.columns):
        if col == "sample_id":
            continue

        if configuration.limit and i >= configuration.limit:
            break

        metadata = {
            "variable_id": col,
            "analysis_type": "CONTINUOUS_VARIABLE",
        }

        # Standardize so that the betas are on the same scale (s.d.).
        y = raw_data[["sample_id", col]].dropna()
        if _normalize:
            y[col] = (y[col] - y[col].mean()) / y[col].std()

        yield (metadata, y)


@analysis_type("ICD10_3CHAR")
def data_generator_icd10_3chars(configuration):
    cache = configuration.get_cache()
    data = cache["diseases"]

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.eid.isin(configuration.subset), :].copy()

    data.diag_icd10 = data.diag_icd10.str[:3]

    data["y"] = 1

    codes = data[["diag_icd10"]].drop_duplicates()

    # Add a column for cancer codes.
    codes["chapter"] = codes.diag_icd10.str.get(0)
    codes["num"] = codes.diag_icd10.str.slice(1, 3).astype(int)
    codes["is_cancer"] = ((codes["chapter"] == "C") |
                          ((codes["chapter"] == "D") & (codes["num"] <= 48)))\
                         .astype(int)

    # We use the cached exclusions because they should include cancers from
    # ICD9 codes.
    cancer_excl_from_controls = pd.Index(cache["cancer_excl_from_controls"])

    for i, row in codes.iterrows():
        if configuration.limit and i >= configuration.limit:
            break

        metadata = {
            "variable_id": row.diag_icd10,
            "analysis_type": "ICD10_3CHAR",
        }

        cur = data.loc[data.diag_icd10 == row.diag_icd10, ["eid", "y"]]

        if cur.shape[0] < configuration.binary_conf.min_num_cases:
            continue

        if row.is_cancer:
            # Exclude other cancer cases from controls.
            to_exclude = pd.DataFrame({
                "eid": cancer_excl_from_controls.difference(cur.eid)
            })
            to_exclude["y"] = np.nan

            cur = pd.concat((cur, to_exclude))

        yield (metadata, cur)


@analysis_type("ICD10_BLOCK")
def data_generator_icd10_block(configuration):
    data = configuration.get_cache()["diseases"]

    # Read the ICD10 block metadata.
    icd10_blocks = pd.read_csv(os.path.join(DATA_ROOT, "icd10_blocks.csv.gz"))

    lr = icd10_blocks.block.str.split("-", expand=True)
    lr.columns = ["left", "right"]

    icd10_blocks = pd.concat((icd10_blocks, lr), axis=1)

    icd10_blocks["left_chapter"] = icd10_blocks.left.str.get(0)
    icd10_blocks["right_chapter"] = icd10_blocks.right.str.get(0)

    icd10_blocks["left_num"] = icd10_blocks.left.str.slice(1, 3).astype(int)
    icd10_blocks["right_num"] = icd10_blocks.right.str.slice(1, 3).astype(int)

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.eid.isin(configuration.subset), :].copy()

    data["y"] = 1

    data["chapter"] = data.diag_icd10.str.get(0)
    data["num"] = data.diag_icd10.str.slice(1, 3).astype(int)

    for i, row in icd10_blocks.iterrows():
        if configuration.limit and i >= configuration.limit:
            break

        # Find all individuals with a code in the block.
        cur = data.loc[
            (data.diag_icd10.str.get(0) == row.left_chapter) & # Match chapter
            ((row.left_num <= data.num) &  # Match number
             (data.num <= row.right_num)), :
        ]

        if cur.shape[0] < configuration.binary_conf.min_num_cases:
            continue

        metadata = {
            "variable_id": row.block,
            "analysis_type": "ICD10_BLOCK",
        }

        yield metadata, cur[["eid", "y"]].reset_index(drop=True)


# Left and right phecodes of the range defining the exclusion.
ExclusionRange = Tuple[float, float]


@functools.lru_cache(maxsize=64)
def phecode_parse_exclusion_range(s: str) -> Tuple[ExclusionRange, ...]:
    """Parse a phecode exclusion range from a string."""
    excl_ranges = []
    for excl_expr in s.split(","):
        left, right = [float(i) for i in excl_expr.strip().split("-")]
        excl_ranges.append((left, right))

    return tuple(excl_ranges)


@functools.lru_cache(maxsize=64)
def phecode_in_exclusion_range(code: str,
                               ranges: Sequence[ExclusionRange]) -> bool:
    code_num = float(code)

    for left, right in ranges:
        if left <= code_num <= right:
            return True

    return False


@analysis_type("PHECODES")
def data_generator_phecodes(configuration):
    # TODO Apply gender restriction if available.
    data = configuration.get_cache()["diseases"]

    # Strip dots from ICD10 codes.
    data.diag_icd10 = data.diag_icd10.str.replace(".", "", regex=False)

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.eid.isin(configuration.subset), :]

    # Get an index of all samples with data.
    all_samples = pd.Index(data.eid)

    # Read the phecode map for ICD10.
    icd10_to_phecode = pd.read_csv(
        os.path.join(DATA_ROOT, "Phecode_map_v1_2_icd10_beta.csv.gz")
    )
    icd10_to_phecode = icd10_to_phecode.iloc[:, :3]
    icd10_to_phecode.columns = ["icd10", "phecode", "exclude_phecode"]

    icd10_to_phecode["icd10"] = icd10_to_phecode["icd10"]\
        .str.replace(".", "", regex=False)

    # Join with the diseases.
    data = pd.merge(data, icd10_to_phecode,
                    left_on="diag_icd10", right_on="icd10")

    # Read the gender exclusions.
    gender_exclusions = pd.read_csv(
        os.path.join(DATA_ROOT, "gender_restriction.csv.gz")
    )
    gender_exclusions = gender_exclusions.set_index(
        "phecode", verify_integrity=True
    )

    data["y"] = 1

    # We build a map of Phecode to exclusion range.
    # There are only 1,572 unique phecodes so this should not be too big.
    phecode_to_excl_ranges = {}
    for i, row in icd10_to_phecode.dropna(subset=["exclude_phecode"]).iterrows():
        phecode_to_excl_ranges[row.phecode] = phecode_parse_exclusion_range(
            row.exclude_phecode
        )

    # Preload the males and females if available to acclerate sex-based
    # exclusions.
    males = None
    females = None
    if configuration.sample_sex_known():
        males = configuration.get_males()
        females = configuration.get_females()

    codes = data.phecode.dropna().drop_duplicates()
    for i, code in enumerate(codes):
        if configuration.limit and i >= configuration.limit:
            break

        metadata = {
            "variable_id": code,
            "analysis_type": "PHECODES",
        }

        cur = data.loc[data.phecode == code, ["eid", "y"]]

        # Check gender exclusion
        if configuration.sample_sex_known():
            sex_excl = gender_exclusions.loc[code, :]

            if sex_excl.male_only:
                cur = pd.concat((
                    cur,
                    pd.DataFrame({"eid": females, "y": np.nan})
                ))
                metadata["sex_subset"] = "MALE_ONLY"

            elif sex_excl.female_only:
                cur = pd.concat((
                    cur,
                    pd.DataFrame({"eid": males, "y": np.nan})
                ))
                metadata["sex_subset"] = "FEMALE_ONLY"

            else:
                metadata["sex_subset"] = "BOTH"

        else:
            metadata["sex_subset"] = "SEX_UNKNOWN"


        # Define samples to exclude.
        exclusion_range = phecode_to_excl_ranges.get(code)

        if exclusion_range is not None:
            # Find codes to exclude in the metadata.
            codes_to_exclude = icd10_to_phecode.apply(
                lambda row: phecode_in_exclusion_range(
                    row.phecode, exclusion_range
                ),
                axis=1
            )

            codes_to_exclude = icd10_to_phecode.loc[
                codes_to_exclude,
                "phecode"
            ].drop_duplicates()

            # Get samples that have excluded codes.
            samples_to_exclude = pd.Index(data.loc[
                data.phecode.isin(codes_to_exclude),
                "eid"
            ].drop_duplicates())

            # We take the samples that are not in the cases to exclude
            # (exclude from the controls).
            samples_to_exclude = samples_to_exclude.difference(cur.eid)
            samples_to_exclude = pd.DataFrame(
                {"eid": samples_to_exclude, "y": np.nan}
            )

            cur = pd.concat((cur, samples_to_exclude), axis=0)

        yield (metadata, cur)
