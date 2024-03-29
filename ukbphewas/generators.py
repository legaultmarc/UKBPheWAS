import functools
import collections
import os
from typing import Tuple, Sequence

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


class MissingDataSourceError(Exception):
    def __init__(self, data_source_name):
        self.data_source_name = data_source_name
        super().__init__("Could not initialize data generator because of "
                        f"missing data source '{data_source_name}'.")


def get_or_raise(conf, ds_name):
    cache = conf.get_cache()
    try:
        return cache[ds_name]
    except KeyError:
        raise MissingDataSourceError(ds_name)


@analysis_type("CONTINUOUS_VARIABLE")
def data_generator_continuous_variables(configuration, _normalize=True):

    raw_data = get_or_raise(configuration, "continuous_variables")

    raw_data["variable"] = "cont_" + raw_data["variable"].astype("string")

    if configuration.only_do is not None:
        raw_data = raw_data.loc[
            raw_data.variable.isin(configuration.only_do),
            :
        ]

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

    n_generated = 0

    columns = raw_data.columns

    print("Py: Data generator (continuous) warmed up.")
    for col in columns:
        if col == "sample_id":
            continue

        if configuration.limit and n_generated >= configuration.limit:
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

        n_generated += 1


@analysis_type("SELF_REPORTED")
def data_generator_self_reported(configuration):
    # sample_id, disease_code, disease
    # which represent sample_id, coding, description
    data = get_or_raise(configuration, "self_reported_diseases")
    data["y"] = 1.0

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.sample_id.isin(configuration.subset), :].copy()

    # Contains coding, meaning, selectable, node, parent
    sr_coding = pd.read_csv(
        os.path.join(DATA_ROOT, "coding6_self_reported.tsv.gz"), sep="\t",
        dtype={
            "coding": str,
        }
    )

    tree = Node.from_adjacency_list([
        # We use Node's _data attribute to store coding.
        (
            row.node_id,
            row.parent_id if row.parent_id != 0 else None,
            row.meaning,
            row.coding
        )
        for _, row in sr_coding.iterrows()
    ])

    if configuration.only_do is not None:
        sr_coding = sr_coding.loc[
            sr_coding.node_id.isin(configuration.only_do),
            :
        ]

    n_generated = 0
    for _, coding in sr_coding.iterrows():
        if coding.coding == "99999":
            # Skip unclassifiable.
            continue

        if configuration.limit and n_generated >= configuration.limit:
            break

        # Find all children and code as cases.
        n = [node for _, node in tree.iter_depth_first()
             if node.code == coding.node_id][0]

        children_codings = [node._data for i, node in n.iter_depth_first()]

        all_relevant_codings = children_codings + [coding.coding, ]

        meta = {
            "variable_id": coding.node_id,
            "analysis_type": "SELF_REPORTED",
        }

        cur = data.loc[
            data.disease_code.isin(all_relevant_codings),
            ["sample_id", "y"]
        ]

        if (cur.y == 1).sum() < configuration.binary_conf.min_num_cases:
            continue

        # Exclude individuals with a parent disease from controls.
        parent_codings = [node._data for node in n.iter_ancestors()]
        parent_cases = pd.Index(data.loc[
            data.disease_code.isin(parent_codings),
            "sample_id"
        ])
        to_exclude = parent_cases.difference(cur.sample_id)
        exclusion_frame = pd.DataFrame({"sample_id": to_exclude})
        exclusion_frame["y"] = np.nan

        cur = pd.concat((cur, exclusion_frame))

        yield (meta, cur.reset_index(drop=True))

        n_generated += 1


@analysis_type("CV_ENDPOINTS")
def data_generator_cv_endpoints(configuration):
    data = get_or_raise(configuration, "cardiovascular_outcomes")

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.sample_id.isin(configuration.subset), :].copy()

    if configuration.only_do is not None:
        data = data.loc[data.variable.isin(configuration.only_do), :]

    phenotypes = data.variable.drop_duplicates()

    n_generated = 0
    for phenotype in phenotypes:
        if configuration.limit and n_generated >= configuration.limit:
            break

        cur = data.loc[data["variable"] == phenotype, ["sample_id", "case"]]
        cur.columns = ["sample_id", "y"]

        metadata = {
            "variable_id": phenotype,
            "analysis_type": "CV_ENDPOINTS"
        }

        yield (metadata, cur)

        n_generated += 1


@analysis_type("ICD10_3CHAR")
def data_generator_icd10_3chars(configuration):
    data = get_or_raise(configuration, "icd10")

    # We use the cached exclusions because they should include cancers from
    # ICD9 codes.
    cancer_excl_from_controls = pd.Index(
        get_or_raise(configuration, "cancer_excl_from_controls")["sample_id"]
    )

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.sample_id.isin(configuration.subset), :].copy()

    data.icd10 = data.icd10.str[:3]

    data["y"] = 1

    codes = data[["icd10"]].dropna().drop_duplicates()

    if configuration.only_do is not None:
        codes = codes.loc[codes.icd10.isin(configuration.only_do), :]

    # Add a column for cancer codes.
    codes["chapter"] = codes.icd10.str.get(0)
    codes["num"] = codes.icd10.str.slice(1, 3).astype(int)
    codes["is_cancer"] = ((codes["chapter"] == "C") |
                          ((codes["chapter"] == "D") & (codes["num"] <= 48)))\
                         .astype(int)

    n_generated = 0
    for _, row in codes.iterrows():
        if configuration.limit and n_generated >= configuration.limit:
            break

        metadata = {
            "variable_id": row.icd10,
            "analysis_type": "ICD10_3CHAR",
        }

        cur = data.loc[data.icd10 == row.icd10, ["sample_id", "y"]]

        if cur.shape[0] < configuration.binary_conf.min_num_cases:
            continue

        if row.is_cancer:
            # Exclude other cancer cases from controls.
            to_exclude = pd.DataFrame({
                "sample_id": cancer_excl_from_controls.difference(cur.sample_id)
            })
            to_exclude["y"] = np.nan

            cur = pd.concat((cur, to_exclude))

        yield (metadata, cur)
        n_generated += 1


@analysis_type("ICD10_BLOCK")
def data_generator_icd10_block(configuration):
    data = get_or_raise(configuration, "icd10")

    # Read the ICD10 block metadata.
    icd10_blocks = pd.read_csv(os.path.join(DATA_ROOT, "icd10_blocks.csv.gz"))

    if configuration.only_do is not None:
        icd10_blocks = icd10_blocks.loc[
            icd10_blocks.block.isin(configuration.only_do),
            :
        ]
        if icd10_blocks.shape[0] == 0:
            return

    lr = icd10_blocks.block.str.split("-", expand=True)
    lr.columns = ["left", "right"]

    icd10_blocks = pd.concat((icd10_blocks, lr), axis=1)

    icd10_blocks["left_chapter"] = icd10_blocks.left.str.get(0)
    icd10_blocks["right_chapter"] = icd10_blocks.right.str.get(0)

    icd10_blocks["left_num"] = icd10_blocks.left.str.slice(1, 3).astype(int)
    icd10_blocks["right_num"] = icd10_blocks.right.str.slice(1, 3).astype(int)

    # Apply subset.
    if configuration.subset:
        data = data.loc[data.sample_id.isin(configuration.subset), :].copy()

    if data.icd10.isnull().sum() >= 1:
        data = data.dropna(subset=["icd10"]).copy()

    data["y"] = 1

    data["chapter"] = data.icd10.str.get(0)
    data["num"] = data.icd10.str.slice(1, 3).astype(int)

    n_generated = 0
    for _, row in icd10_blocks.iterrows():
        if configuration.limit and n_generated >= configuration.limit:
            break

        # Find all individuals with a code in the block.
        cur = data.loc[
            (data.icd10.str.get(0) == row.left_chapter) & # Match chapter
            ((row.left_num <= data.num) &  # Match number
             (data.num <= row.right_num)), :
        ]

        if cur.shape[0] < configuration.binary_conf.min_num_cases:
            continue

        metadata = {
            "variable_id": row.block,
            "analysis_type": "ICD10_BLOCK",
        }

        yield metadata, cur[["sample_id", "y"]].reset_index(drop=True)

        n_generated += 1


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


SexExclusion = collections.namedtuple(
    "SexExclusion",
    ("male_only", "female_only")
)


def _build_phecode_index(data):
    data["first_char"] = data.phecode.str.get(0)
    out = {}
    for char in data["first_char"].drop_duplicates():
        out[char] = data.loc[data["first_char"] == char, :]

    return out


@analysis_type("PHECODES")
def data_generator_phecodes(configuration):
    data = get_or_raise(configuration, "icd10")

    # Strip dots from ICD10 codes.
    data.icd10 = data.icd10.str.replace(".", "", regex=False)

    # Apply subset.
    # This includes sex-stratified analyses.
    if configuration.subset:
        data = data.loc[data.sample_id.isin(configuration.subset), :]

    # Get an index of all samples with data.
    all_samples = pd.Index(data.sample_id)

    # Read the phecode map for ICD10.
    icd10_to_phecode = pd.read_csv(
        os.path.join(DATA_ROOT, "Phecode_map_v1_2_icd10_beta.csv.gz"),
        dtype={"PHECODE": str}
    )
    icd10_to_phecode = icd10_to_phecode.iloc[:, :3]
    icd10_to_phecode.columns = ["icd10", "phecode", "exclude_phecode"]

    icd10_to_phecode["icd10"] = icd10_to_phecode["icd10"]\
        .str.replace(".", "", regex=False)

    # Join with the diseases.
    data = pd.merge(data, icd10_to_phecode,
                    left_on="icd10", right_on="icd10")

    data = data.dropna(subset=["phecode"])
    data = data.drop_duplicates(subset=["sample_id", "phecode"])

    # Read the gender exclusions.
    gender_exclusions = pd.read_csv(
        os.path.join(DATA_ROOT, "gender_restriction.csv.gz"),
        dtype={"phecode": str}
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

    # Preload the males and females if available to accelerate sex-based
    # exclusions.
    males = males_frame = analysis_male_only = None
    females = females_frame = analysis_female_only = None
    if configuration.sample_sex_known():
        males = configuration.get_males()
        males_frame = pd.DataFrame({"sample_id": males, "y": np.nan})

        females = configuration.get_females()
        females_frame = pd.DataFrame({"sample_id": females, "y": np.nan})

        data["male"] = np.nan
        data.loc[data["sample_id"].isin(males), "male"] = 1
        data.loc[data["sample_id"].isin(females), "male"] = 0

        # Check if analysis is already male or female only.
        # If this is the case, we'll skip the exclusions.
        if configuration.get_analyzed_sex() == "BOTH":
            analysis_male_only = analysis_female_only = False
        else:
            analysis_male_only = configuration.get_analyzed_sex() == "MALE"
            analysis_female_only = configuration.get_analyzed_sex() == "FEMALE"

    # Build a dict of 1st phecode char to subset of the data to enhance
    # performance.
    phecode_index = _build_phecode_index(data)

    codes = data.phecode.dropna().drop_duplicates()

    if configuration.only_do is not None:
        only_do = set(configuration.only_do)
        codes = [code for code in codes if code in only_do]

    n_generated = 0
    for code in codes:
        if configuration.limit and n_generated >= configuration.limit:
            break

        metadata = {
            "variable_id": code,
            "analysis_type": "PHECODES",
        }

        # Dict based on first character to accelerate indexing.
        cur = phecode_index.get(code[0])

        # Check gender exclusion
        if configuration.sample_sex_known():
            if code in gender_exclusions.index:
                sex_excl = gender_exclusions.loc[code, :]

            else:
                # If it's not in the metadata file, we assume that both sexes
                # should be analyzed.
                sex_excl = SexExclusion(False, False)

            if sex_excl.male_only:
                if analysis_female_only:
                    continue

                metadata["sex_subset"] = "MALE_ONLY"
                cur = cur.loc[cur.male == 1, :]

            elif sex_excl.female_only:
                if analysis_male_only:
                    continue

                metadata["sex_subset"] = "FEMALE_ONLY"
                cur = cur.loc[cur.male == 0, :]

            else:
                metadata["sex_subset"] = "BOTH"

        else:
            metadata["sex_subset"] = "SEX_UNKNOWN"

        # If code == 278, we should also include 278.*
        # If it's 278.1 we shouldn't...
        #
        # Same goes if code is 008.5 it whould include 008.5*
        # If it's 008.52 it shouldn't include parents.
        cur = cur.loc[
            cur.phecode.str.slice(stop=len(code)) == code,
            ["sample_id", "y"]
        ]

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
            # Unfortunately, codes to exclude don't always have the same 1st
            # character, so we can't use the smaller dataset.
            # This step is not that slow.
            samples_to_exclude = pd.Index(data.loc[
                data.phecode.isin(codes_to_exclude),
                "sample_id"
            ].drop_duplicates())

            # We take the samples that are not in the cases to exclude
            # (exclude from the controls).
            samples_to_exclude = samples_to_exclude.difference(cur.sample_id)
            samples_to_exclude = pd.DataFrame(
                {"sample_id": samples_to_exclude, "y": np.nan}
            )

            cur = pd.concat((cur, samples_to_exclude), axis=0)

        if (cur.y == 1).sum() < configuration.binary_conf.min_num_cases:
            continue

        # Apply sex exclusion.
        if metadata["sex_subset"] == "MALE_ONLY":
            cur = pd.concat((cur, females_frame))
            cur = cur.drop_duplicates()

        elif metadata["sex_subset"] == "FEMALE_ONLY":
            cur = pd.concat((cur, males_frame))
            cur = cur.drop_duplicates()

        yield (metadata, cur)

        n_generated += 1


class Node(object):
    def __init__(self):
        self.is_root = False

        self.parent = None
        self.children = []

        self.code = None
        self.description = None
        self._data = None


    def __repr__(self):
        parent_code = self.parent.code if self.parent else None

        if self.is_root:
            return "<Root Node>"

        return (
            "<Node '{}' - `{}` [parent is '{}' | {} children]>"
            "".format(self.code, self.description, parent_code,
                      len(self.children))
        )

    def iter_depth_first(self, level=0):
        """Depths first tree traversal rooted at this node."""
        if level > 0:
            yield level, self

        for child in self.children:
            yield from child.iter_depth_first(level + 1)


    def iter_ancestors(self):
        """Returns the chain of ancestors up to the root."""
        if self.parent is None:
            return

        else:
            yield self.parent
            yield from self.parent.iter_ancestors()

    @classmethod
    def from_adjacency_list(cls, li):
        """Build the tree from an adjacency list.

        The form of the adjacency list is:
        (code, parent, description=None, data=None)

        """
        root = cls()
        root.is_root = True

        node_dict = {}

        # First create all nodes.
        for code, parent, description, data in li:
            n = cls()

            n.code = code
            n.description = description
            n._data = data
            node_dict[code] = n

        # Set relationships.
        for code, parent, description, data in li:
            n = node_dict[code]

            if parent is None:
                n.parent = root
                root.children.append(n)

            else:
                n.parent = node_dict[parent]
                n.parent.children.append(n)

        return root
