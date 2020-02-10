import pandas as pd

def create_data_cache(configuration):
    # This is todo, for now we will simply use the R pre-built raw data cache.
    raise NotImplemented()


def load_or_create_data_cache(configuration):
    cache = configuration.raw_data_cache

    expected_keys = [
        "diseases",
        "cancer_excl_from_controls",
        "cv_endpoints",
        "self_reported_diseases",
        "continuous",
        "continuous_metadata"
    ]

    if cache is None:
        return create_data_cache(configuration)

    if cache.endswith(".h5") or cache.endswith(".hdf5"):
        import h5py

        print(f"Loading data cache from HDF5 file: '{cache}'")
        out = {}
        for k in expected_keys:
            if k == "cancer_excl_from_controls":
                # This is not supported by pandas (it's just an array), so
                # we need to do it manually.
                f = h5py.File(cache, "r")
                out[k] = pd.DataFrame(
                    {"sample_id": f["cancer_excl_from_controls"][:]}
                )

            else:
                out[k] = pd.read_hdf(cache, key=k)

        return out

    raise ValueError(f"Invalid cache file: '{cache}'")
