"""
Module for native serialization and deserialization of DataFrames.
"""

import pandas as pd
import numpy as np
import json
import struct
import io

from typing import Any 


def _np_dtype_to_dtype(dtype: np.dtype) -> str:
    if dtype.name.startswith("float") or dtype.name.startswith("int"):
        return dtype.name
    
    # Should be a string.
    assert dtype.str.startswith("|S"), dtype
    return "string"


def _dtype_to_np_dtype(dtype: str) -> np.dtype:
    if dtype.startswith("float") or dtype.startswith("int"):
        return np.dtype(dtype)
    
    assert dtype == "string"
    raise ValueError(dtype)


def serialize(df: pd.DataFrame) -> bytes:
    out_buffer = io.BytesIO()
    
    metadata = []
    data = []
    for col in df.columns:
        cur = df[col].to_numpy()
        if cur.dtype.name == "object":
            # Serialize str
            cur = cur.astype("S")

        cur_bytes = cur.tobytes()
            
        meta = {
            "name": col,
            "n": df.shape[0],
            "n_bytes": len(cur_bytes),
            "dtype": _np_dtype_to_dtype(cur.dtype)
        }
        
        metadata.append(meta)
        data.append(cur_bytes)
    
    metadata_bytes = json.dumps(metadata).encode("ascii")
    
    out_buffer.write(struct.pack("I", len(metadata_bytes)))
    out_buffer.write(metadata_bytes)
    
    for b in data:
        out_buffer.write(b)
    
    return out_buffer.getvalue()


def deserialize_vector(metadata: dict, buffer: Any) -> np.array:
    if metadata["dtype"] == "string":
        # Infer n bytes per sample.
        n_bytes_sample = metadata["n_bytes"] / metadata["n"]
        assert n_bytes_sample.is_integer()
        n_bytes_sample = int(n_bytes_sample)

        np_dtype=f"S{n_bytes_sample}"

    else:
        np_dtype = _dtype_to_np_dtype(metadata["dtype"])
    
    v = np.frombuffer(
        buffer.read(metadata["n_bytes"]),
        dtype=np_dtype,
        count=metadata["n"]
    )
    
    return v if metadata["dtype"] != "string" else v.astype(str)


def deserialize(buffer) -> pd.DataFrame:
    # Get the length of metadata.
    metadata_len = struct.unpack("I", buffer.read(4))[0]
    
    metadata_bytes = buffer.read(metadata_len)
    metadata = json.loads(metadata_bytes)
        
    return pd.DataFrame({
        meta["name"]: deserialize_vector(meta, buffer) for meta in metadata
    })
