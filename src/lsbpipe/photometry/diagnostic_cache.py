"""Project-local diagnostic cache. JSON + numeric NPY, never pickle."""
from __future__ import annotations
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile
import numpy as np

SCHEMA="desidera.diagnostic-cache/1"


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()


def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_zip(path,write):
    path=Path(path)
    if path.parent.is_symlink() or path.is_symlink() or (path.exists() and not path.is_file()):
        raise ValueError("Unsafe diagnostic archive path")
    fd,temp=tempfile.mkstemp(prefix=".diagnostic-",suffix=".tmp",dir=path.parent)
    try:
        with os.fdopen(fd,"w+b") as f:
            with zipfile.ZipFile(f,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
                write(z)
            f.flush(); os.fsync(f.fileno())
        os.replace(temp,path)
    finally:
        if os.path.exists(temp): os.unlink(temp)


def write_cache(path,key,payload):
    def write(z):
        arrays=[]
        def enc(v):
            if isinstance(v,np.generic): return enc(v.item())
            if isinstance(v,np.ndarray):
                if v.dtype.hasobject: raise TypeError("Object array prohibited in diagnostic cache")
                name=f"array-{len(arrays):05d}.npy"
                f=io.BytesIO(); np.save(f,v,allow_pickle=False); b=f.getvalue()
                arrays.append(dict(name=name,sha256=hashlib.sha256(b).hexdigest()))
                z.writestr(name,b)
                return {"__cache_type__":"array","name":name}
            if v is None or isinstance(v,(str,bool,int)): return v
            if isinstance(v,float):
                return v if math.isfinite(v) else {"__cache_type__":"float","value":str(v)}
            if isinstance(v,tuple): return {"__cache_type__":"tuple","value":[enc(x) for x in v]}
            if isinstance(v,list): return [enc(x) for x in v]
            if isinstance(v,dict):
                if all(isinstance(k,str) for k in v) and "__cache_type__" not in v:
                    return {k:enc(x) for k,x in v.items()}
                return {"__cache_type__":"mapping","value":[[enc(k),enc(x)] for k,x in v.items()]}
            raise TypeError(f"Unsupported cache value: {type(v).__name__}")
        data=enc(payload)
        manifest=dict(schema=SCHEMA,key=key,payload=data,arrays=arrays)
        z.writestr("manifest.json",json.dumps(manifest,allow_nan=False,ensure_ascii=False,separators=(",",":")))
    atomic_zip(path,write)


def read_cache(path,key):
    path=Path(path)
    if path.is_symlink(): raise ValueError("Unsafe diagnostic cache symlink")
    if not path.is_file(): return None
    try:
        with zipfile.ZipFile(path) as z:
            m=json.loads(z.read("manifest.json"))
            if m.get("schema")!=SCHEMA or m.get("key")!=key: return None
            arrays={x["name"]:x["sha256"] for x in m["arrays"]}
            if len(arrays)!=len(m["arrays"]) or len(z.namelist())!=len(set(z.namelist())):
                raise ValueError("Duplicate cache entry")
            if set(z.namelist())!={"manifest.json",*arrays}: raise ValueError("Unexpected cache entry")
            def dec(v):
                if isinstance(v,list): return [dec(x) for x in v]
                if not isinstance(v,dict): return v
                tag=v.get("__cache_type__")
                if tag=="array":
                    name=v["name"]; b=z.read(name)
                    if hashlib.sha256(b).hexdigest()!=arrays[name]: raise ValueError("Array hash mismatch")
                    a=np.load(io.BytesIO(b),allow_pickle=False)
                    if not isinstance(a,np.ndarray) or a.dtype.hasobject: raise ValueError("Invalid cache array")
                    a.setflags(write=False)
                    return a
                if tag=="float":
                    if v["value"] not in ("nan","inf","-inf"): raise ValueError("Invalid float tag")
                    return float(v["value"])
                if tag=="tuple": return tuple(dec(x) for x in v["value"])
                if tag=="mapping": return {dec(k):dec(x) for k,x in v["value"]}
                if tag is not None: raise ValueError("Unknown cache tag")
                return {k:dec(x) for k,x in v.items()}
            return dec(m["payload"])
    except (OSError,ValueError,KeyError,TypeError,EOFError,zipfile.BadZipFile,RuntimeError):
        return None


def cached(path,key,compute,*,force=False):
    value=None if force else read_cache(path,key)
    if value is not None: return value,"HIT",file_hash(path)
    value=compute()
    write_cache(path,key,value)
    return value,"COMPUTED",file_hash(path)
