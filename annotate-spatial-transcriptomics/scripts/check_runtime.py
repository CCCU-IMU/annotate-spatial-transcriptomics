#!/usr/bin/env python3
"""Report optional Python capabilities without installing or mutating an environment."""
import importlib.util,json,platform,sys
mods=["numpy","pandas","scipy","sklearn","anndata","scanpy","matplotlib","h5py"]
result={"python":sys.version,"platform":platform.platform(),"modules":{m:bool(importlib.util.find_spec(m)) for m in mods}}
print(json.dumps(result,indent=2))
raise SystemExit(0 if all(result["modules"][m] for m in ["numpy","pandas","scipy","sklearn"]) else 2)
