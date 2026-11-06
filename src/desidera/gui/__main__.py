"""Launch DESIDERA Project Hub with the same setup as `desidera-gui`.

Set DESIDERA_WEBENGINE_SOFTWARE=1 on Windows for an opt-in GPU-compositor
diagnostic; this is not a scientifically validated display fix.
"""
from . import main

if __name__ == "__main__":
    raise SystemExit(main())
