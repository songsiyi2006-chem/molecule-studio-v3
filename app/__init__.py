"""Local molecular-property platform, version 2."""

import os

# Set conservative defaults before scientific libraries are imported. Explicit
# user/runtime settings still take precedence.
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
