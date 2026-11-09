#!/bin/bash
# Build the three auxiliary Python environments for the revision round.
# Plain pip venvs, NOT conda: `conda create --clone` deadlocks on this cluster,
# and slr-env must never be mutated (it backs every computed result).
set -euo pipefail
BASE=/home/itg/oleg.vlasovets/slr_example/ige-gut-microbiome-analysis/revision/env
PY=/home/itg/oleg.vlasovets/miniconda3/envs/slr-env/bin/python

# 1. xgb-venv: inherits slr-env's exact numpy/sklearn/pandas via
#    --system-site-packages, adds only xgboost. No numeric drift vs computed results.
if [ ! -x "$BASE/xgb-venv/bin/python" ]; then
  "$PY" -m venv --system-site-packages "$BASE/xgb-venv"
  "$BASE/xgb-venv/bin/pip" install --quiet "xgboost==2.1.*"
fi
"$BASE/xgb-venv/bin/python" -c "import xgboost, sklearn, numpy, pandas; \
print('xgb-venv OK: xgboost', xgboost.__version__, '| sklearn', sklearn.__version__, \
'| numpy', numpy.__version__, '| pandas', pandas.__version__)"
"$BASE/xgb-venv/bin/pip" freeze | grep -i xgboost > "$BASE/xgb-requirements.txt"

# 2. legacy-venv: ISOLATED (no system site packages) - the point is pandas 2.2.3,
#    which still has DataFrame.applymap and the era-faithful semantics of the
#    published notebook. numpy pinned to the 1.26 line of the original stack.
if [ ! -x "$BASE/legacy-venv/bin/python" ]; then
  "$PY" -m venv "$BASE/legacy-venv"
  "$BASE/legacy-venv/bin/pip" install --quiet "pandas==2.2.3" "numpy==1.26.4" \
      "statsmodels==0.14.*" "pyyaml"
fi
"$BASE/legacy-venv/bin/python" -c "import pandas, numpy, statsmodels; \
assert hasattr(pandas.DataFrame, 'applymap'), 'applymap missing'; \
print('legacy-venv OK: pandas', pandas.__version__, '| numpy', numpy.__version__, \
'| statsmodels', statsmodels.__version__)"
"$BASE/legacy-venv/bin/pip" freeze > "$BASE/legacy-requirements.txt"
# 3. classo-venv: slr-env's numpy/sklearn/pandas via --system-site-packages plus
#    c-lasso at the commit that includes C. L. Mueller's C1/C2 classification-path
#    fixes (bio-datascience/c-lasso PRs #54, #55, merged 2026-11-04). Pinned to a
#    SHA: the PyPI release (1.0.11, inherited from slr-env) predates them.
#    --no-deps so nothing in the inherited scientific stack is upgraded.
CLASSO_SHA=c1c33b83b3aac166e017f548681b1fce38fb2a47
# (re)create when pip is missing: an interrupted creation leaves python but no pip
[ -x "$BASE/classo-venv/bin/pip" ] || "$PY" -m venv --clear --system-site-packages "$BASE/classo-venv"
# install unless the venv already has its OWN c-lasso (not the inherited PyPI one)
if ! "$BASE/classo-venv/bin/python" -c "import sys, classo; sys.exit('classo-venv' not in classo.__file__)" 2>/dev/null; then
  "$BASE/classo-venv/bin/pip" install --quiet --no-deps \
      "git+https://github.com/bio-datascience/c-lasso@${CLASSO_SHA}"
fi
"$BASE/classo-venv/bin/python" -c "import classo, sklearn, numpy; \
assert 'classo-venv' in classo.__file__, 'picked up the PyPI c-lasso from slr-env'; \
print('classo-venv OK: c-lasso', classo.__version__, '|', classo.__file__)"
echo "c-lasso @ git+https://github.com/bio-datascience/c-lasso@${CLASSO_SHA}" > "$BASE/classo-requirements.txt"
echo "ALL VENVS READY"
