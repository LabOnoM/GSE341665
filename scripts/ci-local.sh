#!/usr/bin/env bash
# Local code checks: the two jobs of the former GitHub Actions workflow .github/workflows/checks.yml, removed on
# 2026-10-13 (LabOnoM does not use GitHub Actions). Run before every merge: bash scripts/ci-local.sh
# Needs Python 3.12 (set PYTHON, default python3.12) and R 4.3.2 (Rscript), and network access on the first run.
# The Python environment and the isolated R program library live outside the repository, in
# ${GSE341665_CI_CACHE:-~/.cache/gse341665-ci}, and are reused by later runs. No restricted study input is read.
set -euo pipefail
cd "$(dirname "$0")/.."

cache="${GSE341665_CI_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/gse341665-ci}"
log=$(mktemp)
done_steps=()
step() { printf '\n== %s\n' "$1"; current="$1"; }
pass() { done_steps+=("$current"); }
summary() {
  status=$?
  printf '\n== Summary\n'
  for s in "${done_steps[@]}"; do echo "PASS  $s"; done
  grep -E '^(Ran [0-9]+ tests?|OK|FAILED)' "$log" | sed 's/^/      unittest: /' || true
  if [ "$status" -ne 0 ]; then echo "FAIL  ${current:-setup} (exit $status)"; else echo "All ${#done_steps[@]} steps passed."; fi
  rm -f "$log"
}
trap summary EXIT

# --- Job "regression" (Python 3.12, environment/python-lock.txt) ---
py="${PYTHON:-python3.12}"
step "regression: Python 3.12 environment from environment/python-lock.txt"
"$py" -c 'import sys; assert sys.version_info[:2] == (3, 12), sys.version'
[ -x "$cache/venv/bin/python" ] || "$py" -m venv "$cache/venv"
"$cache/venv/bin/python" -m pip install --quiet -r environment/python-lock.txt
pass

step "regression: python -m unittest discover -s tests -v"
"$cache/venv/bin/python" -m unittest discover -s tests -v 2>&1 | tee "$log"
pass

step "regression: python -m compileall"
"$cache/venv/bin/python" -m compileall -q reproduce.py validate.py verify_data_package.py pipeline environment tests
pass

step "regression: bash -n pipeline/04_download_scref.sh"
bash -n pipeline/04_download_scref.sh
pass

# --- Job "program-model" (R 4.3.2, isolated library from environment/program-packages.tsv) ---
export PROGRAM_R_LIBRARY="$cache/program-library"
# Keep the personal R library (~/R/...) off the search path, as on the clean runner: a different Rcpp loaded
# from it makes setup_program.R report a version mismatch.
export R_LIBS_USER="$PROGRAM_R_LIBRARY"
step "program-model: Rscript environment/setup_program.R --install"
Rscript environment/setup_program.R --library "$PROGRAM_R_LIBRARY" --install
pass

step "program-model: Rscript pipeline/check_program_environment.R"
Rscript pipeline/check_program_environment.R
pass

step "program-model: Rscript tests/test_program_models.R"
Rscript tests/test_program_models.R
pass
