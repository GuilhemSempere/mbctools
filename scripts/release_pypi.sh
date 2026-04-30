#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -z "${TWINE_PASSWORD:-}" ]]; then
  echo "TWINE_PASSWORD is not set. Export your PyPI token first." >&2
  echo "Example: export TWINE_USERNAME=__token__" >&2
  echo "         export TWINE_PASSWORD='pypi-<YOUR_REAL_PYPI_TOKEN>'" >&2
  exit 1
fi

export TWINE_USERNAME="${TWINE_USERNAME:-__token__}"

if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree is not clean. Commit/stash changes before releasing." >&2
  exit 1
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$BRANCH" == "HEAD" ]]; then
  echo "Detached HEAD detected. Checkout a branch before releasing." >&2
  exit 1
fi

VERSION="$(python - <<'PY'
import re
from pathlib import Path
txt = Path('mbctools.py').read_text(encoding='utf-8')
m = re.search(r'^__version__\s*=\s*"([^"]+)"', txt, re.M)
if not m:
    raise SystemExit('Unable to parse __version__ from mbctools.py')
print(m.group(1))
PY
)"
TAG="v${VERSION}"

echo "Preparing release ${VERSION} from branch ${BRANCH}"

python -m pip install --upgrade pip build twine setuptools wheel
rm -rf dist build *.egg-info
python -m build
python -m twine check dist/*

if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "Tag ${TAG} already exists locally."
else
  git tag -a "$TAG" -m "mbctools ${VERSION}"
fi

echo "Pushing git branch and tag before publishing to PyPI..."
git push origin "$BRANCH"
git push origin "$TAG"

echo "Uploading artifacts to PyPI..."
python -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/* --verbose

echo "Release complete: ${VERSION}"
