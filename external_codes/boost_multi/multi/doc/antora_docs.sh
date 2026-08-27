#!/bin/bash

set -ex

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_ROOT="$SCRIPT_DIR/.."

cd "$REPO_ROOT"

# Step 1 (optional): Run MrDocs to generate AsciiDoc API reference.
# Requires mrdocs and clang++ on PATH. Skip gracefully if not available.
if command -v mrdocs &> /dev/null && command -v clang++ &> /dev/null; then
  cmake -B build-mrdocs \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DMULTI_BUILD_TESTS=OFF \
    -DMULTI_BUILD_PACKAGE=OFF
  mrdocs mrdocs.yml
else
  echo "mrdocs or clang++ not found — skipping API reference generation."
  echo "Install MrDocs from https://github.com/cppalliance/mrdocs/releases"
fi

# Step 3: Build the Antora site
cd "$SCRIPT_DIR"
npm ci --ignore-scripts
npx antora multi-playbook.yml
