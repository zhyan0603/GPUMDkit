#!/usr/bin/env bash

set -euo pipefail

runtime="${PREFIX}/share/gpumdkit"

mkdir -p "${runtime}"
mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/share/bash-completion/completions"

cp -a gpumdkit.sh "${runtime}/"
cp -a src "${runtime}/"
cp -a Scripts "${runtime}/"

mkdir -p "${runtime}/docs"
cp -a docs/updates.info "${runtime}/docs/"

chmod 755 "${runtime}/gpumdkit.sh"
find "${runtime}/src" "${runtime}/Scripts" -type f -name "*.sh" -exec chmod 755 {} +
find "${runtime}/Scripts" -type f -name "*.py" -exec chmod 644 {} +

install -m 755 conda.recipe/gpumdkit-conda-wrapper.sh "${PREFIX}/bin/gpumdkit"
install -m 755 conda.recipe/gpumdkit-conda-wrapper.sh "${PREFIX}/bin/gpumdkit.sh"
install -m 644 Scripts/utils/completion.sh "${PREFIX}/share/bash-completion/completions/gpumdkit"

