#!/usr/bin/env bash

set -e

bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix="$(cd "${bin_dir}/.." && pwd)"
runtime="${prefix}/share/gpumdkit"

if [ ! -f "${runtime}/gpumdkit.sh" ]; then
    echo " Error: GPUMDkit runtime was not found at ${runtime}."
    exit 1
fi

case "${1:-}" in
    -update|-U)
        echo " GPUMDkit was installed by conda."
        echo " Please update it with: conda update gpumdkit"
        exit 0
        ;;
esac

export GPUMDkit_path="${runtime}"
export PATH="${bin_dir}:${PATH:-}"

exec bash "${runtime}/gpumdkit.sh" "$@"

