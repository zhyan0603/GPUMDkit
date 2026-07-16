#!/usr/bin/env bash

set -euo pipefail

test -f "${PREFIX}/share/gpumdkit/gpumdkit.sh"
test -d "${PREFIX}/share/gpumdkit/src"
test -d "${PREFIX}/share/gpumdkit/Scripts"
test -f "${PREFIX}/share/gpumdkit/Scripts/workflow/cp2k_template.inp"

env -u GPUMDkit_path gpumdkit -h >/tmp/gpumdkit-help.txt
env -u GPUMDkit_path gpumdkit.sh -h >/tmp/gpumdkit-help-sh.txt
env -u GPUMDkit_path gpumdkit -doctor >/tmp/gpumdkit-doctor.txt
env -u GPUMDkit_path gpumdkit -update | grep -q "conda update gpumdkit"

