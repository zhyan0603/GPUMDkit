#!/usr/bin/env bash

set -euo pipefail

test -f "${PREFIX}/share/gpumdkit/gpumdkit.sh"
test -d "${PREFIX}/share/gpumdkit/src"
test -d "${PREFIX}/share/gpumdkit/Scripts"
test -f "${PREFIX}/share/gpumdkit/Scripts/workflow/cp2k_template.inp"
test -f "${PREFIX}/share/gpumdkit/skills/README.md"
for skill in gpumdkit-skill gpumdkit-skill-zh; do
    test -f "${PREFIX}/share/gpumdkit/skills/${skill}/SKILL.md"
    test -f "${PREFIX}/share/gpumdkit/skills/${skill}/references/overview.md"
done

env -u GPUMDkit_path gpumdkit -skill | grep -Fq "${PREFIX}/share/gpumdkit/skills/gpumdkit-skill"

env -u GPUMDkit_path gpumdkit -h >/tmp/gpumdkit-help.txt
env -u GPUMDkit_path gpumdkit.sh -h >/tmp/gpumdkit-help-sh.txt
env -u GPUMDkit_path gpumdkit -doctor >/tmp/gpumdkit-doctor.txt
env -u GPUMDkit_path gpumdkit -update | grep -q "conda update gpumdkit"
