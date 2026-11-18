#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     md_sample_batch_pretreatment_lmp.sh
# Category:   Workflow Scripts
# Purpose:    Batch pretreatment of structures for LAMMPS MD sampling:
#             rename VASP/XYZ files, generate lammps input and data files,
#             and create struct_md directories.
# Usage:      source md_sample_batch_pretreatment_lmp.sh
# Author:     Zihan YAN (yanzihan@westlake.edu.cn)
# Last-modified: 2026-05-16
# =============================================================================

source "${GPUMDkit_path}/Scripts/workflow/md_sample_batch_pretreatment_gpumd.sh"

function f303_md_sample_batch_pretreatment_lmp(){
    local num_xyz_files i
    echo " ------------>>"
    echo " Starting MD sample batch pretreatment..."
    prepare_md_structures || return 1

    for ((i=1; i<=num_xyz_files; i++)); do
        python "${GPUMDkit_path}/Scripts/format_conversion/exyz2lmp.py" \
            "struct_md/model_${i}.xyz" "struct_md/lammps_${i}.data" || return 1
    done
    echo " $num_xyz_files lammps.data files were generated."

    echo " >-------------------------------------------------<"
    echo " | This function calls the script in Scripts       |"
    echo " | Script: md_sample_batch_pretreatment_lmp.sh     |"
    echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
    echo " >-------------------------------------------------<"

    # Create md directory
    mkdir -p md

    # Create one sample directory per structure and set up the links
    for i in $(seq 1 $num_xyz_files); do
        dir_name="sample_${i}"
        mkdir -p "${dir_name}"
        cd "${dir_name}"
        ln -s ../struct_md/lammps_${i}.data ./lammps.data
        ln -s ../md/{lmprun.in,nep.txt} ./
        cd ..
    done

    # Create the presub.sh file for LAMMPS MD sampling
    cat > presub.sh <<-EOF
	#!/bin/bash

	# You can copy this to your submit script.

	for dir in sample_*; do
	    cd \$dir
	    echo "Running MD sample in \$dir..."
	    mpirun -n X lmp -in lmprun.in > log
	    cd ..
	done
	EOF

    # Make presub.sh executable
    chmod +x presub.sh

    echo " >---------------------------------------------------------<"
    echo " | ATTENTION: Place lmprun.in and nep.txt in 'md' Dir.     |"
    echo " >---------------------------------------------------------<"
}
