#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     md_sample_batch_pretreatment_gpumd.sh
# Category:   Workflow Scripts
# Purpose:    Batch pretreatment of structures for GPUMD MD sampling:
#             rename VASP/XYZ files, generate run_*.in and model_*.xyz,
#             and create struct_md directories.
# Usage:      source md_sample_batch_pretreatment_gpumd.sh
# Author:     Zihan YAN (yanzihan@westlake.edu.cn)
# Last-modified: 2026-05-16
# =============================================================================

# Shared by menus 302 and 303. Set num_xyz_files to this run's frame count.
function prepare_md_structures(){
    local -a xyz_files=() vasp_files=() ordered_vasp_files=()
    local file i=1

    for file in ./*.xyz; do
        [ ! -f "$file" ] || xyz_files+=("$file")
    done
    for file in ./*.vasp; do
        [ ! -f "$file" ] || vasp_files+=("$file")
    done
    if [ "${#vasp_files[@]}" -gt 0 ]; then
        if [ "${#xyz_files[@]}" -ne 0 ]; then
            echo " Error: provide one .xyz file or one or more .vasp files."
            return 1
        fi
    elif [ "${#xyz_files[@]}" -ne 1 ]; then
        echo " Error: provide one .xyz file or one or more .vasp files."
        return 1
    fi
    for file in ./sample_[0-9]*; do
        if [[ "$file" =~ ^\./sample_[0-9]+$ ]] && [ -e "$file" ]; then
            echo " Error: $file already exists. Run in a fresh directory."
            return 1
        fi
    done
    if ! mkdir struct_md; then
        echo " Error: cannot create struct_md. Run in a fresh directory."
        return 1
    fi

    if [ "${#vasp_files[@]}" -gt 0 ]; then
        while IFS= read -r file; do
            python "${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py" "$file" "struct_md/model_${i}.xyz" || return 1
            ordered_vasp_files+=("$file")
            i=$((i + 1))
        done < <(printf '%s\n' "${vasp_files[@]}" | sort -V)
        num_xyz_files=$((i - 1))
        for ((i=0; i<num_xyz_files; i++)); do
            mv "${ordered_vasp_files[i]}" "struct_md/POSCAR_$((i + 1)).vasp" || return 1
        done
    else
        file=${xyz_files[0]#./}
        echo " Splitting $file into model_*.xyz files..."
        (cd struct_md && python "${GPUMDkit_path}/Scripts/format_conversion/split_single_xyz.py" "../$file") || return 1
        local -a models=(struct_md/model_*.xyz)
        num_xyz_files=${#models[@]}
        if [ -e "struct_md/$file" ]; then
            echo " Kept the source '$file' in the working directory because its name matches a generated model."
        else
            mv "./$file" struct_md/ || return 1
        fi
    fi
}

function f302_md_sample_batch_pretreatment_gpumd(){
    local num_xyz_files
    echo " ------------>>"
    echo " Starting MD sample batch pretreatment..."
    prepare_md_structures || return 1

    echo " $num_xyz_files model.xyz files were generated."

    echo " >-------------------------------------------------<"
    echo " | This function calls the script in Scripts       |"
    echo " | Script: md_sample_batch_pretreatment_gpumd.sh   |"
    echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
    echo " >-------------------------------------------------<"

    # Create md directory
    mkdir -p md

    # Create one sample directory per structure and set up the links
    for i in $(seq 1 $num_xyz_files); do
        dir_name="sample_${i}"
        mkdir -p "${dir_name}"
        cd "${dir_name}"
        ln -s ../struct_md/model_${i}.xyz ./model.xyz
        ln -s ../md/nep.txt ./
		ln -s ../md/run_${i}.in ./run.in
        cd ..
    done

    # Create the presub.sh file for GPUMD MD sampling
    cat > presub.sh <<-EOF
	#!/bin/bash

	# You can copy this to your submit script.

	for dir in sample_*; do
	    cd \$dir
	    echo "Running MD sample in \$dir..."
	    gpumd > log
	    cd ..
	done
	EOF

    # Make presub.sh executable
    chmod +x presub.sh

    echo " >---------------------------------------------------------<"
    echo " |    ATTENTION: Place run_*.in and nep.txt in 'md' Dir.   |"
    echo " >---------------------------------------------------------<"
    echo " |    You need to provide MD control parameter files in    |"
    echo " |    the format run_*.in (e.g., run_1.in, run_2.in),      |"
    echo " |    each corresponding to a sample (e.g., sample_1,      |"
    echo " |    sample_2) for molecular dynamics simulations.        |"
    echo " >---------------------------------------------------------<"

}
