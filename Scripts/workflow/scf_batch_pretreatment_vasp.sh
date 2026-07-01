#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     scf_batch_pretreatment_vasp.sh
# Category:   Workflow Scripts
# Purpose:    Batch pretreatment of structures for VASP SCF calculations:
#             rename VASP/XYZ files, generate POSCAR files, create
#             struct_fp directories for INCAR/POTCAR/KPOINTS.
# Usage:      source scf_batch_pretreatment_vasp.sh
# Author:     Zihan YAN (yanzihan@westlake.edu.cn)
# Last-modified: 2026-05-16
# =============================================================================

function vasp_scf_batch_pretreatment(){
    echo " ------------>>"
    echo " Starting SCF batch pretreatment..."

    # Find all .vasp and .xyz files in the current directory

    num_vasp_files=$(find . -maxdepth 1 -name "*.vasp" | wc -l)
    num_xyz_files=$(find . -maxdepth 1 -name "*.xyz" | wc -l)

	# Check if there are any .vasp files
	if [ $num_vasp_files -gt 0 ]; then
	    if [ $num_xyz_files -gt 0 ]; then
	        echo " Notice: Found both .vasp and .xyz files in the current directory."
	        echo " This workflow prioritizes .vasp files and will ignore .xyz files."
	    fi
	    # Create the struct directory and move .vasp files into it
	    mkdir -p struct_fp
	    rename_seq=1
		for file in $(ls -v *.vasp); do
		    new_name="POSCAR_${rename_seq}.vasp"
		    mv "$file" ./struct_fp/"$new_name"
		    rename_seq=$((rename_seq + 1))
		done
        num_vasp_files=$(find ./struct_fp -maxdepth 1 -name "*.vasp" | wc -l)
	else
	    # Check available XYZ files
	    if [ $num_xyz_files -ge 1 ]; then
	        if [ $num_xyz_files -eq 1 ]; then
	            xyz_file=$(find . -maxdepth 1 -name "*.xyz" | head -n 1)
	            echo " No .vasp files found, but found one XYZ file: ${xyz_file#./}"
	        else
	            echo " No .vasp files found, but found multiple XYZ files:"
	            select xyz_file in *.xyz; do
	                if [ -n "$xyz_file" ]; then
	                    break
	                fi
	                echo " Invalid selection. Please choose one XYZ file."
	            done
	            [ -n "$xyz_file" ] || { echo " Input closed. Exiting."; return 1; }
	        fi
	        echo " Converting ${xyz_file#./} to POSCAR files using GPUMDkit..."
	        python ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py "$xyz_file"
	        
	        mkdir -p struct_fp
	        mv *.vasp ./struct_fp
	        num_vasp_files=$(find ./struct_fp -maxdepth 1 -name "*.vasp" | wc -l)
	        
	        # Perform additional operations if needed after moving .vasp files
	    else
	        echo " No .vasp files or .xyz files found."
	        exit 1
	    fi
	fi

    echo " Found $num_vasp_files .vasp files."

    # Ask user for directory name prefix
    echo " >-------------------------------------------------<"
    echo " | This function calls the script in Scripts       |"
    echo " | Script: scf_batch_pretreatment_vasp.sh          |"
    echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
    echo " >-------------------------------------------------<"
    echo " We recommend using the prefix to locate the structure."
    echo " The folder name will be added to the second line of XYZ."
    echo " config_type=<prefix>_<ID>"
    echo " ------------>>"
    echo " Please enter the prefix of directory (e.g. FAPBI3_iter01)"
    read_menu_choice prefix || return 1

    # Create fp directory
    mkdir -p fp

    # Create individual directories for each .vasp file and set up the links
    for i in $(seq 1 $num_vasp_files); do
        dir_name="${prefix}_${i}"
        mkdir -p ${dir_name}
        cd ${dir_name}
        ln -s ../struct_fp/POSCAR_${i}.vasp ./POSCAR
        ln -s ../fp/{POTCAR,KPOINTS,INCAR} ./
        cd ..
    done

    # Create the presub.sh file for VASP self-consistency calculations
    cat > presub.sh <<-EOF
	#!/bin/bash

	# You can cat it to your submit script.

	for dir in ${prefix}_*; do
	    cd \$dir
	    echo "Running VASP SCF in \$dir..."
	    mpirun -n X vasp_std > log
	    cd ..
	done
	EOF

    # Make presub.sh executable
    chmod +x presub.sh

    echo " >---------------------------------------------------------<"
    echo " | ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir. |"
    echo " | ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir. |"
    echo " | ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir. |"
    echo " >---------------------------------------------------------<"
}
