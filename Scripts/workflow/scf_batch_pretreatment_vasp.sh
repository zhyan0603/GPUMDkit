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
# Last-modified: 2026-09-26
# =============================================================================

function vasp_scf_batch_pretreatment(){
    local xyz_input=0
    local group_by_element_set=0
    local group_directory group_key matched_group_count
    local -a frame_element_groups
    local -a group_directories
    frame_element_groups=()
    group_directories=()

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
	        xyz_input=1
	        echo " Converting ${xyz_file#./} to POSCAR files using GPUMDkit..."
	        python "${GPUMDkit_path}/Scripts/workflow/exyz2pos_scf.py" "$xyz_file" ./struct_fp || {
            echo " Error: failed to convert ${xyz_file#./} to POSCAR files."
            return 1
        }
	        num_vasp_files=$(find ./struct_fp -type f -name "POSCAR_*.vasp" | wc -l)
	        if [ "$num_vasp_files" -eq 0 ]; then
	            echo " Error: no POSCAR files were created."
	            return 1
	        fi

	        # Perform additional operations if needed after moving .vasp files
	    else
	        echo " No .vasp files or .xyz files found."
	        return 1
	    fi
	fi

    # The Python converter owns composition detection and POSCAR species ordering.
    # Here, use its output directories only to map each frame to its POTCAR key.
    if [ "$xyz_input" -eq 1 ]; then
        for group_directory in ./struct_fp/*/; do
            if [ -d "$group_directory" ]; then
                group_directories+=("${group_directory%/}")
            fi
        done

        if [ "${#group_directories[@]}" -gt 0 ]; then
            group_by_element_set=1
            for i in $(seq 1 "$num_vasp_files"); do
                matched_group_count=0
                for group_directory in "${group_directories[@]}"; do
                    if [ -f "${group_directory}/POSCAR_${i}.vasp" ]; then
                        group_key="${group_directory##*/}"
                        frame_element_groups[$i]="$group_key"
                        matched_group_count=$((matched_group_count + 1))
                    fi
                done
                if [ "$matched_group_count" -ne 1 ]; then
                    echo " Error: could not map POSCAR_${i}.vasp to exactly one element group."
                    return 1
                fi
            done
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
        mkdir -p "${dir_name}"
        cd "${dir_name}"
        if [ "$group_by_element_set" -eq 1 ]; then
            species_key="${frame_element_groups[$i]}"
            ln -s "../struct_fp/${species_key}/POSCAR_${i}.vasp" ./POSCAR
            ln -s "../fp/POTCAR_${species_key}" ./POTCAR
        else
            ln -s "../struct_fp/POSCAR_${i}.vasp" ./POSCAR
            ln -s ../fp/POTCAR ./POTCAR
        fi
        ln -s ../fp/INCAR ./INCAR
        ln -s ../fp/KPOINTS ./KPOINTS
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

    if [ "$group_by_element_set" -eq 1 ]; then
        echo " Prepare the listed POTCAR_<elements> files in fp/."
        echo " All groups share fp/INCAR; match each POTCAR order to POSCAR."
    else
        echo " Prepare the shared fp/POTCAR and fp/INCAR."
    fi
    if [ -f fp/KPOINTS ]; then
        echo " KPOINTS links point to shared fp/KPOINTS."
    else
        echo " KPOINTS links point to fp/KPOINTS, which is currently absent."
        echo " Without that file, VASP uses KSPACING from INCAR."
    fi
}
