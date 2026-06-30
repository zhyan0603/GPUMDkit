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

function f303_md_sample_batch_pretreatment_lmp(){
    echo " ------------>>"
    echo " Starting MD sample batch pretreatment..."

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
	    mkdir -p struct_md
	    rename_seq=1
	    total_vasp_num=$(ls -v *.vasp| wc -l)
		for file in $(ls -v *.vasp); do
		    new_vasp_name="POSCAR_${rename_seq}.vasp"
		    new_xyz_name="model_${rename_seq}.xyz"
			new_lmp_name="lammps_${rename_seq}.data"
		    mv ${file} ./struct_md/${new_vasp_name}
		    python ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py ./struct_md/${new_vasp_name} ./struct_md/${new_xyz_name}
			python ${GPUMDkit_path}/Scripts/format_conversion/exyz2lmp.py ./struct_md/${new_xyz_name} ./struct_md/${new_lmp_name}
		    progress=$((rename_seq * 100 / total_vasp_num))
		    echo -ne " Progress: ["
		    for ((p=0; p<progress/2; p++)); do echo -ne "#"; done
		    for ((p=progress/2; p<50; p++)); do echo -ne "."; done
		    echo -ne "] $progress% ($rename_seq/$total_vasp_num)\r"
		    rename_seq=$((rename_seq + 1))
		done
		num_lmp_files=$(find ./struct_md -maxdepth 1 -name "*.data" | wc -l)
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
	        fi
	        echo " Splitting ${xyz_file#./} to model_*.xyz using GPUMDkit..."
	        python ${GPUMDkit_path}/Scripts/format_conversion/split_single_xyz.py "$xyz_file"
	        
	        mkdir -p struct_md
	        mv *.xyz ./struct_md
			total_xyz_num=$(ls -v ./struct_md/model_*.xyz| wc -l)
			for file in $(ls -v ./struct_md/model_*.xyz); do
				seq=$(basename "$file" .xyz | sed 's/model_//')
    			new_lmp_name="lammps_${seq}.data"
				python ${GPUMDkit_path}/Scripts/format_conversion/exyz2lmp.py ./${file} ./struct_md/${new_lmp_name}
				progress=$((seq * 100 / total_xyz_num))
				echo -ne " Progress: ["
				for ((p=0; p<progress/2; p++)); do echo -ne "#"; done
				for ((p=progress/2; p<50; p++)); do echo -ne "."; done
				echo -ne "] $progress% ($seq/$total_xyz_num)\r"
			done
	        num_lmp_files=$(find ./struct_md -maxdepth 1 -name "*.data" | wc -l)
	        
	        # Perform additional operations if needed after moving .vasp files
	    else
	        echo " No .vasp files or .xyz files found."
	        exit 1
	    fi
	fi

    echo " $num_lmp_files lammps.data files were generated."

    # Ask user for directory name prefix
    echo " >-------------------------------------------------<"
    echo " | This function calls the script in Scripts       |"
    echo " | Script: md_sample_batch_pretreatment_lmp.sh     |"
    echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
    echo " >-------------------------------------------------<"

    # Create md directory
    mkdir -p md

    # Create individual directories for each .vasp file and set up the links
    for i in $(seq 1 $num_lmp_files); do
        dir_name="sample_${i}"
        mkdir -p ${dir_name}
        cd ${dir_name}
        ln -s ../struct_md/lammps_${i}.data ./lammps.data
        ln -s ../md/{lmprun.in,nep.txt} ./
        cd ..
    done

    # Create the presub.sh file for VASP self-consistency calculations
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

    echo " >--------------------------------------------------<"
    echo " ATTENTION: Place lmprun.in and nep.txt in 'md' Dir. "
    echo " ATTENTION: Place lmprun.in and nep.txt in 'md' Dir. "
    echo " ATTENTION: Place lmprun.in and nep.txt in 'md' Dir. "
    echo " >--------------------------------------------------<"
}
