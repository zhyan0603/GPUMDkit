#!/bin/bash

function f302_md_sample_batch_pretreatment(){
    echo " ------------>>"
    echo " Starting MD sample batch pretreatment..."

    # Find all .vasp and .xyz files in the current directory

    num_vasp_files=$(find . -maxdepth 1 -name "*.vasp" | wc -l)
    num_xyz_files=$(find . -maxdepth 1 -name "*.xyz" | wc -l)

	# Check if there are any .vasp files
	if [ $num_vasp_files -gt 0 ]; then
	    # Create the struct directory and move .vasp files into it
	    mkdir -p struct_md
	    rename_seq=1
		for file in *.vasp; do
		    new_vasp_name="POSCAR_${rename_seq}.vasp"
		    new_xyz_name="model_${rename_seq}.xyz"
		    mv "$file" "$new_vasp_name"
		    python ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py ${new_vasp_name} ${new_xyz_name}
		    rename_seq=$((rename_seq + 1))
		done
		mv POSCAR_*.vasp ./struct_md
		mv model_*.xyz ./struct_md
	else
	    # Check if there is exactly one XYZ file
	    if [ $num_xyz_files -eq 1 ]; then
	        echo " No .vasp files found, but found one XYZ file."
	        echo " Converting it to model.xyz using GPUMDkit..."
	        python ${GPUMDkit_path}/Scripts/format_conversion/split_single_xyz.py *.xyz
	        
	        mkdir -p struct_md
	        mv *.xyz ./struct_md
	        num_xyz_files=$(find ./struct_md -maxdepth 1 -name "*.xyz" | wc -l | awk '{print $1-1}')
	        
	        # Perform additional operations if needed after moving .vasp files
	    else
	        echo " No .vasp files found and the XYZ file is not unique."
	        exit 1
	    fi
	fi

    echo " Found $num_xyz_files .xyz files."

    # Ask user for directory name prefix
    echo " >-------------------------------------------------<"
    echo " | This function calls the script in Scripts       |"
    echo " | Script: md_sample_batch_pretreatment.sh         |"
    echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
    echo " >-------------------------------------------------<"

    # Create md directory
    mkdir -p md

    # Create individual directories for each .vasp file and set up the links
    for i in $(seq 1 $num_xyz_files); do
        dir_name="sample_${i}"
        mkdir -p ${dir_name}
        cd ${dir_name}
        ln -s ../struct_md/model_${i}.xyz ./model.xyz
        ln -s ../md/{run.in,nep.txt} ./
        cd ..
    done

    # Create the presub.sh file for VASP self-consistency calculations
    cat > presub.sh <<-EOF
	#!/bin/bash

	# You can cat it to your submit script.

	for dir in sample_*; do
	    cd \$dir
	    echo "Running MD sample in \$dir..."
	    gpumd > log
	    cd ..
	done
	EOF

    # Make presub.sh executable
    chmod +x presub.sh

    echo " >-----------------------------------------------<"
    echo " ATTENTION: Place run.in and nep.txt in 'md' Dir. "
    echo " ATTENTION: Place run.in and nep.txt in 'md' Dir. "
    echo " ATTENTION: Place run.in and nep.txt in 'md' Dir. "
    echo " >-----------------------------------------------<"
}
