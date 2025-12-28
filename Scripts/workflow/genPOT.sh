#!/bin/bash
#
# This script is part of GPUMDkit.
# Repository: https://github.com/zhyan0603/GPUMDkit
#
# Function to create POTCAR file using vaspkit
create_potcar_vaspkit() {
    echo "Creating POTCAR using vaspkit..."

    # Use vaspkit to generate POTCAR from existing POSCAR
    (echo 103) | vaspkit >> /dev/null 2>&1

    if [ -f "POTCAR" ]; then
        echo "POTCAR successfully created using vaspkit"
    else
        echo "Error: Failed to create POTCAR using vaspkit"
        return 1
    fi
}

# Check if POSCAR exists in current directory
if [ -f "POSCAR" ]; then
    echo "Found POSCAR in current directory, processing..."

    # Check if POTCAR already exists
    if [ -f "POTCAR" ]; then
        echo "POTCAR already exists in current directory. Skipping creation."
    else
        create_potcar_vaspkit
    fi
else
    # Ensure fp directory exists
    mkdir -p input_fp
    # Iterate through all iter* folders
    for folder in iter*; do
        if [ -d "$folder" ] && [ -f "$folder/POSCAR" ]; then
            #echo "Processing $folder/POSCAR..."

            # Extract elements line for filename generation
            elements_line=$(head -n 6 "$folder/POSCAR" | tail -n 1)
            elements=($elements_line)

            # Generate POTCAR filename abbreviation
            potcar_short=""
            for element in "${elements[@]}"; do
                potcar_short+="$element"
            done

            # Generate POTCAR filename
            potcar_name="POTCAR_${potcar_short}"
            potcar_path="input_fp/$potcar_name"

            # Check if POTCAR file already exists
            if [ ! -f "$potcar_path" ]; then
                # echo "Creating $potcar_name"

                # Copy POSCAR to input_fp directory and create POTCAR using vaspkit
                cp "$folder/POSCAR" "input_fp/tmp_POSCAR"
                (cd input_fp && mv tmp_POSCAR POSCAR && create_potcar_vaspkit && mv POTCAR "$potcar_name" && rm -f POSCAR)
            fi

            # Create symbolic link
            ln -sf "../$potcar_path" "$folder/POTCAR"
            #echo "Created symbolic link: $folder/POTCAR -> ../$potcar_path"
        fi
    done
    echo "All POTCAR files have been processed in iter* folders."
fi
