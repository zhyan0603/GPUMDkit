#--------------------- clean extra files ----------------------
# This function is used to clean extra files in the current directory
# It will keep some files and delete the rest

function clean_extra_files(){
keep_files=("run.in" "nep.in" "model.xyz" "nep.txt" "train.xyz" "test.xyz")
keep_patterns=("*sub*" "*.sh" "*slurm")
all_files=$(ls)
delete_files=()

for file in $all_files; do
    # Check if the file matches any keep_files or keep_patterns
    keep=false
    for keep_file in "${keep_files[@]}"; do
        if [[ "$file" == "$keep_file" ]]; then
            keep=true
            break
        fi
    done
    for keep_pattern in "${keep_patterns[@]}"; do
        if [[ "$file" == $keep_pattern ]]; then
            keep=true
            break
        fi
    done

    # If the file is not marked to keep, add it to delete_files
    if [ "$keep" == false ]; then
        delete_files+=("$file")
    fi
done

# Display files to delete
if [ ${#delete_files[@]} -eq 0 ]; then
    echo "No files to delete."
    exit 0
else
    echo "The following files will be deleted:"
    echo "---------------------------------------"
    for file in "${delete_files[@]}"; do
        echo -n "$file "
    done
    echo -e "\n---------------------------------------"
fi

# Ask user for confirmation or additional files to keep
echo "Do you want to delete all these files?"
echo "y/yes to delete, n/no to cancel, or input files to keep (separated by spaces):"
read user_input

# Process user input
if [[ "$user_input" == "y" || "$user_input" == "yes" ]]; then
    echo "Deleting the files..."
    for file in "${delete_files[@]}"; do
        rm -f "$file"
    done
    echo "Files deleted."
elif [[ "$user_input" == "n" || "$user_input" == "no" ]]; then
    echo "Operation canceled. No files were deleted."
    exit 0
else
    # Add extra files to keep based on user input
    extra_keep_files=($user_input)
    for extra_file in "${extra_keep_files[@]}"; do
        delete_files=("${delete_files[@]/$extra_file}")
    done

    # Delete remaining files
    if [ ${#delete_files[@]} -eq 0 ]; then
        echo "No files to delete after processing extra keep files."
    else
        echo "Deleting remaining files..."
        for file in "${delete_files[@]}"; do
            if [ -n "$file" ]; then
                rm -f "$file"
            fi
        done
        echo "Files deleted."
    fi
fi
}