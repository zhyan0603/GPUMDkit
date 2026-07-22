#--------------------- clean extra files ----------------------
# This function is used to clean extra files in the current directory
# It will keep some files and delete the rest

function clean_extra_files(){
keep_files=("run.in" "nep.in" "model.xyz" "nep.txt" "train.xyz" "test.xyz")
keep_patterns=("*sub*" "*.sh" "*slurm")
delete_files=()

# Build the deletion list using a glob (avoids parsing `ls` output),
# considering only regular files (matches the net effect of plain `ls`,
# where subdirectory names would have failed `rm -f` anyway).
for file in *; do
    [ -e "$file" ] || continue
    [ -f "$file" ] || continue

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
    echo " No files to delete."
    return 0
else
    echo " The following files will be deleted:"
    echo " ----------------------------------------------------"
    for file in "${delete_files[@]}"; do
        echo "   ${file}"
    done
    echo " ----------------------------------------------------"
fi

# Ask user for confirmation or additional files to keep
echo " +------------------------------------------------------+"
echo " | Do you want to delete all these files?              |"
echo " |   y / yes : delete all                              |"
echo " |   n / no  : cancel                                   |"
echo " |   or input filenames to keep (space-separated)       |"
echo " +------------------------------------------------------+"
if ! IFS= read -r -p " " user_input; then
    echo " Input closed. No files were deleted."
    return 1
fi

# Process user input
if [[ "$user_input" == "y" || "$user_input" == "yes" ]]; then
    echo " Deleting the files..."
    for file in "${delete_files[@]}"; do
        rm -f "$file"
    done
    echo " Files deleted."
elif [[ "$user_input" == "n" || "$user_input" == "no" ]]; then
    echo " Operation canceled. No files were deleted."
    return 0
else
    # Add extra files to keep based on user input.
    # Tokenize and exclude by whole-word match (NOT substring removal,
    # which would corrupt names sharing a substring with the keep input).
    read -r -a extra_keep_files <<< "$user_input"
    new_delete_files=()
    for df in "${delete_files[@]}"; do
        skip=false
        for ek in "${extra_keep_files[@]}"; do
            if [ "$df" = "$ek" ]; then
                skip=true
                break
            fi
        done
        if [ "$skip" = false ]; then
            new_delete_files+=("$df")
        fi
    done
    delete_files=("${new_delete_files[@]}")

    # Delete remaining files
    if [ ${#delete_files[@]} -eq 0 ]; then
        echo " No files to delete after processing extra keep files."
    else
        echo " Deleting remaining files..."
        for file in "${delete_files[@]}"; do
            rm -f "$file"
        done
        echo " Files deleted."
    fi
fi
}