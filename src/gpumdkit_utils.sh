#!/bin/bash

# GPUMDkit Common Utilities
# This file contains common utility functions to reduce code duplication

# Print a formatted header for script execution info
print_script_header() {
    local title="$1"
    local developer="$2"
    local script_path="$3"
    
    echo " >-------------------------------------------------<"
    echo " | Calling the script in Scripts/${title:-general}            |"
    if [ -n "$developer" ]; then
        printf " | Developer: %-46s|\n" "$developer"
    fi
    if [ -n "$script_path" ]; then
        printf " | Script: %-48s|\n" "${script_path##*/}"
    fi
    echo " >-------------------------------------------------<"
}

# Print a formatted footer with code path
print_script_footer() {
    local code_path="$1"
    
    echo " ---------------------------------------------------"
    if [ -n "$code_path" ]; then
        echo " Code path: ${code_path}"
    fi
    echo " ---------------------------------------------------"
}

# Validate that a variable is not empty
validate_not_empty() {
    local value="$1"
    local var_name="$2"
    
    if [ -z "$value" ]; then
        echo "Error: ${var_name} is required but was not provided."
        return 1
    fi
    return 0
}

# Validate that a file exists
validate_file_exists() {
    local filepath="$1"
    
    if [ ! -f "$filepath" ]; then
        echo "Error: File '${filepath}' does not exist."
        return 1
    fi
    return 0
}

# Validate that a directory exists
validate_dir_exists() {
    local dirpath="$1"
    
    if [ ! -d "$dirpath" ]; then
        echo "Error: Directory '${dirpath}' does not exist."
        return 1
    fi
    return 0
}

# Print usage error message with example
print_usage_error() {
    local usage="$1"
    local example="$2"
    local code_path="${3:-}"
    
    echo "Usage: ${usage}"
    if [ -n "$example" ]; then
        echo "Example: ${example}"
    fi
    if [ -n "$code_path" ]; then
        echo "Code path: ${code_path}"
    fi
}

# Read input with validation against a list of allowed values
read_validated_choice() {
    local prompt="$1"
    local allowed_values="$2"  # Space-separated list
    local result_var="$3"      # Variable name to store result
    
    local choice
    read -p "$prompt" choice
    
    while ! echo "$allowed_values" | grep -wq "$choice"; do
        echo "Invalid choice. Please select from: ${allowed_values}"
        read -p "$prompt" choice
    done
    
    eval "$result_var='$choice'"
}

# Print a simple info message box
print_info_box() {
    local message="$1"
    local width=50
    
    local border=$(printf '%*s' "$width" | tr ' ' '-')
    echo " $border"
    printf " | %-${width}s|\n" "$message"
    echo " $border"
}

# Print a warning message
print_warning() {
    local message="$1"
    echo "Warning: ${message}"
}

# Print an error message and optionally exit
print_error() {
    local message="$1"
    local exit_code="${2:-}"
    
    echo "Error: ${message}" >&2
    
    if [ -n "$exit_code" ]; then
        exit "$exit_code"
    fi
}

# Check if running in interactive mode
is_interactive() {
    [ -t 0 ]
}

# Get current timestamp
get_timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}

# Simple logging function (file-based, minimal)
log_message() {
    local level="$1"
    local message="$2"
    local log_file="${3:-gpumdkit.log}"
    
    echo "[$(get_timestamp)] [${level}] ${message}" >> "$log_file"
}
