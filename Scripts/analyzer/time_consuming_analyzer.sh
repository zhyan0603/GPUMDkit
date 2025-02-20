#!/bin/bash

# Check if the first argument ($1) exists and is not empty
log_file="${1:-log}"  # If $1 is not provided, default to "log"

# Check if the log file exists
if [ ! -f "$log_file" ]; then
    echo "Error: $log_file does not exist. Please provide a valid log file path."
    exit 1
fi

# Get the number of atoms (remove the trailing dot if any)
atom_num=$(grep "Number of atoms" "$log_file" | awk '{print $5}' | sed 's/\.$//')

# Get the step per second from the input file
atom_speed=$(grep "step/second" "$log_file" | tail -1 | awk '{print $6}')

# Calculate the frame speed (frames per second)
frame_speed=$(echo "${atom_speed} ${atom_num}" | awk '{print $1 / $2}')

# Get the total number of frames from the "run" file
frames=$(grep run run.in | awk '{sum += $2} END {print sum}')

# Calculate the estimated time in seconds (total frames / frame speed)
time_in_seconds=$(echo "${frames} ${frame_speed}" | awk '{print $1 / $2}')

# Convert the total time to hours, minutes, and seconds using awk
hours=$(echo "$time_in_seconds" | awk '{print int($1 / 3600)}')
minutes=$(echo "$time_in_seconds" | awk '{print int(($1 % 3600) / 60)}')
seconds=$(echo "$time_in_seconds" | awk '{print int($1 % 60)}')

# Attempt to get the current frame number from the "neighbor.out" file (if it exists)
if [ -f "neighbor.out" ]; then
    current_frame=$(tail -1 neighbor.out | awk '{print $5}' | sed 's/\:$//')

    if [ -z "$current_frame" ]; then
        echo "Error reading current frame from neighbor.out"
        current_frame=0
    fi

    # Calculate the remaining time
    remaining_frames=$(awk "BEGIN {print $frames - $current_frame}")
    remaining_time=$(awk "BEGIN {print $remaining_frames / $frame_speed}")

    # Convert remaining time to hours, minutes, seconds using awk
    remaining_hours=$(echo "$remaining_time" | awk '{print int($1 / 3600)}')
    remaining_minutes=$(echo "$remaining_time" | awk '{print int(($1 % 3600) / 60)}')
    remaining_seconds=$(echo "$remaining_time" | awk '{print int($1 % 60)}')

    # Calculate progress as a percentage using awk
    progress=$((current_frame * 100 / frames))
    bar_length=$((progress / 3))
    progress_bar=$(printf "%-${bar_length}s" "#" | tr ' ' '#') 
    progress_bar="${progress_bar}$(printf "%$((33 - bar_length))s" | tr ' ' '.')"


else
    current_frame=0
    remaining_hours=0
    remaining_minutes=0
    remaining_seconds=0
    progress=0
    progress_bar=""
fi

# Output the results
echo "--------------- Time Consuming Results ---------------"
echo "num of atoms: $atom_num"
echo "atom*step/s : $atom_speed"
echo "timesteps/s : $frame_speed"
echo "total frames: $frames"
echo "total time  : ${hours}h ${minutes}min ${seconds}s"

# If neighbor.out exists, show progress, current frame, and remaining time
if [ -f "neighbor.out" ]; then
    echo "time left   : ${remaining_hours}h ${remaining_minutes}min ${remaining_seconds}s"
    echo -n "Progress Bar: ["
    echo "${progress_bar}] ${progress}%"
fi
