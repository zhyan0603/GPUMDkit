#!/bin/bash
#
# This script is part of GPUMDkit.
# Repository: https://github.com/zhyan0603/GPUMDkit
#
# Get the total number of frames from the "run" file
frames=$(grep run run.in | awk '{sum += $2} END {print sum}')

# Initialize variables
current_frame=0
last_frame=0
last_time=$(date +%s.%N)  # Get current time in seconds with nanoseconds
speed=0
first_update=true  # Flag to track first update

# Function to center-align text within a given width
center_text() {
    local text="$1"
    local width="$2"
    local len=${#text}
    if [ $len -ge $width ]; then
        # If text is too long, truncate to fit
        printf "%s" "${text:0:$width}"
        return
    fi
    local padding=$(( (width - len) / 2 ))
    local left_padding=$padding
    local right_padding=$(( width - len - left_padding ))
    printf "%${left_padding}s%s%${right_padding}s" "" "$text" ""
}

# Optimized function to calculate times using awk
calculate_times() {
    remaining_frames=$((frames - current_frame))
    if [ $(awk "BEGIN {print ($speed > 0)}") -eq 1 ]; then
        # Perform all time calculations in a single awk command
        read time_in_seconds hours minutes seconds \
             remaining_time remaining_hours remaining_minutes remaining_seconds \
        < <(awk -v frames="$frames" -v speed="$speed" -v remaining_frames="$remaining_frames" '
            BEGIN {
                # Total time
                time_in_seconds = frames / speed;
                hours = int(time_in_seconds / 3600);
                minutes = int((time_in_seconds % 3600) / 60);
                seconds = int(time_in_seconds % 60);
                
                # Remaining time
                remaining_time = remaining_frames / speed;
                remaining_hours = int(remaining_time / 3600);
                remaining_minutes = int((remaining_time % 3600) / 60);
                remaining_seconds = int(remaining_time % 60);
                
                # Output all values
                print time_in_seconds, hours, minutes, seconds, 
                      remaining_time, remaining_hours, remaining_minutes, remaining_seconds
            }')
        
        total_time="${hours}h ${minutes}m ${seconds}s"
        remaining_time_str="${remaining_hours}h ${remaining_minutes}m ${remaining_seconds}s"

        # Estimated end time
        current_time=$(date +%s)
        end_time_seconds=$(awk -v ct="$current_time" -v rt="$remaining_time" 'BEGIN {print int(ct + rt)}')
        end_time=$(date -d "@$end_time_seconds" "+%Y-%m-%d %H:%M:%S")
    else
        total_time="0h 0m 0s"
        remaining_time_str="0h 0m 0s"
        end_time="N/A"
    fi
}

# Output initial results (printed only once)
echo " ----------------- System Information ----------------"
echo " total frames: $frames"
echo " -----------------------------------------------------"

# Print table header (centered)
printf "%-15s %-12s %-15s %-15s %-20s\n" \
    "$(center_text "Current Frame" 15)" \
    "$(center_text "Speed (steps/s)" 15)" \
    "$(center_text "Total Time" 15)" \
    "$(center_text "Time Left" 15)" \
    "$(center_text "Estimated End" 20)"
printf "%-15s %-12s %-15s %-15s %-20s\n" \
    "$(center_text "-------------" 15)" \
    "$(center_text "-------------" 15)" \
    "$(center_text "-------------" 15)" \
    "$(center_text "-------------" 15)" \
    "$(center_text "-----------------" 20)"

# Check if neighbor.out exists
if [ ! -f "neighbor.out" ]; then
    echo "Error: neighbor.out does not exist. Waiting for file to appear..."
    until [ -f "neighbor.out" ]; do
        sleep 1
    done
fi

# Monitor neighbor.out for updates
# echo "Monitoring neighbor.out for real-time updates..."
tail -f -n 1 neighbor.out | while read -r line; do
    # Extract current frame from the line (remove trailing colon)
    current_frame=$(echo "$line" | awk '{print $5}' | sed 's/\:$//')

    # Validate current_frame
    if [[ ! "$current_frame" =~ ^[0-9]+$ ]]; then
        echo "Error reading current frame from neighbor.out"
        continue
    fi

    # Calculate time difference and speed (assuming each update is 1000 steps)
    current_time=$(date +%s.%N)
    time_diff=$(awk -v ct="$current_time" -v lt="$last_time" 'BEGIN {print ct - lt}')
    if [ $(awk "BEGIN {print ($time_diff > 0)}") -eq 1 ]; then
        steps_diff=$((current_frame - last_frame))
        if [ $steps_diff -gt 0 ]; then
            speed=$(awk -v sd="$steps_diff" -v td="$time_diff" 'BEGIN {print sd / td}')
        fi
    fi

    # Update last frame and time
    last_frame=$current_frame
    last_time=$current_time

    # Calculate times
    calculate_times

    # Skip the first update
    if [ "$first_update" = true ]; then
        first_update=false
        continue
    fi

    # Format speed as a string for centering
    speed_str=$(printf "%.2f" "$speed")

    # Print table row (centered)
    printf "%-15s %-12s %-15s %-15s %-20s\n" \
        "$(center_text "$current_frame" 15)" \
        "$(center_text "$speed_str" 15)" \
        "$(center_text "$total_time" 15)" \
        "$(center_text "$remaining_time_str" 15)" \
        "$(center_text "$end_time" 20)"
done