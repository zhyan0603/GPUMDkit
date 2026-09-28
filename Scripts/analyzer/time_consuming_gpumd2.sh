#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     time_consuming_gpumd.sh
# Category:   Analyzer Scripts
# Purpose:    Monitor GPUMD simulation progress in real time by tracking
#            neighbor.out (or thermo.out), and display speed, total time, and estimated
#            completion time.
# Usage:      ./time_consuming_gpumd.sh
# Output:
#   Real-time table of current frame, speed, total time, time left,
#   and estimated end time
# Author:     Zihan YAN (yanzihan@westlake.edu.cn)
# Last-modified: 2026-09-28
# =============================================================================

# Get the total number of frames from the "run" file
frames=$(awk '$1 == "run" {sum += $2} END {printf "%.0f\n", sum}' run.in)

# Initialize variables
current_frame=0
last_frame=0
speed=0
first_update=true  # Flag to track first update
poll_interval=1       # Seconds between snapshots of the latest progress
averaging_window=30   # Average complete progress intervals over about 30 seconds
warmup_seconds=5      # Avoid estimates based on only a fraction of a second
sample_times=()
sample_frames=()

# Use elapsed time for speed measurement; wall-clock time is only for the ETA.
sample_time() {
    local uptime unused
    if [ -r /proc/uptime ]; then
        read -r uptime unused < /proc/uptime
        printf '%s\n' "$uptime"
    else
        date +%s.%N
    fi
}

# Return success only when there is a new, usable estimate. Keep progress
# endpoints rather than individual log-line arrival times. A long interval
# without output must remain part of the measurement when a batch arrives.
update_speed() {
    local last_index cutoff_time
    last_index=$(( ${#sample_times[@]} - 1 ))
    if [ "$first_update" = true ] || [ "$current_frame" -lt "$last_frame" ] ||
       { [ "$last_index" -ge 0 ] && awk -v now="$current_time" \
            -v previous="${sample_times[$last_index]}" 'BEGIN {exit !(now <= previous)}'; }; then
        sample_times=("$current_time")
        sample_frames=("$current_frame")
        last_frame=$current_frame
        speed=0
        first_update=false
        return 1
    fi
    [ "$current_frame" -gt "$last_frame" ] || return 1

    sample_times+=("$current_time")
    sample_frames+=("$current_frame")
    last_frame=$current_frame
    cutoff_time=$(awk -v now="$current_time" -v window="$averaging_window" \
        'BEGIN {printf "%.6f", now - window}')

    # Retain the endpoint at/before the window boundary. In particular, do
    # not discard the waiting time before a batch that took >30 s to arrive.
    while [ "${#sample_times[@]}" -gt 2 ] &&
          awk -v second="${sample_times[1]}" -v cutoff="$cutoff_time" \
              'BEGIN {exit !(second <= cutoff)}'; do
        sample_times=("${sample_times[@]:1}")
        sample_frames=("${sample_frames[@]:1}")
    done

    speed=$(awk -v now="$current_time" -v start="${sample_times[0]}" \
        -v frame="$current_frame" -v first="${sample_frames[0]}" \
        -v warmup="$warmup_seconds" -v total="$frames" '
        BEGIN {
            elapsed = now - start
            if (elapsed > 0 && (elapsed >= warmup || frame >= total))
                printf "%.6f", (frame - first) / elapsed
        }')
    [ -n "$speed" ]
}

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
    [ "$remaining_frames" -ge 0 ] || remaining_frames=0
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

# Wait for either output file; prefer neighbor.out when both exist.
if [ ! -f "neighbor.out" ] && [ ! -f "thermo.out" ]; then
    echo "Error: neighbor.out does not exist. Waiting for file to appear..."
    until [ -f "neighbor.out" ] || [ -f "thermo.out" ]; do
        sleep 1
    done
fi

if [ -f "neighbor.out" ]; then
    monitor_file="neighbor.out"
else
    monitor_file="thermo.out"
    # Prefer the interval recorded in thermo.out; support older headerless files.
    thermo_interval=$(awk '
        { sub(/\r$/, "") }
        $1 == "#" && $2 == "dump_thermo" && $3 ~ /^[0-9]+$/ && $3 > 0 {
            print $3; exit
        }
    ' thermo.out)
    if [ -z "$thermo_interval" ]; then
        thermo_interval=$(awk '
            { sub(/\r$/, "") }
            $1 == "dump_thermo" && $2 ~ /^[0-9]+$/ && $2 > 0 {
                print $2; exit
            }
        ' run.in)
    fi
    if [ -z "$thermo_interval" ]; then
        echo "Error reading dump_thermo interval from thermo.out or run.in"
        exit 1
    fi
fi

read_current_frame() {
    if [ "$monitor_file" = "neighbor.out" ]; then
        local field latest_frame="" unused
        # Read only a small tail of the file, and use the latest complete
        # progress record. read ignores an unfinished last line during writes.
        while read -r unused unused unused unused field unused; do
            field=${field%$'\r'}
            field=${field%:}
            if [[ "$field" =~ ^[0-9]+$ ]]; then
                latest_frame=$field
            fi
        done < <(tail -n 32 neighbor.out)
        printf '%s\n' "$latest_frame"
    else
        # Count existing data too, so monitoring can start partway through a run.
        awk -v interval="$thermo_interval" '
            { sub(/\r$/, "") }
            NF && $1 !~ /^#/ { rows++ }
            END { printf "%.0f\n", rows * interval }
        ' thermo.out
    fi
}

while true; do
    current_frame=$(read_current_frame)
    current_time=$(sample_time)
    if [[ ! "$current_frame" =~ ^[0-9]+$ ]] || ! update_speed; then
        sleep "$poll_interval"
        continue
    fi

    calculate_times

    # Format speed as a string for centering
    speed_str=$(printf "%.2f" "$speed")

    # Print table row (centered)
    printf "%-15s %-12s %-15s %-15s %-20s\n" \
        "$(center_text "$current_frame" 15)" \
        "$(center_text "$speed_str" 15)" \
        "$(center_text "$total_time" 15)" \
        "$(center_text "$remaining_time_str" 15)" \
        "$(center_text "$end_time" 20)"
    sleep "$poll_interval"
done
