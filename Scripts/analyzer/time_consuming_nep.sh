#!/bin/bash
#
# This script is part of GPUMDkit.
# Repository: https://github.com/zhyan0603/GPUMDkit
#
# Description:
#     Analyze time consumption for NEP training
#
# Usage:
#     bash time_consuming_nep.sh [arguments]
#
# Author: Zihan YAN
# Contact: yanzihan@westlake.edu.cn
# Last Modified: 2025-12-28
#

if [ -f "nep.in" ]; then
    total_steps=$(grep 'generation' "nep.in" | awk '{print $2}')
elif [ -f "gnep.in" ]; then
    total_steps=$(grep 'epoch' "gnep.in" | awk '{print $2}')
else
    echo "Error: Neither nep.in nor gnep.in found"
    exit 1
fi

# Get the last output step from loss.out
current_step=$(tail -1 "loss.out" | awk '{print $1}')

# Initialize variables for time and step tracking
first_time=0
first_step=0
first_read=true

# Print table header with borders
echo "+-----------------+-----------+-----------------+---------------------+"
echo "|       Step      | Time Diff |    Time Left    |    Finish Time      |"
echo "+-----------------+-----------+-----------------+---------------------+"

# Read loss.out file in real-time
tail -1f "loss.out" | while read line; do
  # Extract the step from the current line
  step=$(echo "$line" | awk '{print $1}')

  # Get current timestamp
  timestamp=$(date +%s)

  # If it's the first read, record the time and step
  if $first_read; then
    first_time=$timestamp
    first_step=$step
    first_read=false
    continue
  fi

  # Calculate time difference and step difference
  time_diff=$((timestamp - first_time))
  step_diff=$((step - first_step))

  # If step difference is zero, skip this iteration
  if [ "$step_diff" -eq 0 ]; then
    continue
  fi

  # Update current step and time
  current_step=$step
  first_time=$timestamp
  first_step=$step

  # Calculate remaining steps and time
  remaining_steps=$((total_steps - current_step))
  remaining_time=$(((remaining_steps * time_diff) / step_diff))

  # Convert remaining time to hours, minutes, seconds
  hours=$((remaining_time / 3600))
  minutes=$(((remaining_time % 3600) / 60))
  seconds=$((remaining_time % 60))

  # Calculate finish time (current time + remaining time)
  finish_timestamp=$((timestamp + remaining_time))
  finish_time=$(date -d "@$finish_timestamp" "+%Y-%m-%d %H:%M:%S")

  # Output data in table format
  printf "| %-15s | %-9s | %-15s | %-19s |\n" \
    "$current_step" "$time_diff s" "$hours h $minutes m $seconds s" "$finish_time"
done