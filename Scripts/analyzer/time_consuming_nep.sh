#!/bin/bash

# Get total steps from nep.in (still needed for remaining time calculation)
total_steps=$(grep 'generation' "nep.in" | awk '{print $2}')

# Get the last output step from loss.out (for initializing current step)
current_step=$(tail -1 "loss.out" | awk '{print $1}')

# Initialize variables for time and step tracking
first_time=0
first_step=0
first_read=true

# Print table header with borders
echo "+-----------------+-----------+-----------------+"
echo "|       Step      | Time Diff |    Time Left    |"
echo "+-----------------+-----------+-----------------+"

# Read loss.out file in real-time
tail -1f "loss.out" | while read line; do
  # Extract the step from the current line (assuming step is the first column)
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

  # Prepare centered values
  step_val="$current_step"
  time_val="$time_diff s"
  remain_val="$hours h $minutes m $seconds s"

  # Output data in table format with centered values
  printf "| %-15s | %-9s | %-15s |\n" "$current_step" "$time_diff s" "$hours h $minutes m $seconds s"
done