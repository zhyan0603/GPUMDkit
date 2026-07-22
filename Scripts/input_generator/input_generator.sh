#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     input_generator.sh
# Category:   Input Generator
# Purpose:    Interactively write a commented generic GPUMD run.in template.
# Usage:      gpumdkit.sh -input
# Output:
#   run.in    Commented GPUMD input template for manual editing.
# Author:     Zihan YAN (yanzihan@westlake.edu.cn)
# Last-modified: 2026-07-11
# =============================================================================

function print_usage(){
echo " Usage: gpumdkit.sh -input"
echo ""
echo " This private interactive command writes a generic commented run.in template."
echo " Edit the generated run.in before running GPUMD."
echo ""
}

function read_input_choice(){
    local __target_var="$1"
    local __input
    if ! IFS= read -r __input; then
        echo " Input closed. Exiting."
        return 1
    fi
    printf -v "$__target_var" '%s' "$__input"
}

function write_run_in_header(){
cat > run.in <<'EOF'
# =============================================================================
# GPUMDkit generic run.in template
#
# This file is a template only. Replace every <...> placeholder, then remove
# the leading # from only the commands required by your own protocol.
# Do not run GPUMD with unresolved placeholders or conflicting ensemble lines.
# =============================================================================

# Optional replication. Put it before the potential command when used.
# replicate <n_a> <n_b> <n_c>

# Potential and integration settings.
# potential <potential_filename>
# velocity <initial_temperature_K> seed <seed_number>
# time_step <dt_fs>

EOF
}

function write_minimize_template(){
cat >> run.in <<'EOF'
# =============================================================================
# Optional minimization stage
# =============================================================================
# minimize fire <force_tolerance_eV_per_Angstrom> <max_steps>

EOF
}

function write_ensemble_template(){
cat >> run.in <<'EOF'
# =============================================================================
# MD production stage: uncomment exactly one ensemble for this run block
# =============================================================================
# ensemble nve
# ensemble nvt_nhc <T_initial_K> <T_final_K> <T_coupling_steps>
# ensemble nvt_ber <T_initial_K> <T_final_K> <T_coupling_steps>
# ensemble npt_ber <T_initial_K> <T_final_K> <T_coupling_steps> <P_hydro_GPa> <C_hydro_GPa> <P_coupling_steps>
# ensemble npt_scr <T_initial_K> <T_final_K> <T_coupling_steps> <P_hydro_GPa> <C_hydro_GPa> <P_coupling_steps>

EOF
}

function write_output_template(){
cat >> run.in <<'EOF'
# Standard outputs for the selected production stage.
# dump_thermo <interval_steps>
# dump_exyz <interval_steps> <has_velocity_0_or_1> <has_force_0_or_1>
# dump_restart <interval_steps>

EOF
}

function write_compute_template(){
cat >> run.in <<'EOF'
# Optional analysis commands for the selected production stage.
# compute_msd <sample_interval_steps> <correlation_samples> all_groups <grouping_method>
# compute_rdf <cutoff_Angstrom> <num_bins> <sample_interval_steps>
# compute_adf <sample_interval_steps> <num_bins> <rc_min_Angstrom> <rc_max_Angstrom>

# Run length for this stage.
# run <number_of_steps>

# Copy the MD production block below when a second stage is required. Each
# ensemble and compute command is run-scoped, so write it again in that block.
EOF
}

function f_input_gpumd_run_in(){
echo " Generating the generic GPUMD run.in template."
write_run_in_header
write_minimize_template
write_ensemble_template
write_output_template
write_compute_template
echo " Generated template: run.in"
echo " Edit run.in and uncomment only the commands required by your protocol."
}

function input_generator_menu(){
echo " >-------------------------------------------------<"
echo " | GPUMDkit input generator (private preview)     |"
echo " >-------------------------------------------------<"
echo " | 1) Generic GPUMD run.in template                |"
echo " | 0) Exit                                          |"
echo " >-------------------------------------------------<"
echo " ------------>>"
read_input_choice input_choice || return 1

case "$input_choice" in
    "1") f_input_gpumd_run_in ;;
    "0") echo " Operation canceled." ;;
    *) echo " Error: invalid input generator option."; return 1 ;;
esac
}

case "$1" in
    -h|--help) print_usage ;;
    "") input_generator_menu ;;
    *)
        print_usage
        exit 1 ;;
esac
