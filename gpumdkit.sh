#!/bin/bash

# GPUMDkit - A User-Friendly Toolkit for GPUMD and NEP
# Project URL: https://github.com/zhyan0603/GPUMDkit

# Description: Main entry point for GPUMDkit interactive and command-line interface

# Copyright (c) 2024-2026 Zihan YAN and GPUMDkit contributors
# License: GPL-3.0 License
# Contact Zihan YAN (yanzihan@westlake.edu.cn) if you have any questions or suggestions!

# You need to set the path of GPUMD and GPUMDkit in your ~/.bashrc, for example
# export GPUMDkit_path=/home/yanzihan/software/GPUMDkit

if [ -z "$GPUMDkit_path" ]; then
    echo " Error: GPUMDkit_path is not set."
    echo " Please set it in your ~/.bashrc, e.g.:"
    echo "   export GPUMDkit_path=/home/yanzihan/software/GPUMDkit"
    exit 1
fi

VERSION="1.5.6 (dev) (2026-07-10)"

plt_path="${GPUMDkit_path}/Scripts/plt_scripts"
analyzer_path="${GPUMDkit_path}/Scripts/analyzer"
calc_path="${GPUMDkit_path}/Scripts/calculators"
workflow_path="${GPUMDkit_path}/Scripts/workflow"
format_conv_path="${GPUMDkit_path}/Scripts/format_conversion"
utils_path="${GPUMDkit_path}/Scripts/utils"
sample_path="${GPUMDkit_path}/Scripts/sample_structures"

#--------------------- main script ----------------------
# Show the menu
function menu(){
echo " ---------------------- GPUMD ------------------------"
echo " 1) Format Conversion          2) Sample Structures   "
echo " 3) Workflow                   4) Calculators         "
echo " 5) Analyzer                   6) Visualization       "
echo " 7) Utilities                  8) Help                "
echo " 0) Exit                                              "
}

function read_menu_choice(){
    local __target_var="$1"
    local __input
    if ! IFS= read -r -p " " __input; then
        echo " Input closed. Exiting."
        return 1
    fi
    printf -v "$__target_var" '%s' "$__input"
}

function read_menu_array(){
    if ! read -r -a "$1"; then
        echo " Input closed. Exiting."
        return 1
    fi
}

# Function main
function main(){
    echo " ------------>>"
    echo ' Input the function number:'
    valid_menu_choices=(
        "0" "1" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110"
        "2" "201" "202" "203" "204" "205" "206"
        "3" "301" "302" "303" 
        "4" "401" "402" "403" "404" "405" "406" "407" "408" "409" "410" "411" "412"
        "5" "501" "502" "503" "504" "505" "506" "507" "508"
        "6"
        "7" "701"
        "8"
    ) 
    read_menu_choice choice || return 1
    while ! echo "${valid_menu_choices[@]}" | grep -wq "$choice" 
    do
      echo " ------------>>"
      echo " Please reinput function number:"
      read_menu_choice choice || return 1
    done

    case "${choice:0:1}" in
        "0")
            citation
            exit 0
            ;;
        "1")
            source ${GPUMDkit_path}/src/f1_format_conversions.sh
            case $choice in
                "1") f1_format_conversion ;;
                "101") f101_out2xyz ;;
                "102") f102_mtp2xyz ;;
                "103") f103_cp2k2xyz ;;
                "104") f104_abacus2xyz ;;
                "105") f105_extxyz2poscar ;;
                "106") f106_add_group_labels ;;
                "107") f107_add_weight ;;
                "108") f108_get_frame ;;
                "109") f109_clean_xyz ;;
                "110") f110_replicate_structure ;;
            esac ;;
        "2")
            source ${GPUMDkit_path}/src/f2_sample_structures.sh
            case $choice in
                "2") f2_sample_structures ;;
                "201") f201_sample_structures ;; 
                "202") f202_pynep_sample_structures ;;
                "203") f203_neptrain_sample_structures ;;
                "204") f204_perturb_structure ;;
                "205") f205_select_max_force_deviation_structs ;;
                "206") f206_split_train_test ;;
            esac ;;
        "3")
            source ${GPUMDkit_path}/src/f3_workflows.sh
            case $choice in
                "3") f3_workflow_dev ;;
                "301") 
                    f301_scf_batch_pretreatment ;;
                "302") 
                    source ${workflow_path}/md_sample_batch_pretreatment_gpumd.sh
                    f302_md_sample_batch_pretreatment_gpumd ;;
                "303") 
                    source ${workflow_path}/md_sample_batch_pretreatment_lmp.sh
                    f303_md_sample_batch_pretreatment_lmp ;;
            esac ;;
        "4")
            source ${GPUMDkit_path}/src/f4_calculators.sh
            case $choice in
                "4") f4_calculators ;;
                "401") f401_calc_ionic_conductivity ;;
                "402") f402_calc_properties_with_nep ;;
                "403") f403_calc_descriptors ;;
                "404") f404_calc_doas ;;
                "405") f405_calc_neb ;;
                "406") f406_calc_neighbor_list ;;
                "407") f407_calc_displacement ;;
                "408") f408_calc_averaged_structure ;;
                "409") f409_calc_oct_tilt ;;
                "410") f410_calc_polarization_abo3 ;;
                "411") f411_minimize_structure_by_nep ;;
                "412") f412_calc_msd_from_trajectory ;;
            esac ;;           
        "5")
            source ${GPUMDkit_path}/src/f5_analyzers.sh
            case $choice in
                "5") f5_analyzers ;;
                "501") f501_analyze_composition ;;
                "502") f502_find_outliers ;;
                "503") f503_analyze_chem_species ;;
                "504") f504_charge_balance_check ;;
                "505") f505_energy_force_virial_analyzer ;;
                "506") f506_filter_structures_by_distance ;;
                "507") f507_get_min_dist ;;
                "508") f508_probability_density_analysis ;;
            esac ;;  
        "6")
            source ${GPUMDkit_path}/src/f6_plots.sh
            terminal_width=$(tput cols 2>/dev/null || echo 80)
            if [ "$terminal_width" -ge 98 ]; then
                f6_plots_two_column
            else
                f6_plots_one_column
            fi ;;
        "7")
            source ${GPUMDkit_path}/src/f7_utilities.sh
            case $choice in
                "7") f7_utilities ;;
                "701") f701_time_consuming_analyzer ;;
            esac ;;
        "8")
            echo " ------------>>"
            echo " GPUMDkit Version ${VERSION}"
            echo " ------------------------------------------------------------"
            help_info_table
            ;;
        *)
            echo " Incorrect Options"
            ;;

    esac
}

#--------------------- help info ----------------------
# This function is used to show the help information
# It will show the usage of each function

function help_info_table(){
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                          GPUMDkit ${VERSION} Command Help                               |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                                          MAIN FUNCTIONS                                               |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | -h            Show this help table            | -plt <type>        Plot and visualization tools       |"
    echo " | -calc <type>  Calculator tools                | -time <gpumd|nep>  Time-consuming analyzer            |"
    echo " | -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |"
    echo " | -skill        Show GPUMDkit agent skill info  | -doctor            Check Python environment           |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                                         FORMAT CONVERSION                                             |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | -out2xyz      OUTCAR -> extxyz (shell)        | -out2exyz          OUTCAR -> extxyz (python)          |"
    echo " | -cp2k2xyz     CP2K log -> xyz                 | -xdat2exyz         XDATCAR -> extxyz                  |"
    echo " | -cif2pos      cif -> POSCAR                   | -cif2exyz          cif -> extxyz                      |"
    echo " | -pos2exyz     POSCAR -> extxyz                | -exyz2pos          extxyz -> POSCAR                   |"
    echo " | -pos2lmp      POSCAR -> LAMMPS data           | -lmp2exyz          LAMMPS dump -> extxyz              |"
    echo " | -traj2exyz    ASE traj -> extxyz              | -replicate         Replicate structure                |"
    echo " | -addgroup     Add group labels                | -addweight         Add structure weight in extxyz     |"
    echo " | -clean_xyz    Clean extra info in extxyz      | -get_frame         Extract specific frame             |"
    echo " | -frame_range  Extract frames by range         | -dp2xyz            DeepMD npy -> extxyz               |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                                            ANALYSIS                                                   |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | -range        Energy/force/virial statistics  | -analyze_comp      Analyze composition                |"
    echo " | -chem_species Analyze chemical species        | -cbc               Charge balance check               |"
    echo " | -min_dist     Min distance (no PBC)           | -min_dist_pbc      Min distance with PBC              |"
    echo " | -filter_dist  Filter by min_dist (no PBC)     | -filter_dist_pbc   Filter by min_dist (PBC)           |"
    echo " | -pda          Probability density analysis    | -filter_box        Filter by box-edge length          |"
    echo " | -pynep        Deprecated PyNEP sampling       | -nep_modifier      Modify NEP model interactively     |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | Python option help: gpumdkit.sh -<option> -h    Plot list: gpumdkit.sh -plt -h                        |"
    echo " +-------------------------------------------------------------------------------------------------------+"
}

function calculator_help_table(){
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                                      CALCULATOR TOOLS                                                 |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | Usage: gpumdkit.sh -calc <type> [args...]                                                             |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | ionic-cond <element> <charge>                 Calculate ionic conductivity from MSD data              |"
    echo " | nep <input.xyz> <output.xyz> <nep_model>      Calculate energy/force/virial with a NEP model          |"
    echo " | des <input.xyz> <output.npy> <nep_model> <el> Calculate NEP descriptors for one element               |"
    echo " | doas <input.xyz> <nep_model> <output.txt>     Calculate density of atomistic states                   |"
    echo " | neb <initial.xyz> <final.xyz> <n_images> <nep> Run NEB calculation with a NEP model                   |"
    echo " | minimize <structure> <nep_model> [fmax] [n]   Minimize a structure with a NEP model                   |"
    echo " | msd <trajectory.xyz> <element> <dt_fs> [n]    Calculate MSD from an extxyz trajectory                 |"
    echo " | nlist [script args...]                        Build neighbor lists                                    |"
    echo " | disp [script args...]                         Calculate displacement from trajectory                  |"
    echo " | avg-struct [script args...]                   Calculate averaged structure                            |"
    echo " | oct-tilt [script args...]                     Calculate octahedral tilt                               |"
    echo " | pol-abo3 [script args...]                     Calculate local polarization for ABO3                   |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " | Examples: gpumdkit.sh -calc ionic-cond Li 1                                                           |"
    echo " |           gpumdkit.sh -calc msd dump.xyz Li 10                                                        |"
    echo " |           gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O                                  |"
    echo " +-------------------------------------------------------------------------------------------------------+"
}

function citation(){
echo " +------------------------------------------------------+"
echo " |           THANK YOU FOR USING GPUMDKIT               |"
echo " +------------------------------------------------------+"
echo " | If you find it useful, please cite our paper:        |"
echo " |                                                      |"
echo " | GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP  |"
echo " |           MGE Advances, 2026, 4, e70074              |"
echo " |       (https://doi.org/10.1002/mgea.70074)           |"
echo " |                                                      |"
echo " |     Welcome to join our QQ group (825696376) !       |"
echo " +------------------------------------------------------+"
}

#--------------------- command line ----------------------
# Helper: run a Python calculator script with author/path banner
# Usage: run_python_script "Author Name" /path/to/script.py [args...]
run_python_script() {
    local author="$1"; shift
    local script="$1"; shift
    echo " Calling script by ${author}. "
    echo " Code path: ${script}"
    echo " ------------------------------------------------------------"
    python "${script}" "$@"
}

if [ ! -z "$1" ]; then
    case $1 in
        -h|-help) help_info_table ;;
        -skill)
            source ${utils_path}/skill_info.sh
            skill_info_table ;;
        -doctor)
            GPUMDKIT_BASH_VERSION="${BASH_VERSION}" python "${utils_path}/doctor.py" "${@:2}" ;;
        -clean) 
            source ${utils_path}/clean_extra_files.sh
            clean_extra_files ;;
        -update|-U)
            source ${utils_path}/update_gpumdkit.sh
            update_gpumdkit
            source $GPUMDkit_path/docs/updates.info ;;            
        -time)
            case $2 in
                gpumd) bash ${analyzer_path}/time_consuming_gpumd.sh ;;
                nep|gnep) bash ${analyzer_path}/time_consuming_nep.sh ;;                
                *)
                    echo " See the codes in analyzer folder for more details"
                    echo " Code path: ${analyzer_path}/time_consuming_*.sh"
                    exit 1 ;;
            esac ;;
        -plt)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    "thermo") python ${plt_path}/plt_thermo.py $3 ;;
                    "thermo2") python ${plt_path}/plt_thermo2.py $3 ;;
                    "thermo3") python ${plt_path}/plt_thermo3.py $3 ;;                        
                    "train") python ${plt_path}/plt_train.py $3 ;;
                    "train_density") python ${plt_path}/plt_train_density.py $3 ;;                 
                    "prediction"|"test") python ${plt_path}/plt_prediction.py $3 ;; 
                    "parity_density") python ${plt_path}/plt_parity_density.py $3 ;;
                    "train_test") python ${plt_path}/plt_train_test.py $3 ;;
                    "born_charge"|"bec") python ${plt_path}/plt_born_charge.py $3 ;;
                    "msd") python ${plt_path}/plt_msd.py $3 ;;
                    "msd_all") python ${plt_path}/plt_msd_all.py $3 ${@:4} ;;
                    "msd_conv") python ${plt_path}/plt_msd_convergence_check.py $3 ;;
                    "msd_sdc") python ${plt_path}/plt_msd_sdc.py $3 ;;
                    "sdc") python ${plt_path}/plt_sdc.py $3 ;;
                    "rdf") python ${plt_path}/plt_rdf.py ${@:3} ;;
                    "vac") python ${plt_path}/plt_vac.py $3 ;;
                    "restart") python ${plt_path}/plt_nep_restart.py $3 ;;
                    "dimer") python ${plt_path}/plt_dimer.py $3 $4 $5 $6 ;;
                    "force_errors") python ${plt_path}/plt_force_errors.py $3 ;;
                    "des") python ${plt_path}/plt_descriptors.py $3 ${@:4} ;;
                    "lr") python ${plt_path}/plt_learning_rate.py $3 ;;
                    "doas") python ${plt_path}/plt_doas.py $3 $4 ;;
                    "arrhenius_d"|"D") python ${plt_path}/plt_arrhenius_d.py $3 ;;
                    "D_xyz") python ${plt_path}/plt_arrhenius_d_xyz.py $3 ;;
                    "arrhenius_sigma"|"sigma") python ${plt_path}/plt_arrhenius_sigma.py $3 ;;
                    "sigma_xyz") python ${plt_path}/plt_arrhenius_sigma_xyz.py $3 ;;
                    "net_force") python ${plt_path}/plt_net_force.py ${@:3} ;;
                    "emd") python ${plt_path}/plt_emd.py ${@:3} ;;
                    "nemd") python ${plt_path}/plt_nemd.py ${@:3} ;;
                    "hnemd") python ${plt_path}/plt_hnemd.py ${@:3} ;;
                    "pdos") python ${plt_path}/plt_pdos.py $3 ;;
                    "plane-grid") python ${plt_path}/plt_plane_grid.py ${@:3} ;;
                    "cohesive") python ${plt_path}/plt_cohesive.py ${@:3} ;;
                    "viscosity") python ${plt_path}/plt_viscosity.py ${@:3} ;;
                    "rdf_pmf") python ${plt_path}/plt_rdf_pmf.py ${@:3} ;;
                    "charge")
                        echo " +----------------------------------------------------------+"
                        echo " | Please ensure you are using full batch training process. |"
                        echo " | If not, run the prediction step before plotting to avoid |"
                        echo " | inconsistencies in the atomic order between the training |"
                        echo " | set and charge_train.out.                                |"
                        echo " +----------------------------------------------------------+"
                        python ${plt_path}/plt_charge.py $3 ;;
                    *)
                        echo " Unknown plot type: $2"
                        echo " Available types are listed below."
                        source ${GPUMDkit_path}/src/f6_plots.sh; f6_plots_two_column; exit 1 ;;
                esac
            else
                source ${GPUMDkit_path}/src/f6_plots.sh
                f6_plots_two_column
                echo " See the codes in plt_scripts for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/plt_scripts"
            fi ;;

        -calc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    ionic-cond)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_ion_conductivity.py" "${@:3}" ;;
                    neb)
                        run_python_script "Zhoulin LIU (1776627910@qq.com)" "${calc_path}/neb_calculation.py" "${@:3}" ;;
                    nep)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_properties_with_nep.py" "${@:3}" ;;
                    des)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_descriptors.py" "${@:3}" ;;
                    doas)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_doas.py" "${@:3}" ;;
                    minimize)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_minimize.py" "${@:3}" ;;
                    msd)
                        run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${calc_path}/calc_msd.py" "${@:3}" ;;
                    nlist)
                        run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${calc_path}/calc_neighbor_list.py" "${@:3}" ;;
                    disp)
                        run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${calc_path}/calc_displacement.py" "${@:3}" ;;
                    avg-struct)
                        run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${calc_path}/calc_averaged_structure.py" "${@:3}" ;;
                    oct-tilt)
                        run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${calc_path}/calc_oct_tilt.py" "${@:3}" ;;
                    pol-abo3)
                        run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${calc_path}/calc_polarization_abo3.py" "${@:3}" ;;
                    *)
                        echo " See the codes in calculators folder for more details"
                        echo " Code path: ${calc_path}"; exit 1 ;;
                esac
            else
                calculator_help_table
            fi ;;

        -range)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/energy_force_virial_analyzer.py" "${@:2}" ;;

        -replicate)
            run_python_script "Boyi SITU (situboyi@westlake.edu.cn)" "${format_conv_path}/replicate.py" "${@:2}" ;;

        -out2xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Yanzhou WANG et al. "
                bash ${format_conv_path}/out2xyz.sh $2
                echo " Code path: ${format_conv_path}/out2xyz.sh"
            else
                echo " Usage: -out2xyz dir_name (eg. gpumdkit.sh -out2xyz .)"
                echo " See the source code of out2xyz.sh for more details"
                echo " Code path: ${format_conv_path}/out2xyz.sh"
            fi ;;

        -out2exyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/out2exyz.py" "${@:2}" ;;

        -cp2k2xyz)
            echo " Calling script by Chen HUA "
            python ${format_conv_path}/cp2k_log2xyz.py
            echo " Code path: ${format_conv_path}/cp2k_log2xyz.py" ;;

        -pos2exyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/pos2exyz.py" "${@:2}" ;;

        -cif2pos)
            run_python_script "Boyi SITU (situboyi@westlake.edu.cn)" "${format_conv_path}/cif2pos.py" "${@:2}" ;;

        -cif2exyz)
            run_python_script "Boyi SITU (situboyi@westlake.edu.cn)" "${format_conv_path}/cif2exyz.py" "${@:2}" ;;

        -exyz2pos)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/exyz2pos.py" "${@:2}" ;;

        -xdat2exyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/xdatcar2exyz.py" "${@:2}" ;;

        -pos2lmp)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/pos2lmp.py" "${@:2}" ;;

        -lmp2exyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/lmp2exyz.py" "${@:2}" ;;

        -traj2exyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/traj2exyz.py" "${@:2}" ;;

        -dp2xyz)
            run_python_script "Denan LI (lidenan@westlake.edu.cn)" "${format_conv_path}/dp2xyz.py" "${@:2}" ;;

        -addgroup|-addlabel)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/add_groups.py" "${@:2}" ;;

        -addweight)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/add_weight.py" "${@:2}" ;;

        -get_frame)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/get_frame.py" "${@:2}" ;;

        -clean_xyz)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${format_conv_path}/clean_xyz.py" "${@:2}" ;;

        -min_dist)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/get_min_dist.py" "${@:2}" ;;

        -min_dist_pbc)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/get_min_dist_pbc.py" "${@:2}" ;;

        -filter_dist)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/filter_structures_by_distance.py" "${@:2}" ;;

        -filter_dist_pbc)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/filter_structures_by_distance_pbc.py" "${@:2}" ;;

        -filter_box)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/filter_exyz_by_box.py" "${@:2}" ;;

        -filter_value)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/filter_exyz_by_value.py" "${@:2}" ;;

        -filter_range)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/filter_dist_range.py" "${@:2}" ;;

        -analyze_comp)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/analyze_composition.py" "${@:2}" ;;

        -get_volume)
            python ${analyzer_path}/get_volume.py ;;

        -chem_species)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/analyze_chem_species.py" "${@:2}" ;;

        -pynep)
            source ${GPUMDkit_path}/src/f2_sample_structures.sh
            parallel_pynep_sample_structures ;;

        -nep_modifier)
            python ${GPUMDkit_path}/Scripts/utils/nep_modifier/nep_modifier.py ;;

        -frame_range)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${sample_path}/frame_range.py" "${@:2}" ;;

        -re_atoms)
            run_python_script "Dian HUANG (huangdian@stu.xjtu.edu.cn)" "${utils_path}/renumber_atoms.py" "${@:2}" ;;

        -cbc)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/charge_balance_check.py" "${@:2}" ;;

        -pda)
            run_python_script "Zihan YAN (yanzihan@westlake.edu.cn)" "${analyzer_path}/probability_density_analysis.py" "${@:2}" ;;

        *)
            # echo " Unknown option: $1 "; help_info_table; exit 1 ;;
            alias_key="${1#-}"
            CUSTOM_CONFIG="$HOME/.gpumdkit.in"
            if [ ! -f "$CUSTOM_CONFIG" ]; then
                echo "Unknown: $1"
                echo "Use -h for help, or create $CUSTOM_CONFIG to define custom commands"
                exit 1
            fi
            func_name="custom_${alias_key}"
            source "$CUSTOM_CONFIG"
            if ! type "$func_name" >/dev/null 2>&1; then
                echo "Unknown: $1"
                echo "Tip: Define the function $func_name() { ... } in $CUSTOM_CONFIG to use this alias."
                help_info_table; exit 1
            fi
            shift
            "$func_name" "$@";;
    esac
    exit
fi

## logo
echo -e "\
           ____ ____  _   _ __  __ ____  _    _ _   
          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_ 
         | |  _| |_) | | | | |\/| | | | | |/ / | __|
         | |_| |  __/| |_| | |  | | |_| |   <| | |_ 
          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|
                                           
          GPUMDkit Version ${VERSION}
    Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
 Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA
      "
menu
main
citation
