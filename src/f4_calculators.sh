# ============================================================
# GPUMDkit calculators module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit
#           for GPUMD and NEP, MGE Advances, 2026, 4, e70074
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

# Calculate ionic conductivity
function f401_calc_ionic_conductivity(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_ion_conductivity.py                |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <element> <charge> (eg. Li 1)"
echo " ------------>>"
read -p " " ion_cond_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py ${ion_cond_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py"
echo " ---------------------------------------------------"
}

# Calculate properties with NEP
function f402_calc_properties_with_nep(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_properties_with_nep.py             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <output.xyz> <nep_model> "
echo " Example: input.xyz output.xyz nep.txt"
echo " ------------>>"
read -p " " input_calc_properties
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py ${input_calc_properties}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py"
echo " ---------------------------------------------------"
}

# Calculate descriptors
function f403_calc_descriptors(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_descriptors.py                     |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <output.npy> <nep_model> <element>"
echo " Example: train.xyz des_Li.npy nep.txt Li"
echo " ------------>>"
read -p " " input_calc_descriptors
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_descriptors.py ${input_calc_descriptors}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_descriptors.py"
echo " ---------------------------------------------------"
}

# Calculate density of atomistic states (DOAS)
function f404_calc_doas(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_doas.py                            |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <nep_model> <output_file>"
echo " Example: dump.xyz nep.txt doas.out"
echo " ------------>>"
read -p " " input_calc_doas
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_doas.py ${input_calc_doas}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_doas.py"
echo " ---------------------------------------------------"
}

# Calculate nudged elastic band (NEB) by nep
function f405_calc_neb(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: neb_calculation.py                      |"
echo " | Developer: Zhoulin LIU (1776627910@qq.com)      |"
echo " >-------------------------------------------------<"
echo " Input <initial_struct> <final_struct> <n_image> <nep_model>"
echo " Example: IS.xyz FS.xyz 5 nep.txt"
echo " ------------>>"
read -p " " input_calc_neb
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py ${input_calc_neb}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py"
echo " ---------------------------------------------------"
}

# Build neighbor list for displacement analysis
function f406_calc_neighbor_list(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_neighbor_list.py                   |"
echo " | Developer: Denan LI (lidenan@westlake.edu.cn)   |"
echo " >-------------------------------------------------<"
echo " Input script arguments (eg. -c 4 -n 12 -C Pb Sr -E O)"
echo " See https://gpumdkit.cn/htmls/polar_material_analysis.html for details"
echo " ------------>>"
read -p " " input_calc_neighbor_list
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_neighbor_list.py ${input_calc_neighbor_list}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_neighbor_list.py"
echo " ---------------------------------------------------"
}

# Calculate displacement from trajectory and neighbor list
function f407_calc_displacement(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_displacement.py                    |"
echo " | Developer: Denan LI (lidenan@westlake.edu.cn)   |"
echo " >-------------------------------------------------<"
echo " Input script arguments (eg. -i movie.xyz -n nl-Pb-O.dat -o displacements.dat)"
echo " See https://gpumdkit.cn/htmls/polar_material_analysis.html for details"
echo " ------------>>"
read -p " " input_calc_displacement
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_displacement.py ${input_calc_displacement}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_displacement.py"
echo " ---------------------------------------------------"
}

# Calculate averaged structure from trajectory
function f408_calc_averaged_structure(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_averaged_structure.py              |"
echo " | Developer: Denan LI (lidenan@westlake.edu.cn)   |"
echo " >-------------------------------------------------<"
echo " Input script arguments (eg. -i movie.xyz -l 0.2 -o averaged_structure.xyz)"
echo " See https://gpumdkit.cn/htmls/polar_material_analysis.html for details"
echo " ------------>>"
read -p " " input_calc_averaged_structure
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_averaged_structure.py ${input_calc_averaged_structure}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_averaged_structure.py"
echo " ---------------------------------------------------"
}

# Calculate octahedral tilt from trajectory and B-O neighbor list
function f409_calc_oct_tilt(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_oct_tilt.py                        |"
echo " | Developer: Denan LI (lidenan@westlake.edu.cn)   |"
echo " >-------------------------------------------------<"
echo " Input script arguments (eg. -i model.xyz -n nl-Ti-O.dat -o octahedral_tilt.dat)"
echo " See https://gpumdkit.cn/htmls/polar_material_analysis.html for details"
echo " ------------>>"
read -p " " input_calc_oct_tilt
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_oct_tilt.py ${input_calc_oct_tilt}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_oct_tilt.py"
echo " ---------------------------------------------------"
}

# Calculate local polarization for ABO3 from neighbor lists
function f410_calc_polarization_abo3(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_polarization_abo3.py               |"
echo " | Developer: Denan LI (lidenan@westlake.edu.cn)   |"
echo " >-------------------------------------------------<"
echo " Input script arguments (eg. --nl-ba nl-Ti-A.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Sr=2.0 Ti=4.0 O=-2.0)"
echo " See https://gpumdkit.cn/htmls/polar_material_analysis.html for details"
echo " ------------>>"
read -p " " input_calc_polarization_abo3
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_polarization_abo3.py ${input_calc_polarization_abo3}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_polarization_abo3.py"
echo " ---------------------------------------------------"
}

function f411_minimize_structure_by_nep(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_minimize.py                         |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <structure_file> <nep.txt> [fmax=0.01] [max_steps=1000]"
echo " Supported structure formats: POSCAR, .xyz"
echo " Optional arguments: "
echo "  fmax: Force convergence threshold in eV/Ang (default: 0.01)"
echo "  max_steps: Maximum number of optimization steps (default: 1000)"
echo " Example: POSCAR nep.txt 0.01 1000"
echo " ------------>>"
read -p " " input_minimize_structure
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_minimize.py ${input_minimize_structure}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_minimize.py"
echo " ---------------------------------------------------"
}

function f412_calc_msd_from_trajectory(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_msd.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]"
echo "   Optional argument: max_corr_steps (default: frame number)"
echo " Example: dump.xyz Li 10"
echo " ------------>>"
read -p " " input_calc_msd
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_msd.py ${input_calc_msd}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_msd.py"
echo " ---------------------------------------------------"
}

# main function of calculators
function f4_calculators(){
echo " +----------------------------------------------------------+"
echo " |                     CALCULATOR TOOLS                     |"
echo " +----------------------------------------------------------+"
echo " | 401) Calc ionic conductivity                             |"
echo " | 402) Calc properties by nep                              |"
echo " | 403) Calc descriptors of specific elements               |"
echo " | 404) Calc density of atomistic states (DOAS)             |"
echo " | 405) Calc nudged elastic band (NEB) by nep               |"
echo " | 406) Build neighbor list                                 |"
echo " | 407) Calc displacement from trajectory                   |"
echo " | 408) Calc averaged structure                             |"
echo " | 409) Calc octahedral tilt                                |"
echo " | 410) Calc polarization for ABO3                          |"
echo " | 411) Minimize structure by nep                           |"
echo " | 412) Calc mean square displacement (MSD) from trajectory |"
echo " +----------------------------------------------------------+"
echo " | 000) Return to the main menu                             |"
echo " +----------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "401" "402" "403" "404" "405" "406" "407" "408" "409" "410" "411" "412")
read -p " " num_choice
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
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
    "000") menu; main ;;
esac
}
