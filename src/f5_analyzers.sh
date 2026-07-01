# ============================================================
# GPUMDkit analyzer module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit
#           for GPUMD and NEP, MGE Advances, 2026, 4, e70074
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

# Analyze composition of extxyz
function f501_analyze_composition(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: analyze_composition.py                  |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> you want to analyze"
echo " Example: train.xyz"
echo " ------------>>"
read_menu_choice input_extxyz || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/analyze_composition.py ${input_extxyz}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/analyze_composition.py"
echo " ---------------------------------------------------"
}

# Find outliers of extxyz
function f502_find_outliers(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: find_outliers.py                        |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input the threshold of RMSE to identify outliers"
echo " Required files: train.xyz, energy_train.out, force_train.out, stress_train.out"
echo " The script will ask for energy (meV/atom), force (meV/Angstrom), and stress (GPa) thresholds."
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py"
echo " ---------------------------------------------------"
}

# Analyze chemical species of extxyz
function f503_analyze_chem_species(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: analyze_chem_species.py                 |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> you want to analyze"
echo " Example: train.xyz"
echo " ------------>>"
read_menu_choice input_extxyz || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/analyze_chem_species.py ${input_extxyz}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/analyze_chem_species.py"
echo " ---------------------------------------------------"
}

# Check charge balance of extxyz
function f504_charge_balance_check(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: charge_balance_check.py                 |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> you want to analyze"
echo " Example: train.xyz"
echo " ------------>>"
read_menu_choice input_extxyz || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/charge_balance_check.py ${input_extxyz}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/charge_balance_check.py"
echo " ---------------------------------------------------"
}

# Analyze energy/force/virial range
function f505_energy_force_virial_analyzer(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: energy_force_virial_analyzer.py         |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> and property (energy/force/virial)"
echo " Example: train.xyz force"
echo " ------------>>"
read_menu_choice analyzer_args || return 1
read -r input_extxyz property_name <<< "${analyzer_args}"
echo " Plot histogram? (y/n)"
echo " ------------>>"
read_menu_choice if_hist || return 1
echo " ---------------------------------------------------"
if [ "$if_hist" == "y" ] || [ "$if_hist" == "Y" ]; then
  python ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py ${input_extxyz} ${property_name} hist
else
  python ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py ${input_extxyz} ${property_name}
fi
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py"
echo " ---------------------------------------------------"
}

# Filter structures by minimum distance (without PBC)
function f506_filter_structures_by_distance(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: filter_structures_by_distance.py        |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> and <min_dist>"
echo " Example: dump.xyz 1.5"
echo " ------------>>"
read_menu_choice filter_args || return 1
read -r input_extxyz min_dist <<< "${filter_args}"
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance.py ${input_extxyz} ${min_dist}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance.py"
echo " ---------------------------------------------------"
}

# Get minimum interatomic distance with optional PBC
function f507_get_min_dist(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: get_min_dist.py / get_min_dist_pbc.py   |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> you want to analyze"
echo " Example: dump.xyz"
echo " ------------>>"
read_menu_choice input_extxyz || return 1
echo " Consider PBC? (y/n)"
echo " ------------>>"
read_menu_choice consider_pbc || return 1
echo " ---------------------------------------------------"
if [ "$consider_pbc" == "y" ] || [ "$consider_pbc" == "Y" ]; then
  python ${GPUMDkit_path}/Scripts/analyzer/get_min_dist_pbc.py ${input_extxyz}
  echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist_pbc.py"
else
  python ${GPUMDkit_path}/Scripts/analyzer/get_min_dist.py ${input_extxyz}
  echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist.py"
fi
echo " ---------------------------------------------------"
}

# Probability density analysis for diffusion channels
function f508_probability_density_analysis(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: probability_density_analysis.py         |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <ref_struct> <trajectory_file> <species> <interval>"
echo " Example: LLZO.vasp dump.xyz Li 0.25"
echo " ------------>>"
read_menu_choice pda_args || return 1
read -r ref_struct trajectory_file species interval <<< "${pda_args}"
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/probability_density_analysis.py ${ref_struct} ${trajectory_file} ${species} ${interval}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/probability_density_analysis.py"
echo " ---------------------------------------------------"
}

# main function of analyzers
function f5_analyzers(){
echo " +------------------------------------------------------+"
echo " |                    ANALYZER TOOLS                    |"
echo " +------------------------------------------------------+"
echo " | 501) Analyze composition of extxyz                   |"
echo " | 502) Find outliers of extxyz                         |"
echo " | 503) Analyze chemical species of extxyz              |"
echo " | 504) Check charge balance of extxyz                  |"
echo " | 505) Analyze energy/force/virial range               |"
echo " | 506) Filter structures by minimum distance           |"
echo " | 507) Get minimum interatomic distance                |"
echo " | 508) Probability density analysis                    |"
echo " +------------------------------------------------------+"
echo " | 000) Return to the main menu                         |"
echo " +------------------------------------------------------+"
echo " Input the function number:"

arry_num_choice=("000" "501" "502" "503" "504" "505" "506" "507" "508") 
read_menu_choice num_choice || return 1
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read_menu_choice num_choice || return 1
done

case $num_choice in
    "501") f501_analyze_composition ;;            
    "502") f502_find_outliers ;;
    "503") f503_analyze_chem_species ;;
    "504") f504_charge_balance_check ;;
    "505") f505_energy_force_virial_analyzer ;;
    "506") f506_filter_structures_by_distance ;;
    "507") f507_get_min_dist ;;
    "508") f508_probability_density_analysis ;;
    "000") menu; main ;;
esac
}
