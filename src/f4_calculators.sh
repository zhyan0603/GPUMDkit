#--------------------- function 4 calculators ----------------------
# These functions are used to do the calculators
# See the source codes in Scripts/calculators for more details

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
echo " Examp: input.xyz output.xyz nep.txt"
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
echo " Examp: train.xyz des_Li.npy nep.txt Li"
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
echo " Examp: dump.xyz nep.txt doas.out"
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
echo " Examp: IS.xyz FS.xyz 5 nep.txt"
echo " ------------>>"
read -p " " input_calc_neb
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py ${input_calc_neb}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py"
echo " ---------------------------------------------------"
}

# main function of calculators
function f4_calculators(){
echo " ------------>>"
echo " 401) Calc ionic conductivity"
echo " 402) Calc properties by nep"
echo " 403) Calc descriptors of specific elements"
echo " 404) Calc density of atomistic states (DOAS)"
echo " 405) Calc nudged elastic band (NEB) by nep"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "401" "402" "403" "404" "405") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
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
    "000") menu; main ;;
esac
}