#--------------------- function 2 sample structures ----------------------
# These functions are used to sample structures

# Sample structures from extxyz (uniform or random)
function f201_sample_structures(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/sample_structures |"
echo " | Script: sample_structures.py                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <extxyz_file> <sampling_method> <num_samples> [skip_num]"
echo " [skip_num]: number of initial frames to skip, default value is 0"
echo " Sampling_method: 'uniform' or 'random'"
echo " Examp: train.xyz uniform 50 "
echo " ------------>>"
read -p " " sample_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py ${sample_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py"
echo " ---------------------------------------------------"
}

# Sample structures by pynep
function f202_pynep_sample_structures(){
echo " +-------------------------------------------------+"
echo " |  Function 202 (pynep sampling) is deprecated.   |"
echo " |  PyNEP package is no longer actively maintained.|"
echo " |  Please use 203) neptrain sampling instead.     |"
echo " |  If you still insist on pynep, run directly:    |"
echo " |                 gpumdkit.sh -pynep              |"
echo " +-------------------------------------------------+"
exit 0
}

# Parallel sample structures by pynep
function parallel_pynep_sample_structures(){
# YELLOW="\033[1;33m"; RESET="\033[0m"
# echo -e "${YELLOW}"
echo " +-------------------------------------------------+"
echo " | PyNEP package is no longer actively maintained  |"
echo " |       Recommend using function 203 instead      |"
echo " |  Moreover, qNEP model is not supported for now. |"
echo " +-------------------------------------------------+"
# echo -e "${RESET}"
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/sample_structures |"
echo " | Script: parallel_pynep_select_structs.py        |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <sample.xyz> <train.xyz> <nep_model>"
echo " Examp: dump.xyz train.xyz nep.txt [threads]"
echo " [threads]: number of CPU threads to use, default value is cpu_count"
echo " ------------>>"
read -p " " sample_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/parallel_pynep_select_structs.py ${sample_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/parallel_pynep_select_structs.py"
echo " ---------------------------------------------------"
}

# Sample structures by neptrain
function f203_neptrain_sample_structures(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/sample_structures |"
echo " | Script: neptrain_select_structs.py              |"
echo " | Developer: Benrui TANG (tang070205@proton.me)   |"
echo " >-------------------------------------------------<"
echo " Input <sample.xyz> <train.xyz> <nep_model>"
echo " Examp: dump.xyz train.xyz nep.txt "
echo " ------------>>"
read -p " " sample_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py ${sample_choice}
rm -f dpdispatcher.log
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py"
echo " ---------------------------------------------------"
echo " "
echo " +--------------------------------------------------------+"
echo " | For this function, please cite the NepTrain paper:     |"
echo " | Chen et al. Comput. Phys. Commun., 2025, 317, 109859.  |"
echo " | URL: https://doi.org/10.1016/j.cpc.2025.109859         |"
echo " +--------------------------------------------------------+"
}

# Perturb structure
function f204_perturb_structure(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/sample_structures |"
echo " | Script: perturb_structure.py                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>"
echo " The default paramters for perturb are 20 0.03 0.2 uniform"
echo " Examp: POSCAR 20 0.03 0.2 uniform"
echo " ------------>>"
read -p " " perturb_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py ${perturb_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py"
echo " ---------------------------------------------------"
}

# Select max force deviation structs
function f205_select_max_force_deviation_structs(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/sample_structures |"
echo " | Script: select_max_modev.py                     |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " +----------------------------------------------------+"
echo " | Select max force deviation structs from active.xyz |"
echo " |     generated by the active command in gpumd.      |"
echo " +----------------------------------------------------+"
echo " Input <structs_num> <threshold> (eg. 200 0.15)"
echo " ------------>>"
read -p " " modev_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py ${modev_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py"
echo " ---------------------------------------------------"
}

#--------------------- function 2 ----------------------
function f2_sample_structures(){
echo " +------------------------------------------------------+"
echo " |                 SAMPLE STRUCTURE TOOLS               |"
echo " +------------------------------------------------------+"
echo " | 201) Sample structures from extxyz                   |"
echo " | 202) FPS sampling by PyNEP [deprecated]              |"
echo " | 203) FPS sampling by NepTrain [preferred]            |"
echo " | 204) Perturb structure                               |"
echo " | 205) Select max force deviation structs              |"
echo " +------------------------------------------------------+"
echo " | 000) Return to the main menu                         |"
echo " +------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "201" "202" "203" "204" "205")
read -p " " num_choice
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "201") f201_sample_structures ;;
    "202") f202_pynep_sample_structures ;;
    "203") f203_neptrain_sample_structures ;;
    "204") f204_perturb_structure ;;
    "205") f205_select_max_force_deviation_structs ;;
    "000") menu; main ;;
esac
}
