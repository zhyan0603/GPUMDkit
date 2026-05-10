#--------------------- function 7 utilities ----------------------
# Miscellaneous tools that do not fit other modules

function f701_clean_workdir(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/utils             |"
echo " | Script: clean_extra_files.sh                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " ---------------------------------------------------"
source ${utils_path}/clean_extra_files.sh
clean_extra_files
echo " Code path: ${utils_path}/clean_extra_files.sh"
echo " ---------------------------------------------------"
}

function f702_update_gpumdkit(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/utils             |"
echo " | Script: update_gpumdkit.sh                      |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " ---------------------------------------------------"
source ${utils_path}/update_gpumdkit.sh
update_gpumdkit
if [ -f "${GPUMDkit_path}/docs/updates.info" ]; then
    source ${GPUMDkit_path}/docs/updates.info
fi
echo " Code path: ${utils_path}/update_gpumdkit.sh"
echo " ---------------------------------------------------"
}

function f703_time_consuming_analyzer(){
echo " >-------------------------------------------------<"
echo " | Calling scripts in Scripts/analyzer             |"
echo " | Script: time_consuming_gpumd.sh / nep.sh        |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input mode: gpumd or nep"
echo " Example: gpumd"
echo " ------------>>"
read -p " " time_mode
echo " ---------------------------------------------------"
case $time_mode in
    gpumd) bash ${analyzer_path}/time_consuming_gpumd.sh ;;
    nep|gnep) bash ${analyzer_path}/time_consuming_nep.sh ;;
    *)
        echo " Please input 'gpumd' or 'nep'." ;;
esac
echo " Code path: ${analyzer_path}/time_consuming_*.sh"
echo " ---------------------------------------------------"
}

function f704_renumber_atoms(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/utils             |"
echo " | Script: renumber_atoms.py                       |"
echo " | Developer: Dian HUANG (huangdian@stu.xjtu.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input_file> <output_file>"
echo " Example: dump.in dump_renum.in"
echo " ------------>>"
read -p " " input_file output_file
echo " ---------------------------------------------------"
python ${utils_path}/renumber_atoms.py ${input_file} ${output_file}
echo " Code path: ${utils_path}/renumber_atoms.py"
echo " ---------------------------------------------------"
}

function f7_utilities(){
echo " +------------------------------------------------------+"
echo " |                    UTILITY TOOLS                     |"
echo " +------------------------------------------------------+"
echo " | 701) Clean extra files in current work directory     |"
echo " | 702) Update GPUMDkit from GitHub                     |"
echo " | 703) Time consuming analyzer (gpumd/nep)             |"
echo " | 704) Renumber atom IDs in LAMMPS dump file           |"
echo " +------------------------------------------------------+"
echo " | 000) Return to the main menu                         |"
echo " +------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "701" "702" "703" "704")
read -p " " num_choice
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "701") f701_clean_workdir ;;
    "702") f702_update_gpumdkit ;;
    "703") f703_time_consuming_analyzer ;;
    "704") f704_renumber_atoms ;;
    "000") menu; main ;;
esac
}
