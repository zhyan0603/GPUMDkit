# ============================================================
# GPUMDkit workflow module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit
#           for GPUMD and NEP, MGE Advances, 2026, 4, e70074
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

# CP2K scf batch pretreatment
function cp2k_batch_pretreatment(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/workflow          |"
echo " | Script: scf_batch_pretreatment_cp2k.py          |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " >-------------------- NOTICE  --------------------<"
echo " | Please ensure your cp2k.inp file is configured  |"
echo " | to read coordinates from 'pos.xyz'.             |"
echo " | Refer to Scripts/workflow/cp2k_template.inp     |"
echo " >-------------------------------------------------<"

if ! IFS= read -r -p " Are you ready to continue? [y/N]: " confirm; then
    echo " Input closed. Exiting."
    return 1
fi
if [[ "$confirm" != [yY]* ]]; then
    echo " Operation cancelled. Exiting..."
    return
fi

echo " ---------------------------------------------------"
echo " Input <extxyz_file> <template.inp> <prefix_name>"
echo " Example: dump.xyz template.inp H2O"
echo " ------------>>"
read_menu_choice input_results || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_cp2k.py ${input_results}
echo " Code path: ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_cp2k.py"
echo " ---------------------------------------------------"
}

# SCF batch pretreatment (main function)
function f301_scf_batch_pretreatment(){ 
    echo " +-------------------------------------------------------------+"
    echo " |                 SCF BATCH PRETREATMENT TOOLS                |"
    echo " +-------------------------------------------------------------+"
    echo " | 1) VASP SCF batch pretreatment                              |"
    echo " | 2) CP2K SCF batch pretreatment                              |"
    echo " +-------------------------------------------------------------+"
    echo " Input the function number:"
    valid_menu_choices=("1" "2")
    read_menu_choice num_choice || return 1
    while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
    do
      echo " ------------>>"
      echo " Please reinput function number..."
      read_menu_choice num_choice || return 1
    done   
    case $num_choice in
        "1")
            source ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_vasp.sh
            vasp_scf_batch_pretreatment
            ;;
        "2")
            cp2k_batch_pretreatment
            ;;          
    esac
}

# main function of workflow
function f3_workflow_dev(){
echo " +---------------------------------------------------------+"
echo " |                      WORKFLOW TOOLS                     |"
echo " +---------------------------------------------------------+"
echo " | 301) SCF batch pretreatment                             |"
echo " | 302) MD sample batch pretreatment (gpumd)               |"
echo " | 303) MD sample batch pretreatment (lmp)                 |"
echo " +---------------------------------------------------------+"
echo " | 000) Return to the main menu                            |"
echo " +---------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "301" "302" "303")
read_menu_choice num_choice || return 1
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read_menu_choice num_choice || return 1
done

case $num_choice in
    "301")
        f301_scf_batch_pretreatment
        ;;
    "302")
        source ${GPUMDkit_path}/Scripts/workflow/md_sample_batch_pretreatment_gpumd.sh
        f302_md_sample_batch_pretreatment_gpumd
        ;;
    "303")
        source ${GPUMDkit_path}/Scripts/workflow/md_sample_batch_pretreatment_lmp.sh
        f303_md_sample_batch_pretreatment_lmp
        ;;          
    "000")
        menu
        main
        ;;
esac
}
