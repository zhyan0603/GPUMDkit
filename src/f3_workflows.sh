#--------------------- function 3 workflow ----------------------
# These functions are used to do the workflow
# See the source codes in Scripts/workflow for more details

# CP2K scf batch pretreatment
function cp2k_batch_pretreatment(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/workflow          |"
echo " | Script: scf_batch_pretreatment_cp2k.py          |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <extxyz_file> <template.inp> <prefix_name>"
echo " Examp: dump.xyz template.inp H2O"
echo " ------------>>"
read -p " " input_results
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_cp2k.py ${input_results}
echo " Code path: ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_cp2k.py"
echo " ---------------------------------------------------"
}

# SCF batch pretreatment (main function)
function f301_scf_batch_pretreatment(){ 
    echo " ------------>>"
    echo " 1) VASP scf batch pretreatment"
    echo " 2) CP2K scf batch pretreatment"
    echo " ------------>>"
    echo " Input the function number:"
    arry_num_choice=("1" "2") 
    read -p " " num_choice
    while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
    do
      echo " ------------>>"
      echo " Please reinput function number..."
      read -p " " num_choice
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
echo " ------------>>"
echo " 301) SCF batch pretreatment"
echo " 302) MD sample batch pretreatment (gpumd)"
echo " 303) MD sample batch pretreatment (lmp)"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "301" "302" "303") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
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