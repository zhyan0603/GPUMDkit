#--------------------- function 3 workflow ----------------------
# These functions are used to do the workflow
# See the source codes in Scripts/workflow for more details

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
        source ${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment.sh
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