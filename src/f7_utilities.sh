#--------------------- function 7 utilities ----------------------
# Miscellaneous tools that do not fit other modules

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

function f7_utilities(){
echo " +------------------------------------------------------+"
echo " |                    UTILITY TOOLS                     |"
echo " +------------------------------------------------------+"
echo " | 703) Time consuming analyzer (gpumd/nep)             |"
echo " +------------------------------------------------------+"
echo " | 000) Return to the main menu                         |"
echo " +------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "703")
read -p " " num_choice
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "703") f703_time_consuming_analyzer ;;
    "000") menu; main ;;
esac
}
