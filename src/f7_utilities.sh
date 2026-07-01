# ============================================================
# GPUMDkit utilities module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit
#           for GPUMD and NEP, MGE Advances, 2026, 4, e70074
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

function f701_time_consuming_analyzer(){
echo " >-------------------------------------------------<"
echo " | Calling scripts in Scripts/analyzer             |"
echo " | Script: time_consuming_gpumd.sh/nep.sh          |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input mode: gpumd, nep or gnep"
echo " Example: gpumd"
echo " ------------>>"
read_menu_choice time_mode || return 1
echo " ---------------------------------------------------"
case $time_mode in
    gpumd) bash ${analyzer_path}/time_consuming_gpumd.sh ;;
    nep|gnep) bash ${analyzer_path}/time_consuming_nep.sh ;;
    *)
        echo " Please input 'gpumd', 'nep' or 'gnep'." ;;
esac
echo " Code path: ${analyzer_path}/time_consuming_*.sh"
echo " ---------------------------------------------------"
}

function f7_utilities(){
echo " +------------------------------------------------------+"
echo " |                    UTILITY TOOLS                     |"
echo " +------------------------------------------------------+"
echo " | 701) Time consuming analyzer (gpumd/nep)             |"
echo " +------------------------------------------------------+"
echo " | 000) Return to the main menu                         |"
echo " +------------------------------------------------------+"
echo " Input the function number:"

valid_menu_choices=("000" "701")
read_menu_choice num_choice || return 1
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read_menu_choice num_choice || return 1
done

case $num_choice in
    "701") f701_time_consuming_analyzer ;;
    "000") menu; main ;;
esac
}
