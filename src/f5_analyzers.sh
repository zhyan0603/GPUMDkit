#--------------------- function 5 analyzer ----------------------
function f501_analyze_composition(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: analyze_composition.py                  |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> you want to analyze"
echo " Examp: train.xyz"
echo " ------------>>"
read -p " " input_extxyz
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/analyze_composition.py ${input_extxyz}
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/analyze_composition.py"
echo " ---------------------------------------------------"
}

function f502_find_outliers(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in analyzer      |"
echo " | Script: find_outliers.py                        |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input the threshold of RMSE to identify outliers"
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py
echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py"
echo " ---------------------------------------------------"
}

function f5_analyzers(){
echo " ------------>>"
echo " 501) Analyze composition of extxyz"
echo " 502) Find outliers of extxyz"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "501" "502") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "501") f501_analyze_composition ;;            
    "502") f502_find_outliers ;;
    "000") menu; main ;;
esac
}