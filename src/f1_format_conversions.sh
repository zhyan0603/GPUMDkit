#--------------------- function 1 format conversion ----------------------
# These functions are used to convert the format of the files

function f101_out2xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: out2xyz.sh                              |"
echo " | Developer: Yanzhou WANG (yanzhowang@gmail.com)  |"
echo " >-------------------------------------------------<"
echo " Choose the type of conversion:"
echo " 1) OUTCAR to extxyz"
echo " 2) vasprun.xml to extxyz"
echo " ------------>>"
read -p " " outcar_type
if [ "$outcar_type" == "1" ]; then
    echo " Input the directory containing OUTCARs"
    echo " Examp: ./ "
    echo " ------------>>"
    read -p " " dir_outcars
    echo " ---------------------------------------------------"
    bash ${GPUMDkit_path}/Scripts/format_conversion/out2xyz.sh ${dir_outcars}
    echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/out2xyz.sh"
    echo " ---------------------------------------------------"
elif [ "$outcar_type" == "2" ]; then
    echo " Input the directory containing vasprun.xml files"
    echo " Examp: ./ "
    echo " ------------>>"
    read -p " " dir_vasprun
    echo " ---------------------------------------------------"
    python ${GPUMDkit_path}/Scripts/format_conversion/xml2xyz.py ${dir_vasprun}
    echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/xml2xyz.py"
    echo " ---------------------------------------------------"
else
    echo " Your input is illegal, please try again"
    return
fi
}

function f102_mtp2xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: mtp2xyz.py                              |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>"
echo " Examp: train.cfg Pd Ag"
echo " ------------>>"
read -p " " mtp_variables
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/mtp2xyz.py ${mtp_variables}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/mtp2xyz.py"
echo " ---------------------------------------------------"
}

function f103_cp2k2xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cp2k2xyz.py                             |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input [pos.xyz] [frc.xyz] [cell.cell] [-shifted yes/no] "
echo " ------------>>"
read -p " " cp2k_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/cp2k2xyz.py ${cp2k_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/cp2k2xyz.py"
echo " ---------------------------------------------------"
}

function f104_abacus2xyz(){ 
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: abacus2xyz-scf.sh/abacus2xyz-md.sh      |"
echo " | Developer: Benrui TANG (tang070205@proton.me)   |"
echo " >-------------------------------------------------<"
echo " Choose the type of ABACUS calculation:"
echo " 1) SCF calculation"
echo " 2) MD calculation"
echo " ------------>>"
read -p " " abacus_type
if [ "$abacus_type" == "1" ]; then
    echo " Input the directory containing running_scf.log"
    echo " Examp: ./ "
    echo " ------------>>"
    read -p " " dir_abacus_scf
    echo " ---------------------------------------------------"
    bash ${GPUMDkit_path}/Scripts/format_conversion/abacus2xyz_scf.sh ${dir_abacus_scf}
    echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/abacus2xyz_scf.sh"
    echo " ---------------------------------------------------"
elif [ "$abacus_type" == "2" ]; then
    echo " Input the directory containing running_md.log and MD_dump file"
    echo " ------------>>"
    read -p " " dir_abacus_md
    echo " ---------------------------------------------------"
    bash ${GPUMDkit_path}/Scripts/format_conversion/abacus2xyz_md.sh ${dir_abacus_md}
    echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/abacus2xyz_md.sh"
    echo " ---------------------------------------------------"
else
    echo " Your input is illegal, please try again"
    return
fi
}

function f105_extxyz2poscar(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: exyz2pos.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input the name of extxyz"
echo " Examp: ./train.xyz "
echo " ------------>>"
read -p " " filename
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py ${filename}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py"
echo " ---------------------------------------------------"
}

function f1_format_conversion(){
echo " ------------>>"
echo " 101) Convert VASP to extxyz"
echo " 102) Convert mtp to extxyz"
echo " 103) Convert CP2K to extxyz"
echo " 104) Convert ABACUS to extxyz"
echo " 105) Convert extxyz to POSCAR"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "101" "102" "103" "104" "105") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "101") f101_out2xyz ;;
    "102") f102_mtp2xyz ;;
    "103") f103_cp2k2xyz ;;
    "104") f104_abacus2xyz ;;
    "105") f105_extxyz2poscar ;;       
    "000") menu; main ;;
esac
}