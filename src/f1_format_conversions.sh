# ============================================================
# GPUMDkit format conversion module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit
#           for GPUMD and NEP, MGE Advances, 2026, e70074
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

# Convert VASP OUTCAR to extxyz
function f101_out2xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: out2xyz.sh                              |"
echo " | Developer: Yanzhou WANG (yanzhowang@gmail.com)  |"
echo " >-------------------------------------------------<"
echo " Input the directory containing OUTCARs"
echo " Example: ./ "
echo " ------------>>"
read -p " " dir_outcars
echo " ---------------------------------------------------"
bash ${GPUMDkit_path}/Scripts/format_conversion/out2xyz.sh ${dir_outcars}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/out2xyz.sh"
echo " ---------------------------------------------------"
}

# Convert mtp to extxyz
function f102_mtp2xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: mtp2xyz.py                              |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>"
echo " Example: train.cfg Pd Ag"
echo " ------------>>"
read -p " " mtp_variables
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/mtp2xyz.py ${mtp_variables}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/mtp2xyz.py"
echo " ---------------------------------------------------"
}

# Convert CP2K to extxyz - method by Chen HUA
function cp2k2xyz_chenhua(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cp2k_log2xyz.py                         |"
echo " | Developer: Chen HUA (huachen23@mails.ucas.ac.cn)|"
echo " >-------------------------------------------------<"
echo " This script converts CP2K calculations to extxyz"
echo " Need file: .log for properties and .xyz/.inp for structure"
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/cp2k_log2xyz.py
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/cp2k_log2xyz.py"
echo " ---------------------------------------------------"
}

# Convert CP2K to extxyz - method by Ke XU
function cp2k2xyz_kexu(){
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

#Convert CP2K to extxyz (main function)
function f103_cp2k2xyz(){
echo " >-------------------------------------------------<"
echo " | There are two scripts for CP2K to extxyz:       |"
echo " | 1) cp2k_log2xyz.py by Chen HUA                  |"
echo " | 2) cp2k2xyz.py by Ke XU                         |"
echo " | You can choose either one to use.               |"
echo " >-------------------------------------------------<"
echo " Choose the script to use:"
echo " 1) cp2k_log2xyz.py (from log and inp/xyz files)"
echo " 2) cp2k2xyz.py (from pos.xyz, frc.xyz, cell.cell files)"
echo " ------------>>"
read -p " " cp2k_script_choice
if [ "$cp2k_script_choice" == "1" ]; then
    cp2k2xyz_chenhua
elif [ "$cp2k_script_choice" == "2" ]; then
    cp2k2xyz_kexu
else
    echo " Your input is illegal, please try again"
    return
fi
}

# Convert ABACUS to extxyz
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
    echo " Example: ./ "
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

# Convert extxyz to POSCAR
function f105_extxyz2poscar(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: exyz2pos.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input the name of extxyz"
echo " Example: ./train.xyz "
echo " ------------>>"
read -p " " filename
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py ${filename}
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py"
echo " ---------------------------------------------------"
}

# Add group labels
function f106_add_group_labels(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: add_groups.py                           |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <POSCAR> <element1> <element2> ..."
echo " Example: POSCAR Li Y Cl"
echo " ------------>>"
read -r -a addgroup_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/add_groups.py "${addgroup_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_groups.py"
echo " ---------------------------------------------------"
}

# Add weight to extxyz
function f107_add_weight(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: add_weight.py                           |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <output.xyz> <weight>"
echo " Example: train.xyz train_weighted.xyz 5"
echo " ------------>>"
read -r -a addweight_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/add_weight.py "${addweight_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_weight.py"
echo " ---------------------------------------------------"
}

# Get a single frame from extxyz
function f108_get_frame(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: get_frame.py                            |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <frame_index>"
echo " Example: train.xyz 1"
echo " ------------>>"
read -r -a getframe_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/get_frame.py "${getframe_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/get_frame.py"
echo " ---------------------------------------------------"
}

# Clean extra info in XYZ file
function f109_clean_xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: clean_xyz.py                            |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <output.xyz>"
echo " Example: train.xyz train_clean.xyz"
echo " ------------>>"
read -r -a cleanxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/clean_xyz.py "${cleanxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/clean_xyz.py"
echo " ---------------------------------------------------"
}

# Replicate structure
function f110_replicate_structure(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: replicate.py                            |"
echo " | Developer: Boyi SITU (situboyi@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <inputfile> <outputfile> <a b c / target_num>"
echo " Example 1: POSCAR POSCAR_222.vasp 2 2 2"
echo " Example 2: POSCAR POSCAR_256.vasp 256"
echo " ------------>>"
read -r -a replicate_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/replicate.py "${replicate_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/replicate.py"
echo " ---------------------------------------------------"
}

# Convert OUTCAR to extxyz (Python version)
function out2exyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: out2exyz.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input the directory containing OUTCARs"
echo " Example: ./ "
echo " ------------>>"
read -r -p " " dir_outcars
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/out2exyz.py "${dir_outcars}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/out2exyz.py"
echo " ---------------------------------------------------"
}

# Convert POSCAR to extxyz
function pos2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: pos2exyz.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <POSCAR> <output.xyz>"
echo " Example: POSCAR model.xyz"
echo " ------------>>"
read -r -a pos2extxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py "${pos2extxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py"
echo " ---------------------------------------------------"
}

# Convert CIF to POSCAR
function cif2poscar(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cif2pos.py                              |"
echo " | Developer: Boyi SITU (situboyi@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.cif> <output.vasp>"
echo " Example: input.cif POSCAR.vasp"
echo " ------------>>"
read -r -a cif2pos_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/cif2pos.py "${cif2pos_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/cif2pos.py"
echo " ---------------------------------------------------"
}

# Convert CIF to extxyz
function cif2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cif2exyz.py                             |"
echo " | Developer: Boyi SITU (situboyi@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.cif> <output.xyz>"
echo " Example: input.cif model.xyz"
echo " ------------>>"
read -r -a cif2extxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/cif2exyz.py "${cif2extxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/cif2exyz.py"
echo " ---------------------------------------------------"
}

# Convert XDATCAR to extxyz
function xdatcar2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: xdatcar2exyz.py                         |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <XDATCAR> <output.xyz>"
echo " Example: XDATCAR dump.xyz"
echo " ------------>>"
read -r -a xdatcar_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/xdatcar2exyz.py "${xdatcar_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/xdatcar2exyz.py"
echo " ---------------------------------------------------"
}

# Convert POSCAR to LAMMPS data
function poscar2lammps(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: pos2lmp.py                              |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <POSCAR> <lammps.data>"
echo " Example: POSCAR lammps.data"
echo " ------------>>"
read -r -a pos2lmp_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/pos2lmp.py "${pos2lmp_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2lmp.py"
echo " ---------------------------------------------------"
}

# Convert LAMMPS dump to extxyz
function lammps2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: lmp2exyz.py                             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <dump_file> <element1> <element2> ..."
echo " Example: dump.lammpstrj Li O"
echo " ------------>>"
read -r -a lmp2extxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/lmp2exyz.py "${lmp2extxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/lmp2exyz.py"
echo " ---------------------------------------------------"
}

# Convert ASE traj to extxyz
function traj2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: traj2exyz.py                            |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.traj> <output.xyz>"
echo " Example: input.traj output.xyz"
echo " ------------>>"
read -r -a traj2extxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/traj2exyz.py "${traj2extxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/traj2exyz.py"
echo " ---------------------------------------------------"
}

#--------------------- function 1 ----------------------
function f1_format_conversion(){
    echo " +-------------------------------------------------------------+"
    echo " |                   FORMAT CONVERSION TOOLS                   |"
    echo " +-------------------------------------------------------------+"
    echo " | 101) VASP to extxyz            106) Add group labels        |"
    echo " | 102) MTP to extxyz             107) Add weight to extxyz    |"
    echo " | 103) CP2K to extxyz            108) Extract frame extxyz    |"
    echo " | 104) ABACUS to extxyz          109) Clean XYZ info          |"
    echo " | 105) extxyz to POSCAR          110) Replicate structure     |"
    echo " +-------------------------------------------------------------+"
    echo " | out2exyz) OUTCAR to extxyz     xdat2exyz) XDATCAR to extxyz |"
    echo " | pos2exyz) POSCAR to extxyz     pos2lmp)   POSCAR to LAMMPS  |"
    echo " | cif2pos)  CIF to POSCAR        lmp2exyz)  LAMMPS to extxyz  |"
    echo " | cif2exyz) CIF to extxyz        traj2exyz) ASE traj to extxyz|"
    echo " +-------------------------------------------------------------+"
    echo " | 000) Return to main menu                                    |"
    echo " +-------------------------------------------------------------+"
    echo " Input the function number or converter keyword:"

valid_menu_choices=(
    "000" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" \
    "out2exyz" "pos2exyz" "cif2pos" "cif2exyz" "xdat2exyz" "pos2lmp" "lmp2exyz" "traj2exyz"
) 
read -p " " num_choice
while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number or converter keyword..."
  read -p " " num_choice
done

case $num_choice in
    "101") f101_out2xyz ;;
    "102") f102_mtp2xyz ;;
    "103") f103_cp2k2xyz ;;
    "104") f104_abacus2xyz ;;
    "105") f105_extxyz2poscar ;;
    "106") f106_add_group_labels ;;
    "107") f107_add_weight ;;
    "108") f108_get_frame ;;
    "109") f109_clean_xyz ;;
    "110") f110_replicate_structure ;;
    "out2exyz") out2exyz ;;
    "pos2exyz") pos2extxyz ;;
    "cif2pos") cif2poscar ;;
    "cif2exyz") cif2extxyz ;;
    "xdat2exyz") xdatcar2extxyz ;;
    "pos2lmp") poscar2lammps ;;
    "lmp2exyz") lammps2extxyz ;;
    "traj2exyz") traj2extxyz ;;
    "000") menu; main ;;
esac
}
