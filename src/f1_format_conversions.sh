#--------------------- function 1 format conversion ----------------------
# These functions are used to convert the format of the files

# Convert VASP to extxyz
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

# Convert mtp to extxyz
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

# Convert extxyz to POSCAR
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

# Convert OUTCAR to extxyz (Python version)
function f106_out2exyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: out2exyz.py                             |"
echo " | Developer: Zihan YAN et al.                     |"
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
function f107_pos2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: pos2exyz.py                             |"
echo " | Developer: Zihan YAN                            |"
echo " >-------------------------------------------------<"
echo " Input <POSCAR_pattern> <output.xyz>"
echo " Example: POSCAR model.xyz"
echo " ------------>>"
read -r -a pos2extxyz_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py "${pos2extxyz_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py"
echo " ---------------------------------------------------"
}

# Convert CIF to POSCAR
function f108_cif2poscar(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cif2pos.py                              |"
echo " | Developer: Boyi SITU                            |"
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
function f109_cif2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: cif2exyz.py                             |"
echo " | Developer: Boyi SITU                            |"
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
function f110_xdatcar2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: xdatcar2exyz.py                         |"
echo " | Developer: Zihan YAN                            |"
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
function f111_poscar2lammps(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: pos2lmp.py                              |"
echo " | Developer: Zihan YAN                            |"
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
function f112_lammps2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: lmp2exyz.py                             |"
echo " | Developer: Zihan YAN                            |"
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
function f113_traj2extxyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: traj2exyz.py                            |"
echo " | Developer: Zihan YAN                            |"
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

# Add group labels
function f114_add_group_labels(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: add_groups.py                           |"
echo " | Developer: Zihan YAN                            |"
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
function f115_add_weight(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: add_weight.py                           |"
echo " | Developer: Zihan YAN                            |"
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
function f116_get_frame(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: get_frame.py                            |"
echo " | Developer: Zihan YAN                            |"
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
function f117_clean_xyz(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: clean_xyz.py                            |"
echo " | Developer: Zihan YAN                            |"
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
function f118_replicate_structure(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/format_conversion |"
echo " | Script: replicate.py                            |"
echo " | Developer: Boyi SITU                            |"
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

#--------------------- function 1 ----------------------
function f1_format_conversion(){
echo " ------------>>"
echo " 101) Convert VASP to extxyz"
echo " 102) Convert mtp to extxyz"
echo " 103) Convert CP2K to extxyz"
echo " 104) Convert ABACUS to extxyz"
echo " 105) Convert extxyz to POSCAR"
echo " 106) Add group labels"
echo " 107) Add weight to extxyz"
echo " 108) Extract frame from extxyz"
echo " 109) Clean extra info in XYZ"
echo " 110) Replicate structure"
echo " ---------------- Simple CLI-style converters ----------------"
echo " out2exyz)  Convert OUTCAR to extxyz (also: gpumdkit.sh -out2exyz)"
echo " pos2exyz)  Convert POSCAR to extxyz (also: gpumdkit.sh -pos2exyz)"
echo " cif2pos)   Convert CIF to POSCAR (also: gpumdkit.sh -cif2pos)"
echo " cif2exyz)  Convert CIF to extxyz (also: gpumdkit.sh -cif2exyz)"
echo " xdat2exyz) Convert XDATCAR to extxyz (also: gpumdkit.sh -xdat2exyz)"
echo " pos2lmp)   Convert POSCAR to LAMMPS data (also: gpumdkit.sh -pos2lmp)"
echo " lmp2exyz)  Convert LAMMPS dump to extxyz (also: gpumdkit.sh -lmp2exyz)"
echo " traj2exyz) Convert ASE traj to extxyz (also: gpumdkit.sh -traj2exyz)"
echo " Tip: operations above also support gpumdkit.sh -addgroup/-addweight/-get_frame/-clean_xyz/-replicate"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number or converter keyword:"

array_num_choice=("000" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" "out2exyz" "pos2exyz" "cif2pos" "cif2exyz" "xdat2exyz" "pos2lmp" "lmp2exyz" "traj2exyz") 
read -p " " num_choice
while ! echo "${array_num_choice[@]}" | grep -wq "$num_choice" 
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
    "106") f114_add_group_labels ;;
    "107") f115_add_weight ;;
    "108") f116_get_frame ;;
    "109") f117_clean_xyz ;;
    "110") f118_replicate_structure ;;
    "out2exyz") f106_out2exyz ;;
    "pos2exyz") f107_pos2extxyz ;;
    "cif2pos") f108_cif2poscar ;;
    "cif2exyz") f109_cif2extxyz ;;
    "xdat2exyz") f110_xdatcar2extxyz ;;
    "pos2lmp") f111_poscar2lammps ;;
    "lmp2exyz") f112_lammps2extxyz ;;
    "traj2exyz") f113_traj2extxyz ;;
    "000") menu; main ;;
esac
}
