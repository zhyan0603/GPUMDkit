#!/bin/bash

# You need to set the path of GPUMD and GPUMDkit in your ~/.bashrc, for example
# export GPUMD_path=/d/Westlake/GPUMD
# export GPUMDkit_path=/d/Westlake/Gpumdkit

if [ -z "$GPUMD_path" ] || [ -z "$GPUMDkit_path" ]; then
    echo "Error: GPUMD_path and/or GPUMDkit_path are not set."
    echo "Please set them in your ~/.bashrc, e.g.:"
    echo "  export GPUMD_path=/d/Westlake/GPUMD"
    echo "  export GPUMDkit_path=/d/Westlake/Gpumdkit"
    exit 1
fi

VERSION="1.2.4 (dev) (2025-04-25)"

#--------------------- function 1 format conversion ----------------------
# These functions are used to convert the format of the files

function f101_out2xyz(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: multipleFrames-outcars2nep-exyz.sh      |"
echo " | Developer: Yanzhou WANG (yanzhowang@gmail.com ) |"
echo " >-------------------------------------------------<"
echo " Input the directory containing OUTCARs"
echo " ------------>>"
read -p " " dir_outcars
echo " ---------------------------------------------------"
script_path="${GPUMD_path}/tools/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
alt_path="${GPUMD_path}/tools/Format_Conversion/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
if [[ -f ${script_path} ]]; then
    bash ${script_path} ${dir_outcars}
    echo " Code path: ${script_path}"
    echo " ---------------------------------------------------"
elif [[ -f ${alt_path} ]]; then
    bash ${alt_path} ${dir_outcars}
    echo " Code path: ${alt_path}"
    echo " ---------------------------------------------------"
else
    echo "Error: multipleFrames-outcars2nep-exyz.sh not found"
    return 1
fi
}

function f102_mtp2xyz(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: mtp2xyz.py                              |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>"
echo " Examp: train.cfg Pd Ag"
echo " ------------>>"
read -p " " mtp_variables
echo " ---------------------------------------------------"
script_path="${GPUMD_path}/tools/mtp2xyz/mtp2xyz.py"
alt_path="${GPUMD_path}/tools/Format_Conversion/mtp2xyz/mtp2xyz.py"
if [[ -f ${script_path} ]]; then
    python ${script_path} ${mtp_variables}
    echo " Code path: ${script_path}"
    echo " ---------------------------------------------------"
elif [[ -f ${alt_path} ]]; then
    python ${alt_path} ${mtp_variables}
    echo " Code path: ${alt_path}"
    echo " ---------------------------------------------------"
else
    echo "Error: mtp2xyz.py not found"
    return 1
fi
}

function f103_cp2k2xyz(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: cp2k2xyz.py                             |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input <dir_cp2k> "
echo " Examp: ./cp2k "
echo " ------------>>"
read -p " " dir_cp2k
echo " ---------------------------------------------------"
script_path="${GPUMD_path}/tools/cp2k2xyz/cp2k2xyz.py"
alt_path="${GPUMD_path}/tools/Format_Conversion/cp2k2xyz/cp2k2xyz.py"
if [[ -f ${script_path} ]]; then
    python ${script_path} ${dir_cp2k}
    echo " Code path: ${script_path}"
    echo " ---------------------------------------------------"
elif [[ -f ${alt_path} ]]; then
    python ${alt_path} ${dir_cp2k}
    echo " Code path: ${alt_path}"
    echo " ---------------------------------------------------"
else
    echo "Error: cp2k2xyz.py not found"
    return 1
fi
}

function f104_castep2xyz(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: castep2nep-exyz.sh                      |"
echo " | Developer: Yanzhou WANG (yanzhowang@gmail.com ) |"
echo " >-------------------------------------------------<"
echo " Input <dir_castep>"
echo " Examp: ./castep "
echo " ------------>>"
read -p " " dir_castep
echo " ---------------------------------------------------"
script_path="${GPUMD_path}/tools/castep2exyz/castep2nep-exyz.sh"
alt_path="${GPUMD_path}/tools/Format_Conversion/castep2exyz/castep2nep-exyz.sh"
if [[ -f ${script_path} ]]; then
    bash ${script_path} ${dir_castep}
    echo " Code path: ${script_path}"
    echo " ---------------------------------------------------"
elif [[ -f ${alt_path} ]]; then
    bash ${alt_path} ${dir_castep}
    echo " Code path: ${alt_path}"
    echo " ---------------------------------------------------"
else
    echo "Error: castep2nep-exyz.sh not found"
    return 1
fi
}

function f105_extxyz2poscar(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
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
echo " 101) Convert OUTCAR to extxyz"
echo " 102) Convert mtp to extxyz"
echo " 103) Convert cp2k to extxyz"
echo " 104) Convert castep to extxyz"
echo " 105) Convert extxyz to POSCAR"
echo " 106) Developing ... "
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "101" "102" "103" "104" ) 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "101")
        f101_out2xyz
        ;;
    "102")
        f102_mtp2xyz
        ;;
    "103")
        f103_cp2k2xyz
        ;;
    "104")
        f104_castep2xyz
        ;;
    "105")
        f105_extxyz2poscar
        ;; 
    "106")
        echo " Developing ... "
        ;;
       
    "000")
        menu
        main
        ;;
esac

}

#--------------------- function 2 sample structures ----------------------
# These functions are used to sample structures

function f201_sample_structures(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
echo " | Script: sample_structures.py                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <extxyz_file> <sampling_method> <num_samples> [skip_num]"
echo " [skip_num]: number of initial frames to skip, default value is 0"
echo " Sampling_method: 'uniform' or 'random'"
echo " Examp: train.xyz uniform 50 "
echo " ------------>>"
read -p " " sample_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py ${sample_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py"
echo " ---------------------------------------------------"
}

function f202_pynep_sample_structures(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
echo " | Script: pynep_select_structs.py                 |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <sample.xyz> <train.xyz> <nep_model> <min_dist>"
echo " Examp: dump.xyz train.xyz ./nep.txt 0.01 "
echo " ------------>>"
read -p " " sample_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/pynep_select_structs.py ${sample_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/pynep_select_structs.py"
echo " ---------------------------------------------------"
}

function f203_find_outliers(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: get_max_rmse_xyz.py                     |"
echo " | Developer: Ke XU (kickhsu@gmail.com)            |"
echo " >-------------------------------------------------<"
echo " Input <extxyz_file> <*_train.out> <num_outliers>"
echo " Examp: train.xyz energy_train.out 13 "
echo " ------------>>"
read -p " " maxrmse_choice
echo " ---------------------------------------------------"
script_path="${GPUMD_path}/tools/get_max_rmse_xyz/get_max_rmse_xyz.py"
alt_path="${GPUMD_path}/tools/Analysis_and_Processing/get_max_rmse_xyz/get_max_rmse_xyz.py"
if [[ -f ${script_path} ]]; then
    python ${script_path} ${maxrmse_choice}
    echo " Code path: ${script_path}"
    echo " ---------------------------------------------------"
elif [[ -f ${alt_path} ]]; then
    python ${alt_path} ${maxrmse_choice}
    echo " Code path: ${alt_path}"
    echo " ---------------------------------------------------"
else
    echo "Error: get_max_rmse_xyz.py not found"
    return 1
fi
}

function f204_perturb_structure(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
echo " | Script: perturb_structure.py                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>"
echo " The default paramters for perturb are 20 0.03 0.2 uniform"
echo " Examp: POSCAR 20 0.03 0.2 uniform"
echo " ------------>>"
read -p " " perturb_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py ${perturb_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py"
echo " ---------------------------------------------------"
}

function f205_select_max_force_deviation_structs(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
echo " | Script: select_max_modev.py                     |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <structs_num> <threshold> (eg. 200 0.15)"
echo " ------------>>"
read -p " " modev_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py ${modev_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py"
echo " ---------------------------------------------------"
}

function f2_sample_structures(){
echo " ------------>>"
echo " 201) Sample structures from extxyz"
echo " 202) Sample structures by pynep"
echo " 203) Find the outliers in training set"
echo " 204) Perturb structure"
echo " 205) Select max force deviation structs from active.xyz"
echo " 206) Developing ... "
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "201" "202" "203" "204" "205" "206") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "201")
        f201_sample_structures
        ;;
    "202")
        f202_pynep_sample_structures
        ;;
    "203")
        f203_find_outliers
        ;;
    "204")
        f204_perturb_structure
        ;;
    "205")
        f205_select_max_force_deviation_structs
        ;;
    "206")
        echo " Developing ... "
        ;;
    "000")
        menu
        main
        ;;
esac

}


#--------------------- function 3 workflow ----------------------
# These functions are used to do the workflow
# See the source codes in Scripts/workflow for more details

function f3_workflow_dev(){
echo " ------------>>"
echo " 301) SCF batch pretreatment"
echo " 302) MD sample batch pretreatment (gpumd)"
echo " 303) MD sample batch pretreatment (lmp)"
echo " 304) Developing ... "
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "301" "302" "303" "304") 
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
    "304")
        echo " Developing ... "
        ;;         
    "000")
        menu
        main
        ;;
esac
}

#--------------------- function 4 calculators ----------------------
# These functions are used to do the calculators
# See the source codes in Scripts/calculators for more details

function f401_calc_ionic_conductivity(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_ion_conductivity.py                |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <element> <charge> (eg. Li 1)"
echo " ------------>>"
read -p " " ion_cond_choice
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py ${ion_cond_choice}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py"
echo " ---------------------------------------------------"
}

function f402_calc_properties_with_nep(){
echo " >-------------------------------------------------<"
echo " | This function calls the script in calculators   |"
echo " | Script: calc_properties_with_nep.py             |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
echo " Input <input.xyz> <output.xyz> <nep_model> "
echo " Examp: input.xyz outpt.xyz nep.txt"
echo " ------------>>"
read -p " " input_calc_properties
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py ${input_calc_properties}
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py"
echo " ---------------------------------------------------"
}

function f4_calculators(){
echo " ------------>>"
echo " 401) Calc ionic conductivity"
echo " 402) Calc properties by nep"
echo " 403) Developing ... "
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"

arry_num_choice=("000" "401" "402" "403") 
read -p " " num_choice
while ! echo "${arry_num_choice[@]}" | grep -wq "$num_choice" 
do
  echo " ------------>>"
  echo " Please reinput function number..."
  read -p " " num_choice
done

case $num_choice in
    "401")
        f401_calc_ionic_conductivity
        ;;
    "402")
        f402_calc_properties_with_nep
        ;;
    "403")
        echo " Developing ... "
        ;;            
    "000")
        menu
        main
        ;;
esac
}

#--------------------- clean extra files ----------------------
# This function is used to clean extra files in the current directory
# It will keep some files and delete the rest

function clean_extra_files(){
keep_files=("run.in" "nep.in" "model.xyz" "nep.txt" "train.xyz" "test.xyz")
keep_patterns=("*sub*" "*.sh" "*slurm")
all_files=$(ls)
delete_files=()

for file in $all_files; do
    # Check if the file matches any keep_files or keep_patterns
    keep=false
    for keep_file in "${keep_files[@]}"; do
        if [[ "$file" == "$keep_file" ]]; then
            keep=true
            break
        fi
    done
    for keep_pattern in "${keep_patterns[@]}"; do
        if [[ "$file" == $keep_pattern ]]; then
            keep=true
            break
        fi
    done

    # If the file is not marked to keep, add it to delete_files
    if [ "$keep" == false ]; then
        delete_files+=("$file")
    fi
done

# Display files to delete
if [ ${#delete_files[@]} -eq 0 ]; then
    echo "No files to delete."
    exit 0
else
    echo "The following files will be deleted:"
    echo "---------------------------------------"
    for file in "${delete_files[@]}"; do
        echo -n "$file "
    done
    echo -e "\n---------------------------------------"
fi

# Ask user for confirmation or additional files to keep
echo "Do you want to delete all these files?"
echo "y/yes to delete, n/no to cancel, or input files to keep (separated by spaces):"
read user_input

# Process user input
if [[ "$user_input" == "y" || "$user_input" == "yes" ]]; then
    echo "Deleting the files..."
    for file in "${delete_files[@]}"; do
        rm -f "$file"
    done
    echo "Files deleted."
elif [[ "$user_input" == "n" || "$user_input" == "no" ]]; then
    echo "Operation canceled. No files were deleted."
    exit 0
else
    # Add extra files to keep based on user input
    extra_keep_files=($user_input)
    for extra_file in "${extra_keep_files[@]}"; do
        delete_files=("${delete_files[@]/$extra_file}")
    done

    # Delete remaining files
    if [ ${#delete_files[@]} -eq 0 ]; then
        echo "No files to delete after processing extra keep files."
    else
        echo "Deleting remaining files..."
        for file in "${delete_files[@]}"; do
            if [ -n "$file" ]; then
                rm -f "$file"
            fi
        done
        echo "Files deleted."
    fi
fi

}


#--------------------- main script ----------------------
# Show the menu
function menu(){
echo " ----------------------- GPUMD -----------------------"
echo " 1) Format Conversion          2) Sample Structures   "
echo " 3) Workflow (dev)             4) Calculators         "
echo " 5) Developing ...             6) Developing ...      "
echo " 0) Quit!"
}

# Function main
function main(){
    echo " ------------>>"
    echo ' Input the function number:'
    array_choice=(
        "0" "1" "101" "102" "103" "104" "105" 
        "2" "201" "202" "203" "204" "205" 
        "3" "301" "302" "303" 
        "4" "401" "402" 
        "5"
        "6"
    ) 
    read -p " " choice
    while ! echo "${array_choice[@]}" | grep -wq "$choice" 
    do
      echo " ------------>>"
      echo " Please reinput function number:"
      read -p " " choice
    done

    case $choice in
        "0")
            echo " Thank you for using GPUMDkit. Have a great day!"
            exit 0
            ;;
        "1")
            f1_format_conversion
            ;;
        "101")
            f101_out2xyz
            ;;
        "102")
            f102_mtp2xyz
            ;;
        "103")
            f103_cp2k2xyz
            ;;
        "104")      
            f104_castep2xyz
            ;;
        "105")  
            f105_extxyz2poscar
            ;;
        "2")
            f2_sample_structures
            ;;
        "201")
            f201_sample_structures
            ;; 
        "202")  
            f202_pynep_sample_structures
            ;;
        "203")
            f203_find_outliers
            ;;
        "204")
            f204_perturb_structure
            ;;
        "205")
            f205_select_max_force_deviation_structs
            ;;
        "3")
            f3_workflow_dev
            ;;
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
        "4")
            f4_calculators
            ;;
        "401")
            f401_calc_ionic_conductivity
            ;;
        "402")
            f402_calc_properties_with_nep
            ;;
        "5")
            echo "Developing ..."
            ;;
        "6")
            echo "Developing ..."
            ;;
        *)
            echo "Incorrect Options"
            ;;

    esac
    echo " Thank you for using GPUMDkit. Have a great day!"
}

#--------------------- help info ----------------------
# This function is used to show the help information
# It will show the usage of each function

function help_info_table(){
    echo "+==================================================================================================+"
    echo "|                              GPUMDkit ${VERSION} Usage                             |"
    echo "|                                                                 --- by Zihan YAN                 |"
    echo "+======================================== Conversions =============================================+"
    echo "| -outcar2exyz   Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           |"
    echo "| -castep2exyz   Convert castep to extxyz       | -pos2lmp      Convert POSCAR to LAMMPS           |"
    echo "| -cp2k2exyz     Convert cp2k output to extxyz  | -lmp2exyz     Convert LAMMPS-dump to extxyz      |"
    echo "| -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz |"
    echo "| Developing...                                 | Developing...                                    |"
    echo "+========================================= Analysis ===============================================+"
    echo "| -range         Print range of energy etc.     | -max_rmse     Get max RMSE from XYZ              |"
    echo "| -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       |"
    echo "| -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       |"
    echo "| -filter_dist   Filter struct by min_dist      | Developing...                                    |"
    echo "+=========================================    Misc  ==============+================================+"
    echo "| -plt           Plot scripts                   | -get_frame     Extract the specified frame       |"
    echo "| -calc          Calculators                    | -clear_xyz     Clear extra info in XYZ file      |"
    echo "| -clean         Clear files for work_dir       | -time          Time consuming Analyzer           |"
    echo "| Developing...                                 | Developing...                                    |"    
    echo "+==================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   |"
    echo "+==================================================================================================+"
}

if [ ! -z "$1" ]; then
    case $1 in
        -h|-help)
            help_info_table
            ;;
        -clean)
            clean_extra_files
            ;;
        -time)
            case $2 in
                gpumd)
                    bash ${GPUMDkit_path}/Scripts/analyzer/time_consuming_gpumd.sh $3
                    ;;
                nep)
                    bash ${GPUMDkit_path}/Scripts/analyzer/time_consuming_nep.sh
                    ;;                
                *)
                    echo " See the codes in analyzer folder for more details"
                    echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/time_consuming_*.sh"
                    exit 1
                    ;;
            esac
            ;;
        -plt)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    "thermo")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_nep_thermo.py $3
                        ;;
                    "train")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_nep_train_results.py $3
                        ;;  
                    "prediction"|"valid"|"test")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_nep_prediction_results.py $3
                        ;;
                    "train_test"|"tt")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_nep_train_test.py $3
                        ;;
                    "msd")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_msd.py $3
                        ;;
                    "sdc")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_sdc.py $3
                        ;;
                    "rdf")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_rdf.py $3 $4
                        ;;
                    "vac")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_vac.py $3
                        ;;                
                    "restart")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_nep_restart.py $3
                        ;;
                    "dimer")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_dimer.py $3 $4 $5 $6
                        ;;
                    "force_errors"|"force_error"|"force")
                        python ${GPUMDkit_path}/Scripts/plt_scripts/plt_force_errors.py $3 
                        ;;    
                    *)
                        echo "Usage: -plt thermo/train/prediction/train_test/msd/sdc/rdf/vac/restart/dimer/force [save]"
                        echo "Examp: gpumdkit.sh -plt thermo save"
                        exit 1
                        ;;
                esac
            else
                echo " Usage: -plt thermo/train/prediction/train_test/msd/vac/sdc/rdf/vac/restart/dimer/force [save] (eg. gpumdkit.sh -plt thermo)"
                echo " See the codes in plt_scripts for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/plt_scripts"
            fi
            ;;

        -calc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    ion-cond|ionic-cond|ionic-conductivity)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] ; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py"
                            python ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py $3 $4
                        else
                            echo " Usage: -calc ion-cond <element> <charge>"
                            echo " Examp: gpumdkit.sh -calc ion-cond Li 1"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_ion_conductivity.py"
                            exit 1
                        fi
                        ;;
                    nep)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ]; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py"
                            python ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py $3 $4 $5
                        else
                            echo " Usage: -calc nep <input.xyz> <output.xyz> <nep_model>"
                            echo " Examp: gpumdkit.sh -calc nep input.xyz output.xyz nep.txt"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_properties_with_nep.py"
                            exit 1
                        fi
                        ;;
                    des)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ] && [ ! -z "$6" ]; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_descriptors.py"
                            python ${GPUMDkit_path}/Scripts/calculators/calc_descriptors.py $3 $4 $5 $6
                        else
                            echo " Usage: -calc des <input.xyz> <output.npy> <nep_model> <element>"
                            echo " Examp: gpumdkit.sh -calc des train.xyz des_Li.npy nep.txt Li"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_descriptors.py"
                            exit 1
                        fi
                        ;;                  
                    *)
                        echo " See the codes in calculators folder for more details"
                        echo " Code path: ${GPUMDkit_path}/Scripts/calculators"
                        exit 1
                        ;;
                esac
            else
                echo " See the codes in calculators folder for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/calculators"
            fi
            ;;

        -range)
            if [ ! -z "$2" ] && [ ! -z "$3" ] && [ "$2" != "-h" ]  ; then
                echo " Calling script by Zihan YAN. "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py"
                python ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py $2 $3 ${@:4}
            else
                echo " Usage: -range <exyzfile> <property> [hist] (eg. gpumdkit.sh -range train.xyz energy hist)" 
                echo " See the source code of energy_force_virial_analyzer.py for more details"
                echo " Code path: Code path: ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py"
            fi
            ;;

        -out2xyz|-outcar2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Yanzhou WANG et al. "
                script_path="${GPUMD_path}/tools/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
                alt_path="${GPUMD_path}/tools/Format_Conversion/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
                if [[ -f ${script_path} ]]; then
                    echo " Code path: ${script_path}"
                    bash ${script_path} $2
                elif [[ -f ${alt_path} ]]; then
                    echo " Code path: ${alt_path}"
                    bash ${alt_path} $2
                else
                    echo "Error: multipleFrames-outcars2nep-exyz.sh not found"
                    return 1
                fi
            else
                echo " Usage: -out2xyz|-outcar2exyz dir_name (eg. gpumdkit.sh -outcar2exyz .)"
                echo " See the source code of multipleFrames-outcars2nep-exyz.sh for more details"
                script_path="${GPUMD_path}/tools/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
                alt_path="${GPUMD_path}/tools/Format_Conversion/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
                if [[ -f "${script_path}" ]]; then
                    echo " Code path: ${script_path}"
                elif [[ -f "${alt_path}" ]]; then
                    echo " Code path: ${alt_path}"
                else
                    echo "Error: multipleFrames-outcars2nep-exyz.sh not found"
                    return 1
                fi
            fi
            ;;

        -cast2xyz|-castep2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Yanzhou WANG et al. "
                script_path="${GPUMD_path}/tools/castep2exyz/castep2nep-exyz.sh"
                alt_path="${GPUMD_path}/tools/Format_Conversion/castep2exyz/castep2nep-exyz.sh"
                if [[ -f ${script_path} ]]; then
                    echo " Code path: ${script_path}"
                    bash ${script_path} $2
                elif [[ -f ${alt_path} ]]; then
                    echo " Code path: ${alt_path}"
                    bash ${alt_path} $2
                else
                    echo "Error: castep2nep-exyz.sh not found"
                    return 1
                fi
            else
                echo " Usage: -cast2xyz|-castep2exyz dir_name (eg. gpumdkit.sh -castep2exyz .)"
                echo " See the source code of castep2nep-exyz.sh for more details"
                script_path="${GPUMD_path}/tools/castep2exyz/castep2nep-exyz.sh"
                alt_path="${GPUMD_path}/tools/Format_Conversion/castep2exyz/castep2nep-exyz.sh"
                if [[ -f "${script_path}" ]]; then
                    echo " Code path: ${script_path}"
                elif [[ -f "${alt_path}" ]]; then
                    echo " Code path: ${alt_path}"
                else
                    echo "Error: castep2nep-exyz.sh not found"
                    return 1
                fi
            fi
            ;;

        -cp2k2xyz|-cp2k2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Ke XU et al. "
                script_path="${GPUMD_path}/tools/cp2k2xyz/cp2k2xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/cp2k2xyz/cp2k2xyz.py"
                if [[ -f ${script_path} ]]; then
                    echo " Code path: ${script_path}"
                    python ${script_path} $2
                elif [[ -f ${alt_path} ]]; then
                    echo " Code path: ${alt_path}"
                    python ${alt_path} $2
                else
                    echo "Error: cp2k2xyz.py not found"
                    return 1
                fi
            else
                echo " Usage: -cp2k2xyz|-cp2k2exyz dir_name (eg. gpumdkit.sh -cp2k2exyz .)"
                echo " See the source code of cp2k2xyz.py for more details"
                script_path="${GPUMD_path}/tools/cp2k2xyz/cp2k2xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/cp2k2xyz/cp2k2xyz.py"                
                if [[ -f "${script_path}" ]]; then
                    echo " Code path: ${script_path}"
                elif [[ -f "${alt_path}" ]]; then
                    echo " Code path: ${alt_path}"
                else
                    echo "Error: cp2k2xyz.py not found"
                    return 1
                fi
            fi
            ;;

        -mtp2xyz|-mtp2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Ke XU et al. "
                script_path="${GPUMD_path}/tools/mtp2xyz/mtp2xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/mtp2xyz/mtp2xyz.py"
                if [[ -f ${script_path} ]]; then
                    echo " Code path: ${script_path}"
                    python ${script_path} train.cfg $2 ${@:3}
                elif [[ -f ${alt_path} ]]; then
                    echo " Code path: ${alt_path}"
                    python ${alt_path} train.cfg $2 ${@:3}
                else
                    echo "Error: mtp2xyz.py not found"
                    return 1
                fi
            else
                echo " Usage: -mtp2xyz|-mtp2exyz train.cfg Symbol1 Symbol2 Symbol3 ..."
                echo "   Examp: gpumdkit.sh -mtp2exyz train.cfg Pd Ag"
                echo " See the source code of mtp2xyz.py for more details"
                script_path="${GPUMD_path}/tools/mtp2xyz/mtp2xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/mtp2xyz/mtp2xyz.py"
                if [[ -f "${script_path}" ]]; then
                    echo " Code path: ${script_path}"
                elif [[ -f "${alt_path}" ]]; then
                    echo " Code path: ${alt_path}"
                else
                    echo "Error: mtp2xyz.py not found"
                    return 1
                fi
            fi
            ;;

        -pos2exyz)
            if [ ! -z "$2" ] && [ ! -z "$3" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py $2 $3
            else
                echo " Usage: -pos2exyz POSCAR model.xyz"
                echo " See the source code of pos2exyz.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2exyz.py"
            fi
            ;;

        -exyz2pos)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py $2
            else
                echo " Usage: -exyz2pos model.xyz"
                echo " See the source code of exyz2pos.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/exyz2pos.py"
            fi
            ;;

        -pos2lmp)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2lmp.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/pos2lmp.py $2 $3
            else
                echo " Usage: -pos2lmp POSCAR lammps.data"
                echo " See the source code of pos2lmp.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/pos2lmp.py"
            fi
            ;;

        -lmp2exyz|-lmpdump2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/lmp2exyz.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/lmp2exyz.py $2 ${@:3}
            else
                echo " Usage: -lmp2exyz <dump_file> <element1> <element2> ..."
                echo " See the source code of lmp2exyz.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/lmp2exyz.py"
            fi
            ;;

        -addgroup|-addlabel)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_groups.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/add_groups.py $2 ${@:3}
            else
                echo " Usage: -addgroup <POSCAR> <element1> <element2> ..."
                echo " See the source code of add_groups.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_groups.py"
            fi
            ;;

        -addweight)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] && [ ! -z "$4" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_weight.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/add_weight.py $2 $3 $4
            else
                echo " Usage: -addweight <input.xyz> <output.xyz> <weight> "
                echo " See the source code of add_groups.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/add_weight.py"
            fi
            ;;

        -max_rmse|-get_max_rmse_xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] && [ ! -z "$4" ]; then
                echo " Calling script by Ke XU "
                script_path="${GPUMD_path}/tools/get_max_rmse_xyz/get_max_rmse_xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/get_max_rmse_xyz/get_max_rmse_xyz.py"
                if [[ -f ${script_path} ]]; then
                    echo " Code path: ${script_path}"
                    python ${script_path} $2 $3 $4
                elif [[ -f ${alt_path} ]]; then
                    echo " Code path: ${alt_path}"
                    python ${alt_path} $2 $3 $4
                else
                    echo "Error: get_max_rmse_xyz.py not found"
                    return 1
                fi                
            else
                echo " Usage: -getmax|-get_max_rmse_xyz train.xyz force_train.out 13"
                echo " See the source code of get_max_rmse_xyz.py for more details"
                script_path="${GPUMD_path}/tools/get_max_rmse_xyz/get_max_rmse_xyz.py"
                alt_path="${GPUMD_path}/tools/Format_Conversion/get_max_rmse_xyz/get_max_rmse_xyz.py"                
                if [[ -f "${script_path}" ]]; then
                    echo " Code path: ${script_path}"
                elif [[ -f "${alt_path}" ]]; then
                    echo " Code path: ${alt_path}"
                else
                    echo "Error: get_max_rmse_xyz.py not found"
                    return 1
                fi
            fi
            ;;

        -min_dist)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist.py"
                python ${GPUMDkit_path}/Scripts/analyzer/get_min_dist.py $2
            else
                echo " Usage: -min_dist <exyzfile>"
                echo " See the source code of get_min_dist.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist.py"
            fi
            ;;

        -min_dist_pbc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist_pbc.py"
                python ${GPUMDkit_path}/Scripts/analyzer/get_min_dist_pbc.py $2
            else
                echo " Usage: -min_dist_pbc <exyzfile>"
                echo " See the source code of get_min_dist_pbc.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/get_min_dist_pbc.py"
            fi
            ;;

        -filter_dist)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance.py"
                python ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance.py $2 $3
            else
                echo " Usage: -filter_xyz <exyzfile> <min_dist>"
                echo " See the source code of filter_structures_by_distance.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance.py"
            fi
            ;;

        -filter_dist_pbc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance_pbc.py"
                python ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance_pbc.py $2 $3
            else
                echo " Usage: -filter_xyz_pbc <exyzfile> <min_dist>"
                echo " See the source code of filter_structures_by_distance_pbc.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_structures_by_distance_pbc.py"
            fi
            ;;

        -filter_box)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_box.py"
                python ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_box.py $2 $3
            else
                echo " Usage: -filter_box <exyzfile> <lattice limit>"
                echo " See the source code of filter_exyz_by_box.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_box.py"
            fi
            ;;

        -filter_value)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_value.py"
                python ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_value.py $2 $3 $4
            else
                echo " Usage: -filter_value <exyzfile> <property> <value>"
                echo " See the source code of filter_exyz_by_value.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/analyzer/filter_exyz_by_value.py"
            fi
            ;;

        -get_frame)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/get_frame.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/get_frame.py $2 $3
            else
                echo " Usage: -get_frame <exyzfile> <frame_index>"
                echo " See the source code of get_frame.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/get_frame.py"
            fi
            ;;
        -clear_xyz|-clean_xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/clear_xyz.py"
                python ${GPUMDkit_path}/Scripts/format_conversion/clear_xyz.py $2 $3
            else
                echo " Usage: -clear_xyz <input.xyz> <output.xyz>"
                echo " See the source code of clear_xyz.py for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/clear_xyz.py"
            fi
            ;;

        *)
            echo " Unknown option: $1 "
            help_info_table
            exit 1
            ;;
    esac
    exit
fi

## logo
echo -e "\
         ____ ____  _   _ __  __ ____  _    _ _   
        / ___|  _ \| | | |  \/  |  _ \| | _(_) |_ 
       | |  _| |_) | | | | |\/| | | | | |/ / | __|
       | |_| |  __/| |_| | |  | | |_| |   <| | |_ 
        \____|_|    \___/|_|  |_|____/|_|\_\_|\__|
                                          
        GPUMDkit Version ${VERSION}
     Developer: Zihan YAN (yanzihan@westlake.edu.cn)
      "
menu
main