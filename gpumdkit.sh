#!/bin/bash

# You need to set the path of GPUMD and GPUMDkit in your ~/.bashrc, for example
# export GPUMDkit_path=/home/yanzihan/software/GPUMDkit

if [ -z "$GPUMDkit_path" ]; then
    echo "Error: GPUMDkit_path is not set."
    echo "Please set it in your ~/.bashrc, e.g.:"
    echo "  export GPUMDkit_path=/home/yanzihan/software/GPUMDkit"
    exit 1
fi

VERSION="1.4.3 (dev) (2025-12-28)"

plt_path="${GPUMDkit_path}/Scripts/plt_scripts"
analyzer_path="${GPUMDkit_path}/Scripts/analyzer"
calc_path="${GPUMDkit_path}/Scripts/calculators"
workflow_path="${GPUMDkit_path}/Scripts/workflow"
format_conv_path="${GPUMDkit_path}/Scripts/format_conversion"
utils_path="${GPUMDkit_path}/Scripts/utils"
sample_path="${GPUMDkit_path}/Scripts/sample_structures"

#--------------------- main script ----------------------
# Show the menu
function menu(){
echo " ----------------------- GPUMD -----------------------"
echo " 1) Format Conversion          2) Sample Structures   "
echo " 3) Workflow                   4) Calculators         "
echo " 5) Analyzer                   6) Developing ...      "
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
        "4" "401" "402" "403" "404" "405"
        "5" "501" "502"
        "6"
    ) 
    read -p " " choice
    while ! echo "${array_choice[@]}" | grep -wq "$choice" 
    do
      echo " ------------>>"
      echo " Please reinput function number:"
      read -p " " choice
    done

    case "${choice:0:1}" in
        "0")
            echo " Thank you for using GPUMDkit. Have a great day!"
            exit 0
            ;;
        "1")
            source ${GPUMDkit_path}/src/f1_format_conversions.sh
            case $choice in
                "1") f1_format_conversion ;;
                "101") f101_out2xyz ;;
                "102") f102_mtp2xyz ;;
                "103") f103_cp2k2xyz ;;
                "104") f104_abacus2xyz ;;
                "105") f105_extxyz2poscar ;;
            esac ;;
        "2")
            source ${GPUMDkit_path}/src/f2_sample_structures.sh
            case $choice in
                "2") f2_sample_structures ;;
                "201") f201_sample_structures ;; 
                "202") f202_pynep_sample_structures ;;
                "203") f203_neptrain_sample_structures ;;
                "204") f204_perturb_structure ;;
                "205") f205_select_max_force_deviation_structs ;;
            esac ;;
        "3")
            source ${GPUMDkit_path}/src/f3_workflows.sh
            case $choice in
                "3") f3_workflow_dev ;;
                "301") 
                    f301_scf_batch_pretreatment ;;
                "302") 
                    source ${workflow_path}/md_sample_batch_pretreatment_gpumd.sh
                    f302_md_sample_batch_pretreatment_gpumd ;;
                "303") 
                    source ${workflow_path}/md_sample_batch_pretreatment_lmp.sh
                    f303_md_sample_batch_pretreatment_lmp ;;
            esac ;;
        "4")
            source ${GPUMDkit_path}/src/f4_calculators.sh
            case $choice in
                "4") f4_calculators ;;
                "401") f401_calc_ionic_conductivity ;;
                "402") f402_calc_properties_with_nep ;;
                "403") f403_calc_descriptors ;;
                "404") f404_calc_doas ;;
                "405") f405_calc_neb ;;                
            esac ;;           
        "5")
            source ${GPUMDkit_path}/src/f5_analyzers.sh
            case $choice in
                "5") f5_analyzers ;;
                "501") f501_analyze_composition ;;
                "502") f502_find_outliers ;;
            esac ;;  
        "6")
            echo " Developing ..." ;;
        *)
            echo " Incorrect Options"
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
    echo "+======================================== Conversions =============================================+"
    echo "| -out2xyz       Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           |"
    echo "| -cif2pos       Convert cif to POSCAR          | -pos2lmp      Convert POSCAR to LAMMPS           |"
    echo "| -cif2exyz      Convert cif to extxyz          | -lmp2exyz     Convert LAMMPS-dump to extxyz      |"
    echo "| -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz |"
    echo "| Developing...                                 | Developing...                                    |"
    echo "+========================================= Analysis ===============================================+"
    echo "| -range         Print range of energy etc.     | -max_rmse     Get max RMSE from extxyz           |"
    echo "| -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       |"
    echo "| -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       |"
    echo "| -filter_dist   Filter struct by min_dist      | -analyze_comp Analyze composition of extxyz      |"
    echo "| -pynep         Sample struct by pynep         | Developing...                                    |"
    echo "+====================================== Misc Utilities ============================================+"
    echo "| -plt           Plot scripts                   | -get_frame     Extract the specified frame       |"
    echo "| -calc          Calculators                    | -frame_range   Extract frames by fraction range  |"
    echo "| -clean         Clear files for work_dir       | -clean_xyz     Clean extra info in XYZ file      |"
    echo "| -time          Time consuming Analyzer        | -update        Update GPUMDkit                   |"    
    echo "+==================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   |"
    echo "+==================================================================================================+"
}

function plot_info_table(){
    echo "+=====================================================================================================+"
    echo "|                              GPUMDkit ${VERSION} Plotting Usage                       |"
    echo "+=============================================== Plot Types ==========================================+"
    echo "| thermo          Plot thermo info                   | train          Plot NEP train results          |"
    echo "| prediction      Plot NEP prediction results        | train_test     Plot NEP train and test results |"
    echo "| msd             Plot mean square displacement      | msd_conv       Plot the convergence of MSD     |"
    echo "| msd_all         Plot MSD of all species            | sdc            Plot self diffusion coefficient |"
    echo "| rdf             Plot radial distribution function  | vac            Plot velocity autocorrelation   |"
    echo "| restart         Plot parameters in nep.restart     | dimer          Plot dimer plot                 |"
    echo "| force_errors    Plot force errors                  | des            Plot descriptors                |"
    echo "| charge          Plot charge distribution           | lr             Plot learning rate              |"
    echo "| doas            Plot density of atomistic states   | arrhenius_d    Plot Arrhenius diffusivity      |"
    echo "| arrhenius_sigma Plot Arrhenius sigma               | net_force      Plot net force distribution     |"
    echo "| emd             Plot EMD results                   | nemd           Plot NEMD results               |"
    echo "| hnemd           Plot HNEMD results                 |                                                |"
    echo "+=====================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -plt <plot_type> -h                               |"
    echo "+=====================================================================================================+"
}

#--------------------- command line ----------------------
if [ ! -z "$1" ]; then
    case $1 in
        -h|-help) help_info_table ;;
        -clean) 
            source ${utils_path}/clean_extra_files.sh
            clean_extra_files ;;
        -update|-U)
            source ${utils_path}/update_gpumdkit.sh
            update_gpumdkit
            source $GPUMDkit_path/docs/updates.info ;;            
        -time)
            case $2 in
                gpumd) bash ${analyzer_path}/time_consuming_gpumd.sh ;;
                nep|gnep) bash ${analyzer_path}/time_consuming_nep.sh ;;                
                *)
                    echo " See the codes in analyzer folder for more details"
                    echo " Code path: ${analyzer_path}/time_consuming_*.sh"
                    exit 1 ;;
            esac ;;
        -plt)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    "thermo") python ${plt_path}/plt_thermo.py $3 ;;
                    "thermo2") python ${plt_path}/plt_thermo2.py $3 ;;
                    "thermo3") python ${plt_path}/plt_thermo3.py $3 ;;                        
                    "train") python ${plt_path}/plt_train.py $3 ;;                 
                    "prediction"|"test") python ${plt_path}/plt_prediction.py $3 ;; 
                    "parity_density") python ${plt_path}/plt_parity_density.py $3 ;;
                    "train_test") python ${plt_path}/plt_train_test.py $3 ;;
                    "msd") python ${plt_path}/plt_msd.py $3 ;;
                    "msd_all") python ${plt_path}/plt_msd_all.py $3 ${@:4} ;;
                    "msd_conv") python ${plt_path}/plt_msd_convergence_check.py $3 ;;
                    "sdc") python ${plt_path}/plt_sdc.py $3 ;;
                    "rdf") python ${plt_path}/plt_rdf.py $3 $4 ;;
                    "vac") python ${plt_path}/plt_vac.py $3 ;;
                    "restart") python ${plt_path}/plt_nep_restart.py $3 ;;
                    "dimer") python ${plt_path}/plt_dimer.py $3 $4 $5 $6 ;;
                    "force_errors") python ${plt_path}/plt_force_errors.py $3 ;;
                    "des") python ${plt_path}/plt_descriptors.py $3 ${@:4} ;;
                    "lr") python ${plt_path}/plt_learning_rate.py $3 ;;
                    "doas") python ${plt_path}/plt_doas.py $3 $4 ;;
                    "arrhenius_d"|"D") python ${plt_path}/plt_arrhenius_d.py $3 ;;
                    "arrhenius_sigma"|"sigma") python ${plt_path}/plt_arrhenius_sigma.py $3 ;;
                    "net_force") python ${plt_path}/plt_net_force.py ${@:3} ;;
                    "emd") python ${plt_path}/plt_emd.py ${@:3} ;;
                    "nemd") python ${plt_path}/plt_nemd.py ${@:3} ;;
                    "hnemd") python ${plt_path}/plt_hnemd.py ${@:3} ;;
                    "charge")
                        echo " +----------------------------------------------------------+"
                        echo " | Please ensure you are using full batch training process. |"
                        echo " | If not, run the prediction step before plotting to avoid |"
                        echo " | inconsistencies in the atomic order between the training |"
                        echo " | set and charge_train.out.                                |"
                        echo " +----------------------------------------------------------+"
                        python ${plt_path}/plt_charge.py $3 ;;
                    *) plot_info_table; exit 1 ;;
                esac
            else
                plot_info_table
                echo " See the codes in plt_scripts for more details"
                echo " Code path: ${GPUMDkit_path}/Scripts/plt_scripts"
            fi ;;

        -calc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                case $2 in
                    ionic-cond)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] ; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${calc_path}/calc_ion_conductivity.py"
                            python ${calc_path}/calc_ion_conductivity.py $3 $4
                        else
                            echo " Usage: -calc ion-cond <element> <charge>"
                            echo " Examp: gpumdkit.sh -calc ion-cond Li 1"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${calc_path}/calc_ion_conductivity.py"
                            exit 1
                        fi ;;
                    nep)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ]; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${calc_path}/calc_properties_with_nep.py"
                            python ${calc_path}/calc_properties_with_nep.py $3 $4 $5
                        else
                            echo " Usage: -calc nep <input.xyz> <output.xyz> <nep_model>"
                            echo " Examp: gpumdkit.sh -calc nep input.xyz output.xyz nep.txt"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${calc_path}/calc_properties_with_nep.py"
                            exit 1
                        fi ;;
                    des)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ] && [ ! -z "$6" ]; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${calc_path}/calc_descriptors.py"
                            python ${calc_path}/calc_descriptors.py $3 $4 $5 $6
                        else
                            echo " Usage: -calc des <input.xyz> <output.npy> <nep_model> <element>"
                            echo " Examp: gpumdkit.sh -calc des train.xyz des_Li.npy nep.txt Li"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${calc_path}/calc_descriptors.py"
                            exit 1
                        fi ;; 
                    doas)
                        if [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ]; then
                            echo " Calling script by Zihan YAN. "
                            echo " Code path: ${calc_path}/calc_doas.py"
                            python ${calc_path}/calc_doas.py $3 $4 $5
                        else
                            echo " Usage: -calc doas <input.xyz> <nep_model> <output_file>"
                            echo " Examp: gpumdkit.sh -calc doas dump.xyz nep.txt doas.out"
                            echo " See the codes in calculators folder for more details"
                            echo " Code path: ${calc_path}/calc_doas.py"
                            exit 1
                        fi ;;                                          
                    *)
                        echo " See the codes in calculators folder for more details"
                        echo " Code path: ${calc_path}"; exit 1 ;;
                esac
            else
                echo " See the codes in calculators folder for more details"
                echo " Code path: ${calc_path}"
            fi ;;

        -range)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]  ; then
                echo " Calling script by Zihan YAN. "
                echo " Code path: ${analyzer_path}/energy_force_virial_analyzer.py"
                python ${analyzer_path}/energy_force_virial_analyzer.py $2 $3 ${@:4}
            else
                echo " Usage: -range <exyzfile> <property> [hist] (eg. gpumdkit.sh -range train.xyz energy hist)" 
                echo " See the source code of energy_force_virial_analyzer.py for more details"
                echo " Code path: Code path: ${analyzer_path}/energy_force_virial_analyzer.py"
            fi ;;

        -replicate)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Boyi SITU. "
                echo " Code path: ${format_conv_path}/replicate.py"
                python ${format_conv_path}/replicate.py $2 $3 ${@:4}
            else
                echo " Please give the file name suffix (e.g. input.xyz or output.vasp)"
                echo " Usage 1: -replicate <inputfile> <outputfile> a b c" 
                echo " Usage 2: -replicate <inputfile> <outputfile> target_num"
                echo " See the source code of replicate.py for more details"
                echo " Code path: Code path: ${format_conv_path}/replicate.py"
            fi ;;

        -out2xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Yanzhou WANG et al. "
                bash ${format_conv_path}/out2xyz.sh $2
                echo " Code path: ${format_conv_path}/out2xyz.sh"
            else
                echo " Usage: -out2xyz dir_name (eg. gpumdkit.sh -out2xyz .)"
                echo " See the source code of out2xyz.sh for more details"
                echo " Code path: ${format_conv_path}/out2xyz.sh"
            fi ;;

        -xml2xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zezhu Zeng et al. "
                python ${format_conv_path}/xml2xyz.py $2
                echo " Code path: ${format_conv_path}/xml2xyz.py"
            else
                echo " Usage: -xml2xyz dir_name (eg. gpumdkit.sh -xml2xyz .)"
                echo " See the source code of xml2xyz.py for more details"
                echo " Code path: ${format_conv_path}/xml2xyz.py"               
            fi ;;

        -pos2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/pos2exyz.py"
                python ${format_conv_path}/pos2exyz.py $2 $3
            else
                echo " Usage: -pos2exyz POSCAR model.xyz"
                echo " See the source code of pos2exyz.py for more details"
                echo " Code path: ${format_conv_path}/pos2exyz.py"
            fi ;;

        -cif2pos)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Boyi SITU "
                echo " Code path: ${format_conv_path}/cif2pos.py"
                python ${format_conv_path}/cif2pos.py $2 $3
            else
                echo " Usage: -cif2pos input.cif POSCAR.vasp"
                echo " See the source code of cif2pos.py for more details"
                echo " Code path: ${format_conv_path}/cif2pos.py"
            fi ;;

        -cif2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Boyi SITU "
                echo " Code path: ${format_conv_path}/cif2exyz.py"
                python ${format_conv_path}/cif2exyz.py $2 $3
            else
                echo " Usage: -cif2exyz input.cif model.xyz"
                echo " See the source code of cif2exyz.py for more details"
                echo " Code path: ${format_conv_path}/cif2exyz.py"
            fi ;;

        -exyz2pos)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/exyz2pos.py"
                python ${format_conv_path}/exyz2pos.py $2
            else
                echo " Usage: -exyz2pos model.xyz"
                echo " See the source code of exyz2pos.py for more details"
                echo " Code path: ${format_conv_path}/exyz2pos.py"
            fi ;;

        -pos2lmp)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/pos2lmp.py"
                python ${format_conv_path}/pos2lmp.py $2 $3
            else
                echo " Usage: -pos2lmp POSCAR lammps.data"
                echo " See the source code of pos2lmp.py for more details"
                echo " Code path: ${format_conv_path}/pos2lmp.py"
            fi ;;

        -lmp2exyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/lmp2exyz.py"
                python ${format_conv_path}/lmp2exyz.py $2 ${@:3}
            else
                echo " Usage: -lmp2exyz <dump_file> <element1> <element2> ..."
                echo " See the source code of lmp2exyz.py for more details"
                echo " Code path: ${format_conv_path}/lmp2exyz.py"
            fi ;;

        -addgroup|-addlabel)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/add_groups.py"
                python ${format_conv_path}/add_groups.py $2 ${@:3}
            else
                echo " Usage: -addgroup <POSCAR> <element1> <element2> ..."
                echo " See the source code of add_groups.py for more details"
                echo " Code path: ${format_conv_path}/add_groups.py"
            fi ;;

        -addweight)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] && [ ! -z "$4" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/add_weight.py"
                python ${format_conv_path}/add_weight.py $2 $3 $4
            else
                echo " Usage: -addweight <input.xyz> <output.xyz> <weight> "
                echo " See the source code of add_groups.py for more details"
                echo " Code path: ${format_conv_path}/add_weight.py"
            fi ;;

        -get_frame)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/get_frame.py"
                python ${format_conv_path}/get_frame.py $2 $3
            else
                echo " Usage: -get_frame <exyzfile> <frame_index>"
                echo " See the source code of get_frame.py for more details"
                echo " Code path: ${format_conv_path}/get_frame.py"
            fi ;;

        -clean_xyz)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${format_conv_path}/clean_xyz.py"
                python ${format_conv_path}/clean_xyz.py $2 $3
            else
                echo " Usage: -clean_xyz <input.xyz> <output.xyz>"
                echo " See the source code of clean_xyz.py for more details"
                echo " Code path: ${format_conv_path}/clean_xyz.py"
            fi ;;

        -min_dist)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/get_min_dist.py"
                python ${analyzer_path}/get_min_dist.py $2
            else
                echo " Usage: -min_dist <exyzfile>"
                echo " See the source code of get_min_dist.py for more details"
                echo " Code path: ${analyzer_path}/get_min_dist.py"
            fi ;;

        -min_dist_pbc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/get_min_dist_pbc.py"
                python ${analyzer_path}/get_min_dist_pbc.py $2
            else
                echo " Usage: -min_dist_pbc <exyzfile>"
                echo " See the source code of get_min_dist_pbc.py for more details"
                echo " Code path: ${analyzer_path}/get_min_dist_pbc.py"
            fi ;;

        -filter_dist)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/filter_structures_by_distance.py"
                python ${analyzer_path}/filter_structures_by_distance.py $2 $3
            else
                echo " Usage: -filter_xyz <exyzfile> <min_dist>"
                echo " See the source code of filter_structures_by_distance.py for more details"
                echo " Code path: ${analyzer_path}/filter_structures_by_distance.py"
            fi ;;

        -filter_dist_pbc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/filter_structures_by_distance_pbc.py"
                python ${analyzer_path}/filter_structures_by_distance_pbc.py $2 $3
            else
                echo " Usage: -filter_xyz_pbc <exyzfile> <min_dist>"
                echo " See the source code of filter_structures_by_distance_pbc.py for more details"
                echo " Code path: ${analyzer_path}/filter_structures_by_distance_pbc.py"
            fi ;;

        -filter_box)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/filter_exyz_by_box.py"
                python ${analyzer_path}/filter_exyz_by_box.py $2 $3
            else
                echo " Usage: -filter_box <exyzfile> <lattice limit>"
                echo " See the source code of filter_exyz_by_box.py for more details"
                echo " Code path: ${analyzer_path}/filter_exyz_by_box.py"
            fi ;;

        -filter_value)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/filter_exyz_by_value.py"
                python ${analyzer_path}/filter_exyz_by_value.py $2 $3 $4
            else
                echo " Usage: -filter_value <exyzfile> <property> <value>"
                echo " See the source code of filter_exyz_by_value.py for more details"
                echo " Code path: ${analyzer_path}/filter_exyz_by_value.py"
            fi ;;

        -filter_range)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/filter_dist_range.py"
                python ${analyzer_path}/filter_dist_range.py $2 $3 $4 $5 $6
            else
                echo " Usage: -filter_range <exyzfile> <element1> <element2> <min_dist> <max_dist>"
                echo " See the source code of filter_dist_range.py for more details"
                echo " Code path: ${analyzer_path}/filter_dist_range.py"
            fi ;;

        -analyze_comp)
            if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/analyze_composition.py"
                python ${analyzer_path}/analyze_composition.py $2
            else
                echo " Usage: -analyze_composition <exyzfile>"
                echo " See the source code of analyze_composition.py for more details"
                echo " Code path: ${analyzer_path}/analyze_composition.py"
                exit 1
            fi ;;

        -get_volume)
            python ${analyzer_path}/get_volume.py
            ;;

        -pynep)
            parallel_pynep_sample_structures ;;

        -frame_range)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] && [ ! -z "$4" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${sample_path}/frame_range.py"
                python ${sample_path}/frame_range.py $2 $3 $4
            else
                echo " Usage: -frame_range <exyzfile> <start_frac> <end_frac>"
                echo " Examp: gpumdkit.sh -frame_range dump.xyz 0.2 0.5"
                echo " See the source code of frame_range.py for more details"
                echo " Code path: ${sample_path}/frame_range.py"
            fi ;;

        -re_atoms)
            echo " Calling script by Dian HUANG et al. "
            python ${utils_path}/renumber_atoms.py $2 $3 ;;

        -cbc)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] ; then
                echo " Calling script by Zihan YAN "
                python ${analyzer_path}/charge_balance_check.py $2
            else
                echo " Usage: -cbc <exyzfile>"
                echo " See the source code of charge_balance_check.py for more details"
                echo " Code path: ${analyzer_path}/charge_balance_check.py"
            fi ;;

        -pda)
            if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ] && [ ! -z "$4" ] && [ ! -z "$5" ] ; then
                echo " Calling script by Zihan YAN "
                echo " Code path: ${analyzer_path}/probability_density_analysis.py"
                python ${analyzer_path}/probability_density_analysis.py $2 $3 $4 $5
            else
                echo " Usage: -pda <ref_struct> <trajectory_file> <species> <interval>"
                echo " Examp: gpumdkit.sh -pda LLZO.vasp dump.xyz Li 0.25"
                echo " See the source code of probability_density_analysis.py for more details"
                echo " Code path: ${analyzer_path}/probability_density_analysis.py"
            fi ;;

        -hbond)
            echo " Calling script by Zherui CHEN "
            python ${GPUMD_path}/tools/Analysis_and_Processing/hydrogen_bond_analysis/Hydrogen-bond-analysis.py ${@:2}
            echo " Code path: ${GPUMD_path}/tools/Analysis_and_Processing/hydrogen_bond_analysis/Hydrogen-bond-analysis.py" ;;

        *)
            # echo " Unknown option: $1 "; help_info_table; exit 1 ;;
            alias_key="${1#-}"
            CUSTOM_CONFIG="$HOME/.gpumdkit.in"
            if [ ! -f "$CUSTOM_CONFIG" ]; then
                echo "Unknown: $1"
                echo "Use -h for help, or create $CUSTOM_CONFIG to define custom commands"
                exit 1
            fi
            func_name="custom_${alias_key}"
            source "$CUSTOM_CONFIG"
            if ! type "$func_name" >/dev/null 2>&1; then
                echo "Unknown: $1"
                echo "Tip: Define the function $func_name() { ... } in $CUSTOM_CONFIG to use this alias."
                help_info_table; exit 1
            fi
            shift
            "$func_name" "$@";;
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
  Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
      "
menu
main