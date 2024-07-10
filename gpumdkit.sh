#!/bin/bash

# You need to set the path of GPUMDtools
GPUMDtools_path=/d/Westlake/GPUMD/GPUMD/tools
VERSION="0.0.0 (dev) (2024-07-10)"


#--------------------- main script ----------------------
# Show the menu
function menu(){
echo " ------------------------- GPUMD --------------------------"
echo "  1) Developing ...              2) Developing ...         "
echo "  0) Quit!"
}

# Function main
function main(){
    echo " ------------>>"
    echo ' Input the function number:'
    array_choice=("0" "1" "2") 
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
            echo "Developing ..."
            ;;
        "2")
            echo "Developing ..."
            ;;
        
    esac
    echo " Thank you for using GPUMDkit. Have a great day!"
}


######### Custom functional area ###############
function help_info(){
echo " GPUMDkit ${VERSION}"
echo " Usage: GPUMDkit -[options]"
echo " Options:
    -outcar2exyz    Convert OUTCAR to nep-exyz file
                    Usage: -outcar2exyz dir_name
                      Examp: gpumdkit.sh -outcar2exyz .

    -castep2exyz    Convert castep to nep-exyz file
                    Usage: -castep2exyz dir_name
                      Examp: gpumdkit.sh -castep2exyz .

    -cp2k2exyz    Convert cp2k output to nep-exyz file
                    Usage: -cp2k2exyz dir_name
                      Examp: gpumdkit.sh -cp2k2exyz .

    -max_rmse         get_max_rmse_xyz
                    Usage: -getmax|-get_max_rmse_xyz train.xyz force_train.out 13

    -h,-help    Show this help message"
}

if [ ! -z "$1" ]; then
    case $1 in
        -h|-help)
            help_info
            ;;

        -out2xyz|-outcar2exyz)
            if [ ! -z "$2" ]; then
                echo "Calling script by Yanzhou WANG et al. "
                echo "Code path: ${GPUMDtools_path}/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
                bash ${GPUMDtools_path}/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh $2
            else
                echo "Missing argument"
                echo "Usage: -out2xyz|-outcar2exyz dir_name (eg. gpumdkit.sh -outcar2exyz .)"
                echo "See the source code of multipleFrames-outcars2nep-exyz.sh for more details"
                echo "Code path: ${GPUMDtools_path}/vasp2xyz/outcar2xyz/multipleFrames-outcars2nep-exyz.sh"
            fi
            ;;

        -cast2xyz|-castep2exyz)
            if [ ! -z "$2" ]; then
                echo "Calling script by Yanzhou WANG et al. "
                echo "Code path: ${GPUMDtools_path}/castep2exyz/castep2nep-exyz.sh"
                bash ${GPUMDtools_path}/castep2exyz/castep2nep-exyz.sh $2
            else
                echo "Missing argument"
                echo "Usage: -cast2xyz|-castep2exyz dir_name (eg. gpumdkit.sh -castep2exyz .)"
                echo "See the source code of castep2nep-exyz.sh for more details"
                echo "Code path: ${GPUMDtools_path}/castep2exyz/castep2nep-exyz.sh"
            fi
            ;;

        -cp2k2xyz|-cp2k2exyz)
            if [ ! -z "$2" ]; then
                echo "Calling script by Ke XU et al. "
                echo "Code path: ${GPUMDtools_path}/cp2k2xyz/cp2k2xyz.py"
                python ${GPUMDtools_path}/cp2k2xyz/cp2k2xyz.py $2
            else
                echo "Missing argument"
                echo "Usage: -cp2k2xyz|-cp2k2exyz dir_name (eg. gpumdkit.sh -cp2k2exyz .)"
                echo "See the source code of cp2k2xyz.py for more details"
                echo "Code path: ${GPUMDtools_path}/cp2k2xyz/cp2k2xyz.py"
            fi
            ;;

        -mtp2xyz|-mtp2exyz)
            if [ ! -z "$2" ]; then
                echo "Calling script by Ke XU et al. "
                echo "Code path: ${GPUMDtools_path}/mtp2xyz/mtp2xyz.py"
                python ${GPUMDtools_path}/mtp2xyz/mtp2xyz.py train.cfg $2 ${@:3}
            else
                echo "Missing argument"
                echo "Usage: -mtp2xyz|-mtp2exyz train.cfg Symbol1 Symbol2 Symbol3 ..."
                echo "  Examp: gpumdkit.sh -mtp2exyz train.cfg Pd Ag"
                echo "See the source code of mtp2xyz.py for more details"
                echo "Code path: ${GPUMDtools_path}/mtp2xyz/mtp2xyz.py"
            fi
            ;;

        -max_rmse|-get_max_rmse_xyz)
            if [ ! -z "$2" ] && [ ! -z "$3" ] && [ ! -z "$4" ]; then
                echo "Calling script by Ke XU "
                echo "Code path: ${GPUMDtools_path}/get_max_rmse_xyz/get_max_rmse_xyz.py"
                python ${GPUMDtools_path}/get_max_rmse_xyz/get_max_rmse_xyz.py $2 $3 $4
            else
                echo "Missing argument"
                echo "Usage: -getmax|-get_max_rmse_xyz train.xyz force_train.out 13"
                echo "See the source code of get_max_rmse_xyz.py for more details"
                echo "Code path: ${GPUMDtools_path}/get_max_rmse_xyz/get_max_rmse_xyz.py"
            fi
            ;;

        *)
            echo " Unknown option: $1 "
            help_info
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