#!/bin/bash

# GPUMDkit Help System
# This file contains all help information to keep main script clean

# Get the version from main script or set default
GPUMDKIT_VERSION="${VERSION:-1.5.3 (dev)}"

# Terminal width for formatting
TERM_WIDTH=${COLUMNS:-80}

# Calculate dynamic padding for centering
center_text() {
    local text="$1"
    local width="$2"
    local text_len=${#text}
    local padding=$(( (width - text_len) / 2 ))
    printf "%*s%s\n" $padding "" "$text"
}

# Print a horizontal line
print_line() {
    local char="${1:-=}"
    local width="${2:-$TERM_WIDTH}"
    printf '%*s\n' "$width" | tr ' ' "$char"
}

# Print a section header
print_section_header() {
    local title="$1"
    local width="${2:-$TERM_WIDTH}"
    local char="${3:-=}"
    
    local title_len=${#title}
    local total_decor_width=$((width - 2))
    local side_padding=$(( (total_decor_width - title_len) / 2 ))
    
    printf "+%s+\n" "$(printf '%*s' "$((width-2))" | tr ' ' "$char")"
    printf "|%*s%*s%*s|\n" "$side_padding" "" "$title_len" "$title" "$side_padding" ""
    printf "+%s+\n" "$(printf '%*s' "$((width-2))" | tr ' ' "$char")"
}

# Main help table function
help_info_table() {
    local width=100
    
    echo "+==================================================================================================+"
    echo "|                              GPUMDkit ${GPUMDKIT_VERSION} Usage                             |"
    echo "+======================================== Conversions =============================================+"
    echo "| -out2xyz       Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           |"
    echo "| -cif2pos       Convert cif to POSCAR          | -pos2lmp      Convert POSCAR to LAMMPS           |"
    echo "| -cif2exyz      Convert cif to extxyz          | -lmp2exyz     Convert LAMMPS-dump to extxyz      |"
    echo "| -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz |"
    echo "| -cp2k2xyz      Convert CP2K file to extxyz    | -traj2exyz    Convert ASE traj to extxyz         |"
    echo "| -xdat2exyz     Convert XDATCAR to extxyz      | -abacus2xyz   Convert ABACUS to extxyz           |"
    echo "| -extxyz2pos    Convert extxyz to POSCAR       | -mtp2xyz      Convert MTP to extxyz              |"
    echo "+========================================= Analysis ===============================================+"
    echo "| -range         Print range of energy etc.     | -max_rmse     Get max RMSE from extxyz           |"
    echo "| -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       |"
    echo "| -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       |"
    echo "| -filter_dist   Filter struct by min_dist      | -analyze_comp Analyze composition of extxyz      |"
    echo "| -pynep         Sample struct by pynep         | -chem_species Analyze chemical species           |"
    echo "| -cbc           Charge balance check           | -pda          Probability density analysis       |"
    echo "+====================================== Misc Utilities ============================================+"
    echo "| -plt           Plot scripts                   | -get_frame     Extract the specified frame       |"
    echo "| -calc          Calculators                    | -frame_range   Extract frames by fraction range  |"
    echo "| -clean         Clear files for work_dir       | -clean_xyz     Clean extra info in XYZ file      |"
    echo "| -time          Time consuming Analyzer        | -update        Update GPUMDkit                   |"
    echo "| -replicate     Replicate structures           | -re_atoms      Renumber atoms                    |"
    echo "+========================================== Workflow ==============================================+"
    echo "| -scf_batch     SCF batch pretreatment         | -md_batch_gpumd MD batch pretreatment (GPUMD)   |"
    echo "| -md_batch_lmp  MD batch pretreatment (LAMMPS) | Developing...                                    |"
    echo "+==================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   |"
    echo "| Example: gpumdkit.sh -out2xyz -h                                                                   |"
    echo "+==================================================================================================+"
}

# Plotting help table function
plot_info_table() {
    echo "+=====================================================================================================+"
    echo "|                              GPUMDkit ${GPUMDKIT_VERSION} Plotting Usage                       |"
    echo "+=============================================== Plot Types ==========================================+"
    echo "| thermo          Plot thermo info                   | train          Plot NEP train results          |"
    echo "| prediction      Plot NEP prediction results        | train_test     Plot NEP train and test results |"
    echo "| msd             Plot mean square displacement      | msd_conv       Plot the convergence of MSD     |"
    echo "| msd_all         Plot MSD of all species            | sdc            Plot self diffusion coefficient |"
    echo "| msd_sdc         Plot MSD and SDC together          | cohesive       Plot cohsive energy             |"
    echo "| rdf             Plot radial distribution function  | vac            Plot velocity autocorrelation   |"
    echo "| restart         Plot parameters in nep.restart     | dimer          Plot dimer plot                 |"
    echo "| force_errors    Plot force errors                  | des            Plot descriptors                |"
    echo "| charge          Plot charge distribution           | lr             Plot learning rate              |"
    echo "| doas            Plot density of atomistic states   | net_force      Plot net force distribution     |"
    echo "| sigma           Plot Arrhenius sigma               | D              Plot Arrhenius diffusivity      |"
    echo "| sigma_xyz       Plot directional Arrhenius sigma   | D_xyz          Plot directional Arrhenius D    |"
    echo "| emd             Plot EMD results                   | nemd           Plot NEMD results               |"
    echo "| hnemd           Plot HNEMD results                 | pdos           Plot VAC and PDOS               |"
    echo "| plane-grid      Plot displacement plane grid       | parity_density Plot parity plot density        |"
    echo "| rdf_pmf         Plot potential of mean force (PMF) | viscosity      Plot visconsity                 |"
    echo "+=====================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -plt <plot_type> -h                               |"
    echo "| Example: gpumdkit.sh -plt thermo -h                                                                  |"
    echo "+=====================================================================================================+"
}

# Calculator help table function
calc_info_table() {
    echo "+=====================================================================================================+"
    echo "|                              GPUMDkit ${GPUMDKIT_VERSION} Calculators Usage                     |"
    echo "+============================================= Calculator Types =======================================+"
    echo "| ionic-cond      Calculate ionic conductivity       | nep            Calculate properties with NEP   |"
    echo "| des             Calculate descriptors              | doas           Calculate DOAS                  |"
    echo "| nlist           Calculate neighbor list            | disp           Calculate displacement          |"
    echo "| avg-struct      Calculate averaged structure       | oct-tilt       Calculate octahedral tilt       |"
    echo "| pol-abo3        Calculate polarization (ABO3)      | Developing...                                  |"
    echo "+=====================================================================================================+"
    echo "| For detailed usage and examples, use: gpumdkit.sh -calc <calc_type> -h                              |"
    echo "| Example: gpumdkit.sh -calc ionic-cond -h                                                             |"
    echo "+=====================================================================================================+"
}

# Function-specific help messages (can be extended)
show_option_help() {
    local option="$1"
    
    case "$option" in
        -out2xyz)
            echo "Usage: gpumdkit.sh -out2xyz <directory>"
            echo ""
            echo "Description:"
            echo "  Convert VASP OUTCAR or vasprun.xml files to extended XYZ format."
            echo ""
            echo "Arguments:"
            echo "  directory    Path to directory containing OUTCAR or vasprun.xml files"
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -out2xyz ."
            echo "  gpumdkit.sh -out2xyz ./vasp_runs/"
            ;;
        -plt)
            echo "Usage: gpumdkit.sh -plt <plot_type> [arguments]"
            echo ""
            echo "Description:"
            echo "  Generate various plots from simulation data."
            echo ""
            echo "Arguments:"
            echo "  plot_type    Type of plot to generate (see -plt -h for full list)"
            echo "  arguments    Additional arguments specific to the plot type"
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -plt thermo log.txt"
            echo "  gpumdkit.sh -plt msd msd.out"
            echo "  gpumdkit.sh -plt -h    # Show all available plot types"
            ;;
        -calc)
            echo "Usage: gpumdkit.sh -calc <calculator_type> [arguments]"
            echo ""
            echo "Description:"
            echo "  Perform various calculations on simulation data."
            echo ""
            echo "Arguments:"
            echo "  calculator_type  Type of calculation (see -calc -h for full list)"
            echo "  arguments        Additional arguments specific to the calculation"
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -calc ionic-cond Li 1"
            echo "  gpumdkit.sh -calc nep input.xyz output.xyz nep.txt"
            echo "  gpumdkit.sh -calc -h    # Show all available calculators"
            ;;
        -range)
            echo "Usage: gpumdkit.sh -range <file> <property> [hist]"
            echo ""
            echo "Description:"
            echo "  Print the range of energy, force, or virial from an extended XYZ file."
            echo ""
            echo "Arguments:"
            echo "  file         Path to extended XYZ file"
            echo "  property     Property to analyze: energy, force, or virial"
            echo "  hist         Optional: show histogram if specified"
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -range train.xyz energy"
            echo "  gpumdkit.sh -range train.xyz force hist"
            ;;
        -clean)
            echo "Usage: gpumdkit.sh -clean"
            echo ""
            echo "Description:"
            echo "  Clean up unnecessary files from the working directory."
            echo "  Removes temporary files, backup files, and other clutter."
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -clean"
            ;;
        -update|-U)
            echo "Usage: gpumdkit.sh -update"
            echo "       gpumdkit.sh -U"
            echo ""
            echo "Description:"
            echo "  Update GPUMDkit to the latest version from the repository."
            echo ""
            echo "Examples:"
            echo "  gpumdkit.sh -update"
            echo "  gpumdkit.sh -U"
            ;;
        *)
            return 1
            ;;
    esac
    return 0
}

# Interactive menu help (for submenu contexts)
show_menu_help() {
    local menu_type="$1"
    
    case "$menu_type" in
        main)
            echo "Main Menu Options:"
            echo "  1) Format Conversion    - Convert between different file formats"
            echo "  2) Sample Structures    - Generate and sample atomic structures"
            echo "  3) Workflow             - Batch processing workflows"
            echo "  4) Calculators          - Various property calculators"
            echo "  5) Analyzer             - Analyze simulation results"
            echo "  6) Developing           - Under development features"
            echo "  0) Quit                 - Exit GPUMDkit"
            ;;
        format_conversion)
            echo "Format Conversion Options:"
            echo "  101) OUTCAR to extxyz           - Convert VASP OUTCAR files"
            echo "  102) MTP to extxyz              - Convert MTP training files"
            echo "  103) CP2K to extxyz             - Convert CP2K output files"
            echo "  104) ABACUS to extxyz           - Convert ABACUS output files"
            echo "  105) extxyz to POSCAR           - Convert extended XYZ to VASP"
            echo "  000) Return to main menu"
            ;;
        *)
            echo "Help not available for this menu type."
            ;;
    esac
}

# Citation message
citation() {
    echo " +------------------------------------------------------+"
    echo " |           THANK YOU FOR USING GPUMDKIT               |"
    echo " +------------------------------------------------------+"
    echo " | If you find it useful, please cite our paper:        |"
    echo " |                                                      |"
    echo " | GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP  |"
    echo " |       (https://arxiv.org/abs/2603.17367)             |"
    echo " |                                                      |"
    echo " |     Welcome to join our QQ group (825696376) !       |"
    echo " +------------------------------------------------------+"
}

# Print GPUMDkit logo/banner
print_banner() {
    echo ""
    echo "           ____ ____  _   _ __  __ ____  _    _ _   "
    echo "          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_ "
    echo "         | |  _| |_) | | | | |\/| | | | | |/ / | __|"
    echo "         | |_| |  __/| |_| | |  | | |_| |   <| | |_ "
    echo "          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|"
    echo "                                           "
    echo "          GPUMDkit Version ${GPUMDKIT_VERSION}"
    echo "    Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)"
    echo " Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA"
    echo ""
}
