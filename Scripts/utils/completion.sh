#!/usr/bin/env bash

# Check if 'complete' and 'compgen' commands are available
if ! command -v complete >/dev/null 2>&1 || ! command -v compgen >/dev/null 2>&1; then
    # If either command is not found, exit silently
    return 0 2>/dev/null || exit 0
fi

# Bash completion function for gpumdkit.sh, supporting primary and secondary parameters
_gpumdkit_completions() {
    local cur prev
    cur="${COMP_WORDS[COMP_CWORD]}"    # Current word being typed
    prev="${COMP_WORDS[COMP_CWORD-1]}" # Previous word

    # List of primary options (extracted from gpumdkit.sh)
    local opts="-h -update -help -clean -time -plt -calc -cbc -range -out2xyz -xml2xyz -pos2exyz -cif2pos -cif2exyz -exyz2pos -pos2lmp -lmp2exyz -cp2k2xyz -addgroup -addlabel -addweight -min_dist -min_dist_pbc -filter_dist -filter_dist_pbc -filter_box -filter_value -get_frame -clean_xyz -get_volume -analyze_comp -replicate -hbond -pda -pynep -frame_range -chem_species"

    # Provide secondary completion based on the previous word
    case "$prev" in
        # Secondary options for -time
        -time)
            COMPREPLY=($(compgen -W "gpumd nep gnep" -- "$cur")) ;;

        # Secondary options for -plt
        -plt)
            COMPREPLY=($(compgen -W "thermo thermo2 thermo3 train prediction test train_test msd msd_all msd_conv sdc rdf vac restart dimer force_error des doas charge lr parity_density arrhenius_d arrhenius_sigma sigma D sigma_xyz D_xyz net_force emd nemd hnemd pdos" -- "$cur")) ;;

        # Secondary options for -calc
        -calc)
            COMPREPLY=($(compgen -W "ionic-cond nep des doas neb" -- "$cur")) ;;
            
        # Options requiring files or directories, complete with filenames
        -out2xyz|-cp2k2xyz|-exyz2pos|-min_dist|-min_dist_pbc|-filter_dist|-filter_dist_pbc|-filter_box|-get_frame|-clean_xyz|-mtp2xyz|-pos2exyz|-cif2exyz|-cif2pos|-pos2lmp|-lmp2exyz|-addgroup|-addlabel|-addweight|-analyze_comp|-replicate|-pda|-cbc|-frame_range|-xml2xyz|-chem_species)
            COMPREPLY=($(compgen -f -- "$cur")) ;;

        # Default case: complete primary options
        *)
            COMPREPLY=($(compgen -W "$opts" -- "$cur")) ;;
    esac
}

# Register the completion function
complete -F _gpumdkit_completions gpumdkit.sh