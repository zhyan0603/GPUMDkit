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
    local opts="-h -help -clean -time -plt -calc -range -out2xyz -outcar2exyz -cast2xyz -castep2exyz -cp2k2xyz -cp2k2exyz -mtp2xyz -mtp2exyz -pos2exyz -exyz2pos -pos2lmp -lmp2exyz -addgroup -addlabel -addweight -max_rmse -get_max_rmse_xyz -min_dist -min_dist_pbc -filter_dist -filter_dist_pbc -filter_box -filter_value -get_frame -clear_xyz -clean_xyz"

    # Provide secondary completion based on the previous word
    case "$prev" in
        # Secondary options for -time
        -time)
            COMPREPLY=($(compgen -W "gpumd nep" -- "$cur"))
            ;;

        # Secondary options for -plt
        -plt)
            COMPREPLY=($(compgen -W "thermo train prediction valid test train_test msd sdc rdf vac restart dimer force_error" -- "$cur"))
            ;;

        # Secondary options for -calc
        -calc)
            COMPREPLY=($(compgen -W "ionic-cond nep" -- "$cur"))
            ;;

        # Secondary options for -range
        -range)
            COMPREPLY=($(compgen -W "energy force virial" -- "$cur"))  # Example properties, assumed supported
            ;;

        # Options requiring files or directories, complete with filenames
        -out2xyz|-outcar2exyz|-cast2xyz|-castep2exyz|-cp2k2xyz|-cp2k2exyz|-exyz2pos|-min_dist|-min_dist_pbc|-filter_dist|-filter_dist_pbc|-filter_box|-get_frame|-clear_xyz|-clean_xyz|-mtp2xyz|-mtp2exyz|-pos2exyz|-pos2lmp|-lmp2exyz|-addgroup|-addlabel|-addweight|-max_rmse|-get_max_rmse_xyz)
            COMPREPLY=($(compgen -f -- "$cur"))  # Complete files or directories
            ;;

        # Secondary completion for -filter_value
        -filter_value)
            if [ "$COMP_CWORD" -eq 2 ]; then
                COMPREPLY=($(compgen -f -- "$cur"))  # Complete extxyz file
            else 
                COMPREPLY=($(compgen -W "energy force virial" -- "$cur"))  # Properties
            fi
            ;;

        # Default case: complete primary options
        *)
            COMPREPLY=($(compgen -W "$opts" -- "$cur"))
            ;;
    esac
}

# Register the completion function
complete -F _gpumdkit_completions gpumdkit.sh