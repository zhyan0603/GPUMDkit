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
    local opts="-h -help -skill -doctor -update -U -clean -time -plt -calc -cbc -range -out2xyz -out2exyz -cp2k2xyz -pos2exyz -cif2pos -cif2exyz -exyz2pos -xdat2exyz -pos2lmp -lmp2exyz -traj2exyz -dp2xyz -addgroup -addlabel -addweight -min_dist -min_dist_pbc -filter_dist -filter_dist_pbc -filter_box -filter_value -filter_range -get_frame -clean_xyz -analyze_comp -chem_species -replicate -pda -pynep -frame_range -nep_modifier"

    # Calculator subcommands
    local calc_subcmds="ionic-cond nep des doas neb minimize msd nlist disp avg-struct oct-tilt pol-abo3"

    # If we are completing a subcommand right after -calc
    if [[ "$prev" == "-calc" ]]; then
        COMPREPLY=($(compgen -W "$calc_subcmds" -- "$cur"))
        return
    fi

    # Provide secondary completion based on the previous word
    case "$prev" in
        # Secondary options for -time
        -time)
            COMPREPLY=($(compgen -W "gpumd nep gnep" -- "$cur")) ;;

        # Secondary options for -plt
        -plt)
            COMPREPLY=($(compgen -W "thermo thermo2 thermo3 train train_density prediction test train_test parity_density born_charge bec msd msd_all msd_conv msd_sdc sdc rdf rdf_pmf vac restart dimer force_errors des doas charge lr arrhenius_d arrhenius_sigma sigma D sigma_xyz D_xyz net_force emd nemd hnemd pdos plane-grid cohesive viscosity" -- "$cur")) ;;

        # Options requiring files or directories, complete with filenames
        -out2xyz|-out2exyz|-cp2k2xyz|-exyz2pos|-min_dist|-min_dist_pbc|-filter_dist|-filter_dist_pbc|-filter_box|-filter_value|-filter_range|-get_frame|-clean_xyz|-pos2exyz|-cif2exyz|-cif2pos|-pos2lmp|-lmp2exyz|-traj2exyz|-dp2xyz|-addgroup|-addlabel|-addweight|-analyze_comp|-replicate|-pda|-cbc|-frame_range|-chem_species|-xdat2exyz)
            COMPREPLY=($(compgen -f -- "$cur")) ;;

        # Default case: complete primary options
        *)
            COMPREPLY=($(compgen -W "$opts" -- "$cur")) ;;
    esac
}

# Register the completion function for gpumdkit.sh and common invocations
complete -F _gpumdkit_completions gpumdkit.sh
complete -F _gpumdkit_completions ./gpumdkit.sh
complete -F _gpumdkit_completions gpumdkit
