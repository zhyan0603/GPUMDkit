#!/usr/bin/env bash

# Enable bash-style completion under zsh; bash is unaffected (block guarded by $ZSH_VERSION).
if [ -n "$ZSH_VERSION" ]; then
    case $- in
        *i*) ;;
        *) return 0 2>/dev/null || exit 0 ;;
    esac
    # Initialize zsh completion only if not already done (e.g. by oh-my-zsh)
    if ! typeset -p _comps >/dev/null 2>&1; then
        autoload -Uz +X compinit 2>/dev/null && compinit -i 2>/dev/null
    fi
    autoload -Uz +X bashcompinit 2>/dev/null && bashcompinit 2>/dev/null
fi

# Check if 'complete' and 'compgen' commands are available
if ! command -v complete >/dev/null 2>&1 || ! command -v compgen >/dev/null 2>&1; then
    # If either command is not found, exit silently
    return 0 2>/dev/null || exit 0
fi

# Bash completion function for gpumdkit.sh, supporting primary and secondary parameters
_gpumdkit_completions() {
    local cur prev candidate
    cur="${COMP_WORDS[COMP_CWORD]}"    # Current word being typed
    prev="${COMP_WORDS[COMP_CWORD-1]}" # Previous word

    # List of primary options (extracted from gpumdkit.sh)
    local opts="-h -help -skill -doctor -server -update -U -clean -time -plt -calc -prediction -prediction_dpa -cbc -range -shift_energy -out2xyz -out2xyz_bec -out2exyz -cp2k2xyz -pos2exyz -cif2pos -cif2exyz -exyz2pos -xdat2exyz -pos2lmp -lmp2exyz -traj2exyz -dp2xyz -xyz2dp -addgroup -addlabel -addweight -min_dist -min_dist_pbc -filter_dist -filter_dist_pbc -filter_box -filter_value -filter_range -get_frame -clean_xyz -analyze_comp -chem_species -replicate -pda -pynep -frame_range -nep_modifier"

    # Calculator subcommands
    local calc_subcmds="ionic-cond nep des doas neb minimize msd nlist disp avg-struct oct-tilt pol-abo3"

    # The first word after gpumdkit.sh is always a primary option
    if [ "$COMP_CWORD" -le 1 ]; then
        COMPREPLY=($(compgen -W "$opts" -- "$cur"))
        return
    fi

    # Later words: a secondary keyword right after -calc, -time or -plt,
    # otherwise a file or directory name (input files, output names, models)
    case "$prev" in
        -calc)
            COMPREPLY=($(compgen -W "$calc_subcmds" -- "$cur")) ;;

        -time)
            COMPREPLY=($(compgen -W "gpumd wpe nep gnep" -- "$cur")) ;;

        -plt)
            COMPREPLY=($(compgen -W "thermo thermo2 train train_density prediction test train_test parity_density born_charge bec msd msd_all msd_conv msd_sdc sdc rdf rdf_pmf xrd xrd_comp vac restart dimer force_errors des doas charge arrhenius_d arrhenius_sigma sigma D D_PT sigma_PT sigma_xyz D_xyz net_force emd emd2 nemd hnemd pdos phonon phonon_comp plane-grid cohesive viscosity" -- "$cur")) ;;

        *)
            # Mark the matches as file names so directories get a trailing slash
            # (compopt is missing in bash 3.2; completion still works without it)
            compopt -o filenames 2>/dev/null
            COMPREPLY=()
            while IFS= read -r candidate; do
                COMPREPLY+=("$candidate")
            done < <(compgen -f -- "$cur") ;;
    esac
}

# Register the completion function for gpumdkit.sh and common invocations
complete -F _gpumdkit_completions gpumdkit.sh
complete -F _gpumdkit_completions ./gpumdkit.sh
complete -F _gpumdkit_completions gpumdkit
