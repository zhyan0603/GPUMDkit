#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     out2xyz_bec.sh
# Category:   Format Conversion Scripts
# Purpose:    Convert VASP OUTCAR files to extxyz with Born effective charges,
#             energy, forces, and optional virial information.
# Usage:      gpumdkit.sh -out2xyz_bec <directory>
#             bash out2xyz_bec.sh <directory>
# Arguments:
#   directory  Root directory to search recursively for OUTCAR files
# Output:
#   NEPdataset/train.xyz  (one final configuration per OUTCAR)
# Author:     Yanzhou WANG, Shunda CHEN; Benrui TANG
# Last-modified: 2026-09-05
# =============================================================================

show_usage() {
    echo " Usage: gpumdkit.sh -out2xyz_bec <directory>"
    echo "    or: bash out2xyz_bec.sh <directory>"
    echo ""
    echo " Arguments:"
    echo "   directory   Root directory to search recursively for OUTCAR files"
    echo ""
    echo " Output:"
    echo "   NEPdataset/train.xyz   Extxyz data with optional bec:R:9 columns"
    echo ""
    echo " Notes:"
    echo "   One final configuration is written for each OUTCAR."
    echo "   An existing NEPdataset directory in the current directory is replaced."
    echo ""
}

if [ "$#" -lt 1 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_usage
    if [ "$#" -gt 0 ]; then
        exit 0
    fi
    exit 1
fi

read_dire="$1"
if [ ! -d "$read_dire" ]; then
    echo " Error: input directory '$read_dire' does not exist."
    exit 1
fi

writ_dire="NEPdataset"
writ_file="train.xyz"
isol_ener=0

list_file=""
temp_dir=""
cleanup() {
    if [ -n "$list_file" ]; then
        rm -f "$list_file"
    fi
    if [ -n "$temp_dir" ]; then
        rm -rf "$temp_dir"
    fi
}
trap cleanup EXIT

list_file=$(mktemp "${TMPDIR:-/tmp}/gpumdkit-out2xyz-bec.XXXXXX")
if [ ! -f "$list_file" ]; then
    echo " Error: failed to create a temporary OUTCAR list."
    exit 1
fi
if ! find -L "$read_dire" -type f -name "OUTCAR" -print | LC_ALL=C sort > "$list_file"; then
    echo " Error: failed to search for OUTCAR files under '$read_dire'."
    exit 1
fi

total_outcar=$(wc -l < "$list_file" | tr -d '[:space:]')
if [ "$total_outcar" -eq 0 ]; then
    echo " Error: no OUTCAR files were found under '$read_dire'."
    exit 1
fi

temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/gpumdkit-out2xyz-bec.XXXXXX")
if [ ! -d "$temp_dir" ]; then
    echo " Error: failed to create a temporary conversion directory."
    exit 1
fi
symbol_file="$temp_dir/symbols"
force_file="$temp_dir/forces"
bec_file="$temp_dir/bec"

rm -rf "$writ_dire"
if ! mkdir -p "$writ_dire"; then
    echo " Error: failed to create output directory '$writ_dire'."
    exit 1
fi
output_file="$writ_dire/$writ_file"
if ! : > "$output_file"; then
    echo " Error: failed to create output file '$output_file'."
    exit 1
fi

echo " Found $total_outcar OUTCAR file(s)."
echo " Converting the final configuration from each OUTCAR ..."

processed=0
frames_written=0
bec_configurations=0

while IFS= read -r file; do
    if [ ! -r "$file" ]; then
        echo " Error: OUTCAR '$file' is not readable."
        exit 1
    fi

    syst_numb_atom=$(awk '/number of ions/ {value=$12} END {print value}' "$file")
    case "$syst_numb_atom" in
        ''|*[!0-9]*)
            echo " Error: failed to parse the atom count from '$file'."
            exit 1
            ;;
    esac
    if [ "$syst_numb_atom" -le 0 ]; then
        echo " Error: the atom count in '$file' is not positive."
        exit 1
    fi

    latt=$(awk '
        /VOLUME and BASIS-vectors are now/ {
            in_cell=1
            vector_count=0
            lattice=""
            next
        }
        in_cell && $1 ~ /^[-+0-9.]+([eE][-+0-9]+)?$/ && \
            $2 ~ /^[-+0-9.]+([eE][-+0-9]+)?$/ && \
            $3 ~ /^[-+0-9.]+([eE][-+0-9]+)?$/ {
            vector[vector_count++]=$1 " " $2 " " $3
            if (vector_count == 3) {
                lattice=vector[0] " " vector[1] " " vector[2]
                in_cell=0
            }
        }
        END {
            if (lattice != "") print lattice
        }
    ' "$file" | tail -n 1)
    if [ -z "$latt" ]; then
        echo " Error: failed to parse lattice vectors from '$file'."
        exit 1
    fi

    ener=$(grep -F "free  energy   TOTEN" "$file" | tail -n 1 | awk -v atoms="$syst_numb_atom" -v shift="$isol_ener" '{printf "%.10f\n", $5 - atoms * shift}')
    if [ -z "$ener" ]; then
        echo " Error: failed to parse the final energy from '$file'."
        exit 1
    fi

    force_header=$(grep -n -F "TOTAL-FORCE (eV/Angst)" "$file" | tail -n 1 | cut -d: -f1)
    if [ -z "$force_header" ]; then
        echo " Error: no final force block was found in '$file'."
        exit 1
    fi
    sed -n "$((force_header + 2)),$((force_header + syst_numb_atom + 1))p" "$file" > "$force_file"
    force_rows=$(wc -l < "$force_file" | tr -d '[:space:]')
    if [ "$force_rows" -ne "$syst_numb_atom" ] || ! awk '
        $1 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ || \
        $2 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ || \
        $3 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ || \
        $4 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ || \
        $5 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ || \
        $6 !~ /^[-+0-9.]+([eE][-+0-9]+)?$/ {bad=1}
        END {exit bad}
    ' "$force_file"; then
        echo " Error: the final force block in '$file' is incomplete or malformed."
        exit 1
    fi

    ion_numb_arra=($(awk -F= '/ions per type/ {value=$2} END {print value}' "$file"))
    ion_symb_arra=($(awk '
        /TITEL/ {
            symbol=$4
            sub(/_.*/, "", symbol)
            if (symbol != "" && !seen[symbol]++) print symbol
        }
    ' "$file"))
    if [ "${#ion_numb_arra[@]}" -eq 0 ] || [ "${#ion_numb_arra[@]}" -ne "${#ion_symb_arra[@]}" ]; then
        echo " Error: failed to parse element types from '$file'."
        exit 1
    fi

    : > "$symbol_file"
    species_total=0
    for ((j=0; j<${#ion_numb_arra[@]}; j++)); do
        case "${ion_numb_arra[$j]}" in
            ''|*[!0-9]*)
                echo " Error: failed to parse element counts from '$file'."
                exit 1
                ;;
        esac
        species_total=$((species_total + ion_numb_arra[j]))
        for ((k=0; k<ion_numb_arra[j]; k++)); do
            printf '%s\n' "${ion_symb_arra[$j]}" >> "$symbol_file"
        done
    done
    if [ "$species_total" -ne "$syst_numb_atom" ]; then
        echo " Error: element counts do not match NIONS in '$file'."
        exit 1
    fi

    has_bec=0
    : > "$bec_file"
    if grep -q -F "BORN EFFECTIVE CHARGES" "$file"; then
        has_bec=1
        if ! awk -v atoms="$syst_numb_atom" '
            /BORN EFFECTIVE CHARGES/ {
                in_bec=1
                value_count=0
                next
            }
            in_bec && value_count < atoms * 9 && $1 ~ /^[123]$/ && NF >= 4 {
                for (column=2; column<=4; column++) values[++value_count]=$column
            }
            END {
                if (value_count != atoms * 9) exit 1
                for (i=1; i<=value_count; i++) {
                    printf "%15.6f", values[i]
                    if (i % 9 == 0) printf "\n"
                }
            }
        ' "$file" > "$bec_file"; then
            echo " Error: the BEC block in '$file' is incomplete or malformed."
            exit 1
        fi
        bec_rows=$(wc -l < "$bec_file" | tr -d '[:space:]')
        if [ "$bec_rows" -ne "$syst_numb_atom" ]; then
            echo " Error: BEC rows do not match NIONS in '$file'."
            exit 1
        fi
    fi

    configuration=$(basename "$(dirname "$file")")
    isif=$(awk '/ISIF/ {for (i=1; i<=NF; i++) if ($i == "=") value=$(i+1)} END {print value}' "$file")
    viri=""
    case "$isif" in
        ''|*[!0-9]*) ;;
        *)
            if [ "$isif" -ne 0 ]; then
                viri=$(grep -A 20 -F "FORCE on cell =-STRESS" "$file" | grep -F "Total " | tail -n 1 | awk '{print $2,$5,$7,$5,$3,$6,$7,$6,$4}')
            fi
            ;;
    esac

    if [ "$has_bec" -eq 1 ]; then
        properties="species:S:1:pos:R:3:forces:R:3:bec:R:9"
    else
        properties="species:S:1:pos:R:3:forces:R:3"
    fi

    printf '%s\n' "$syst_numb_atom" >> "$output_file"
    if [ -n "$viri" ]; then
        printf '%s\n' "Config_type=$configuration Weight=1.0 Lattice=\"$latt\" Energy=$ener Virial=\"$viri\" pbc=\"T T T\" Properties=$properties" >> "$output_file"
    else
        printf '%s\n' "Config_type=$configuration Weight=1.0 Lattice=\"$latt\" Energy=$ener pbc=\"T T T\" Properties=$properties" >> "$output_file"
    fi
    if [ "$has_bec" -eq 1 ]; then
        paste "$symbol_file" "$force_file" "$bec_file" >> "$output_file"
        bec_configurations=$((bec_configurations + 1))
    else
        paste "$symbol_file" "$force_file" >> "$output_file"
    fi

    processed=$((processed + 1))
    frames_written=$((frames_written + 1))
    progress=$((processed * 100 / total_outcar))
    printf " Progress: %d%% (%d/%d)\r" "$progress" "$processed" "$total_outcar"
done < "$list_file"

printf "\n Conversion complete: %d configuration(s) written to %s.\n" "$frames_written" "$output_file"
if [ "$bec_configurations" -lt "$frames_written" ]; then
    echo " Warning: $((frames_written - bec_configurations)) configuration(s) had no BEC block; those rows omit bec:R:9."
fi
