#!/bin/bash
# =============================================================================
# GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
# Repository: https://github.com/zhyan0603/GPUMDkit
# Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
#           MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
# =============================================================================
# Script:     abacus2xyz_md.sh
# Category:   Format Conversion Scripts
# Purpose:    Convert ABACUS MD trajectory files (running_md.log and MD_dump)
#             to extended XYZ format for NEP training.
# Usage:      ./abacus2xyz_md.sh <dire_name>
# Arguments:
#   dire_name  Directory containing running_md.log and MD_dump files
# Output:
#   NEPdataset/train.xyz  (converted dataset in extxyz format)
# Author:     Benrui Tang (tang070205@proton.me)
# Last-modified: 2026-09-14
# =============================================================================

#--- DEFAULT ASSIGNMENTS ---------------------------------------------------------------------
isol_ener=0     # Shifted energy, specify the value?
viri_logi=1     # Logical value for virial, true=1, false=0
#--------------------------------------------------------------------------------------------
read_dire=$1
if [ -z "$read_dire" ]; then
        echo " Your syntax is illegal, please try again"
        exit 1
fi

md_log_file=$(find -L "$read_dire" -name "running_md.log" 2>/dev/null | head -n 1)
md_dump_file=$(find -L "$read_dire" -name "MD_dump" 2>/dev/null | head -n 1)
input_file=$(find -L "$read_dire" -name "INPUT" 2>/dev/null | head -n 1)
if [ -z "$md_log_file" ] || [ -z "$md_dump_file" ] || [ -z "$input_file" ]; then
        echo " Error: could not find running_md.log, MD_dump, and INPUT under '$read_dire'."
        echo " Please check the directory."
        exit 1
fi

writ_dire="NEPdataset"; writ_file="train.xyz";
rm -rf "$writ_dire"; mkdir "$writ_dire"

configuration=$(cd "$(dirname "$md_log_file")" && pwd | awk -F'/' '{print $(NF-2)"/"$(NF-1)"/"$NF}')
syst_numb_atom=$(grep "TOTAL ATOM NUMBER" "$md_log_file" |awk '{print $5}')
ener_values=($(grep 'etot' "$md_log_file" |awk '{print $4}'))
scf_lines=($(grep -n 'STEP OF MOLECULAR DYNAMICS' "$md_log_file" | awk -F: '{print $1}'))
scf_last=$(wc -l < "$md_log_file")
scf_lines+=($scf_last)
scf_nmax=$(grep 'scf_nmax' "$input_file" | awk '{print $2}')
mdstep_lines=($(grep -n 'MDSTEP' "$md_dump_file" | awk -F: '{print $1}'))
mdstep_last=$(wc -l < "$md_dump_file")
mdstep_lines+=($mdstep_last)
N_counts=$(( ${#mdstep_lines[@]} - 2 ))


for ((i=1; i<$(( ${#mdstep_lines[@]} - 1 )); i++)); do
    ener=${ener_values[i]}

    scf_start=${scf_lines[i]}
    scf_end=${scf_lines[i+1]}
    scf_act=$(sed -n "${scf_start},${scf_end}p" "$md_log_file" | grep -c "ALGORITHM")
    if [ "$scf_act" -eq "$scf_nmax" ]; then
        echo " Skipping the $i structure due to non convergence"
        echo -ne " Process: ${i}/${N_counts}\r"
        continue
    fi

    md_start=${mdstep_lines[i]}
    md_end=${mdstep_lines[i+1]}
    sed -n "${md_start},${md_end}p" "$md_dump_file" > temp.file
    echo "$syst_numb_atom" >> "$writ_dire/$writ_file"
    latt=$(grep -A 3 "LATTICE_VECTORS" temp.file | tail -n 3 | awk '{for (i = 1; i <= NF; i++) {printf "%.8f ", $i}}' |xargs)
    conversion_value=$(echo "$latt" | awk '{a1=$1; a2=$2; a3=$3; b1=$4; b2=$5; b3=$6; c1=$7; c2=$8; c3=$9;
        V=a1*(b2*c3 - b3*c2) + a2*(b3*c1 - b1*c3) + a3*(b1*c2 - b2*c1); if (V < 0) V=-V; printf "%.8f", V/1602.1766208}')
    if [[ $viri_logi -eq 1 ]]; then
        viri=$(grep -A 3 "VIRIAL (kbar)" temp.file | tail -n 3 | awk '{for (i = 1; i <= NF; i++) {printf "%.8f ", $i * '$conversion_value'}}' |xargs)
        echo "Energy=$ener Lattice=\"$latt\" Virial=\"$viri\" Config_type=$configuration-$i Weight=1.0 Properties=species:S:1:pos:R:3:forces:R:3" >> "$writ_dire/$writ_file"
    else
        echo "Energy=$ener Lattice=\"$latt\" Config_type=$configuration-$i Weight=1.0 Properties=species:S:1:pos:R:3:forces:R:3" >> "$writ_dire/$writ_file"
    fi
    grep -A $syst_numb_atom "INDEX" temp.file | tail -n $syst_numb_atom | awk '{print $2,$3,$4,$5,$6,$7,$8}' >> $writ_dire/$writ_file
    echo -ne " Process: ${i}/${N_counts}\r"
    rm -f temp.file
done

echo
dos2unix "$writ_dire/$writ_file"
echo " All done."

