"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     nep_modifier.py
Category:   Model Modification Scripts
Purpose:    Interactively modify NEP models using calorine.
Usage:      python nep_modifier.py <input.nep>

Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-06-17
=============================================================================
"""

import os
import sys
from pathlib import Path

try:
    from calorine.nep import read_model, read_nepfile, write_nepfile
except ImportError:
    print(" Error: calorine package (>=3.4) is required.")
    print(" Install via pip:")
    print("   pip install 'calorine>=3.4'")
    print(" Or install the latest version from GitLab:")
    print("   pip install git+https://gitlab.com/materials-modeling/calorine.git@master")
    sys.exit(1)


def print_dependency_notice():
    print(" This function requires the calorine package.")
    print(" If you use this function, we recommend citing:")
    print(" Lindgren et al., J. Open Source Softw. 9, 6264 (2024).")
    print(" https://doi.org/10.21105/joss.06264")


# ─── helpers ──────────────────────────────────────────────────────────────────

def banner(title):
    print(f" +---------------------------------------------------+")
    print(f" | {title:<51}|")
    print(f" +---------------------------------------------------+")


def separator():
    print(" -----------------------------------------------------")


def ask(prompt, hint=None):
    if hint:
        print(f" {prompt} ({hint})")
    else:
        print(f" {prompt}")
    print(" ------------>>")
    return input(" ").strip()


def confirm(msg):
    print(f" {msg} (y/n)")
    print(" ------------>>")
    return input(" ").strip().lower() in ("y", "yes")


def safe_output_path(out_dir, base_name, ext):
    path = out_dir / f"{base_name}{ext}"
    if not path.exists():
        return path
    i = 1
    while True:
        path = out_dir / f"{base_name}_{i}{ext}"
        if not path.exists():
            return path
        i += 1


# ─── core functions ───────────────────────────────────────────────────────────

def show_model_info(model):
    separator()
    print(f" Version          : {model.version}")
    print(f" Model type       : {model.model_type}")
    print(f" Types            : {model.types}")
    print(f" Radial cutoff    : {model.radial_cutoff}")
    print(f" Angular cutoff   : {model.angular_cutoff}")
    print(f" n_max (rad/ang)  : {model.n_max_radial} / {model.n_max_angular}")
    print(f" basis_size       : {model.n_basis_radial} / {model.n_basis_angular}")
    print(f" l_max (3b/4b/5b) : {model.l_max_3b} / {model.l_max_4b} / {model.l_max_5b}")
    print(f" n_neuron         : {model.n_neuron}")
    print(f" n_parameters     : {model.n_parameters}")
    separator()


def op_augment(model):
    banner("augment: Expand Model")
    print(f" Current: neuron={model.n_neuron}, l_max_5b={model.l_max_5b}")
    separator()
    print(" Choose what to augment:")
    print(" 1) Increase n_neuron")
    print(" 2) Increase l_max_5b")
    print(" 3) Add charge_head")
    print(" 0) Cancel")
    print(" ------------>>")
    choice = input(" ").strip()

    kwargs = {}
    if choice == "1":
        val = ask(" Input new n_neuron:", str(model.n_neuron + 10))
        kwargs["n_neuron"] = int(val)
    elif choice == "2":
        val = ask(" Input new l_max_5b:", "1")
        kwargs["l_max_5b"] = int(val)
    elif choice == "3":
        kwargs["charge_head"] = True
    else:
        print(" Cancelled.")
        return model

    separator()
    for k, v in kwargs.items():
        print(f"   {k} = {v}")

    if not confirm(" Confirm augment?"):
        return model

    new_model = model.augment(**kwargs)
    print(f" Done: {model.n_parameters} -> {new_model.n_parameters} parameters")
    return new_model


def op_prune(model):
    banner("prune: Shrink Model")
    print(f" Current: neuron={model.n_neuron}, l_max_4b={model.l_max_4b}")
    separator()
    print(" Choose what to prune:")
    print(" 1) Reduce n_neuron")
    print(" 2) Reduce l_max_4b")
    print(" 0) Cancel")
    print(" ------------>>")
    choice = input(" ").strip()

    kwargs = {}
    if choice == "1":
        val = ask(" Input new n_neuron:", str(max(1, model.n_neuron - 10)))
        kwargs["n_neuron"] = int(val)
    elif choice == "2":
        val = ask(" Input new l_max_4b (0 to remove):", "0")
        kwargs["l_max_4b"] = int(val)
    else:
        print(" Cancelled.")
        return model

    separator()
    for k, v in kwargs.items():
        print(f"   {k} = {v}")

    if not confirm(" Confirm prune?"):
        return model

    new_model = model.prune(**kwargs)
    print(f" Done: {model.n_parameters} -> {new_model.n_parameters} parameters")
    return new_model


def op_add_species(model):
    banner("add_species: Add Elements")
    print(f" Current types: {model.types}")
    val = ask(" Input elements to add (space-separated):", "Bi Sn")
    species = val.split()
    if not species:
        print(" Cancelled.")
        return model

    separator()
    print(f" Will add: {species}")
    if not confirm(" Confirm?"):
        return model

    new_model = model.add_species(species)
    print(f" Done: {model.types} -> {new_model.types}")
    print(f" Params: {model.n_parameters} -> {new_model.n_parameters}")
    return new_model


def op_remove_species(model):
    banner("remove_species: Remove Elements")
    print(f" Current types: {model.types}")
    val = ask(" Input elements to remove (space-separated):")
    species = val.split()
    if not species:
        print(" Cancelled.")
        return model

    separator()
    print(f" Will remove: {species}")
    if not confirm(" Confirm?"):
        return model

    new_model = model.remove_species(species)
    print(f" Done: {model.types} -> {new_model.types}")
    print(f" Params: {model.n_parameters} -> {new_model.n_parameters}")
    return new_model


def op_keep_species(model):
    banner("keep_species: Keep Elements Only")
    print(f" Current types: {model.types}")
    val = ask(" Input elements to keep (space-separated):", "Li La Zr O")
    species = val.split()
    if not species:
        print(" Cancelled.")
        return model

    separator()
    print(f" Will keep: {species}")
    if not confirm(" Confirm?"):
        return model

    new_model = model.keep_species(species)
    print(f" Done: {model.types} -> {new_model.types}")
    print(f" Params: {model.n_parameters} -> {new_model.n_parameters}")
    return new_model


def op_export(model, original_path):
    banner("Export Model Files")
    out_dir = ask(" Output directory:", ".")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = Path(original_path).stem
    nep_out = safe_output_path(out_dir, base + "_modified", ".txt")
    restart_out = out_dir / (nep_out.stem + ".restart")

    model.write(str(nep_out))
    print(f" Written: {nep_out}")

    try:
        model.write_restart(str(restart_out))
        print(f" Written: {restart_out}")
    except Exception:
        print(" Skipped: nep.restart (no restart data loaded)")

    try:
        if hasattr(model, "training_parameters"):
            nep_in = out_dir / (nep_out.stem + ".in")
            params = model.training_parameters
            if Path("nep.in").exists():
                old = read_nepfile("nep.in")
                old.update(params)
                params = old
            write_nepfile(params, str(out_dir))
            default_in = out_dir / "nep.in"
            if default_in.exists() and default_in != nep_in:
                default_in.rename(nep_in)
            print(f" Written: {nep_in}")
    except Exception:
        pass

    separator()


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    print("")
    print(" +----------------------------------------------------+")
    print(" |         NEP Model Modifier (calorine)              |")
    print(" +----------------------------------------------------+")
    print(" |  Powered by calorine package. If you use this tool,|")
    print(" |  please cite the following papers:                 |")
    print(" |                                                    |")
    print(" |   Calorine: A Python package for constructing and  |")
    print(" |      sampling neuroevolution potential models.     |")
    print(" |   Journal of Open Source Software, 2024, 9, 6264.  |")
    print(" |      https://doi.org/10.21105/joss.06264           |")
    print(" +----------------------------------------------------+")
    print(" [!] WARNING: This feature is currently in TESTING.")
    print(" [!] The docs were generated by AI and may contain errors.")
    print(" [!] Report issues or join the QQ group: 825696376")
    print(" +----------------------------------------------------+")
    print("")

    val = ask(" Input <nep.txt> <nep.restart> (Enter to use defaults):")
    parts = val.split()
    nep_path = parts[0] if parts else "nep.txt"
    restart_path = parts[1] if len(parts) > 1 else "nep.restart"

    if not Path(nep_path).exists():
        print(f" Error: file not found: {nep_path}")
        return
    if not Path(restart_path).exists():
        print(f" Warning: {restart_path} not found, skipping restart.")
        restart_path = None

    separator()
    print(" Loading model...")
    if restart_path:
        model = read_model(nep_path, restart_file=restart_path)
    else:
        model = read_model(nep_path)
    print(" Model loaded successfully!")

    show_model_info(model)

    while True:
        print("")
        print(" Choose an operation:")
        print(" 1) augment          (expand neurons/descriptors/charge)")
        print(" 2) prune            (shrink neurons/descriptors)")
        print(" 3) add_species      (add element)")
        print(" 4) remove_species   (remove element)")
        print(" 5) keep_species     (keep selected elements only)")
        print(" 6) Show model info")
        print(" 7) Export model")
        print(" 0) Exit")
        print(" ------------>>")
        choice = input(" ").strip()

        if choice == "0":
            print(" Bye!")
            break
        elif choice == "1":
            model = op_augment(model)
        elif choice == "2":
            model = op_prune(model)
        elif choice == "3":
            model = op_add_species(model)
        elif choice == "4":
            model = op_remove_species(model)
        elif choice == "5":
            model = op_keep_species(model)
        elif choice == "6":
            show_model_info(model)
        elif choice == "7":
            op_export(model, nep_path)
        else:
            print(" Invalid option, try again.")


if __name__ == "__main__":
    print_dependency_notice()
    main()
