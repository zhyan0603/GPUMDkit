"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     nep_modifier.py
Category:   Model Modification Scripts
Purpose:    Safely inspect and interactively modify NEP4 models with calorine.
Usage:      gpumdkit.sh -nep_modifier [nep.txt] [nep.restart|-] [nep.in|-]
            python nep_modifier.py [nep.txt] [nep.restart|-] [nep.in|-]
Arguments:
  nep.txt      NEP model file (default: ./nep.txt)
  nep.restart  Matching SNES restart file; use '-' to skip
  nep.in       Training input whose non-architecture settings are preserved;
               use '-' to generate architecture fields only
Output:
  *_modified.txt          Modified model
  *_modified.restart      Matching restart, when restart data are available
  *_modified.in           Updated training input
  *_modified.changes.txt  Modification and provenance summary
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-07-12
=============================================================================
"""

import inspect
import shutil
import sys
import tempfile
import textwrap
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


MIN_CALORINE_VERSION = (3, 4)
BOX_WIDTH = 79


def print_usage():
    """Print command-line usage without importing optional dependencies."""
    print(" Usage: gpumdkit.sh -nep_modifier [nep.txt] [nep.restart|-] [nep.in|-]")
    print("    or: python nep_modifier.py [nep.txt] [nep.restart|-] [nep.in|-]")
    print("")
    print(" Arguments:")
    print("   nep.txt      NEP model file (default: ./nep.txt)")
    print("   nep.restart  Matching restart file; use '-' to load no restart")
    print("   nep.in       Training input to update; use '-' for architecture only")
    print("")
    print(" The command starts an interactive model editor and never overwrites inputs.")
    print("")


CLI_ARGS = sys.argv[1:]
if CLI_ARGS and CLI_ARGS[0] in ("-h", "--help"):
    print_usage()
    sys.exit(0)
if len(CLI_ARGS) > 3:
    print_usage()
    sys.exit(1)


try:
    from calorine.nep import read_model, read_nepfile, write_nepfile
except ImportError:
    print(" Error: calorine package (>=3.4) is required.")
    print(" Install it with: pip install 'calorine>=3.4'")
    sys.exit(1)


class InputClosed(Exception):
    """Raised when interactive input is closed or interrupted."""


def version_tuple(value):
    """Return the leading numeric components of a package version."""
    parts = []
    for part in value.split("."):
        digits = "".join(character for character in part if character.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


try:
    CALORINE_VERSION = version("calorine")
except PackageNotFoundError:
    CALORINE_VERSION = "unknown"

if CALORINE_VERSION != "unknown" and version_tuple(CALORINE_VERSION) < MIN_CALORINE_VERSION:
    print(f" Error: calorine >=3.4 is required; found {CALORINE_VERSION}.")
    sys.exit(1)


def print_box(title, lines, width=BOX_WIDTH):
    """Print an ASCII information box with aligned borders."""
    print(" +" + "-" * width + "+")
    print(f" |{title:^{width}}|")
    print(" +" + "-" * width + "+")
    for line in lines:
        wrapped_lines = textwrap.wrap(
            line,
            width=width - 2,
            subsequent_indent="  ",
            break_long_words=False,
            break_on_hyphens=False,
        ) or [""]
        for wrapped_line in wrapped_lines:
            print(f" | {wrapped_line:<{width - 1}}|")
    print(" +" + "-" * width + "+")


def read_input(prompt, default=None, required=False):
    """Read one line safely, applying an explicit default when requested."""
    while True:
        if default is None:
            print(f" {prompt}")
        else:
            print(f" {prompt} [default: {default}]")
        print(" ------------>>")
        try:
            value = input(" ").strip()
        except (EOFError, KeyboardInterrupt) as error:
            raise InputClosed from error
        if value:
            return value
        if default is not None:
            return str(default)
        if not required:
            return ""
        print(" Error: a value is required.")


def ask_yes_no(prompt, default=None):
    """Ask a validated yes/no question."""
    if default is True:
        suffix = "Y/n"
    elif default is False:
        suffix = "y/N"
    else:
        suffix = "y/n"
    while True:
        value = read_input(f"{prompt} ({suffix})")
        if not value and default is not None:
            return default
        if value.lower() in ("y", "yes"):
            return True
        if value.lower() in ("n", "no"):
            return False
        print(" Error: please input 'y' or 'n'.")


def ask_int(prompt, default=None, minimum=None, maximum=None):
    """Read and validate an integer."""
    while True:
        value = read_input(prompt, default=default, required=default is None)
        try:
            result = int(value)
        except ValueError:
            print(" Error: please input an integer.")
            continue
        if minimum is not None and result < minimum:
            print(f" Error: the value must be at least {minimum}.")
            continue
        if maximum is not None and result > maximum:
            print(f" Error: the value must not exceed {maximum}.")
            continue
        return result


def ask_float(prompt, default=None, minimum=None, maximum=None):
    """Read and validate a floating-point value."""
    while True:
        value = read_input(prompt, default=default, required=default is None)
        try:
            result = float(value)
        except ValueError:
            print(" Error: please input a number.")
            continue
        if minimum is not None and result < minimum:
            print(f" Error: the value must be at least {minimum}.")
            continue
        if maximum is not None and result > maximum:
            print(f" Error: the value must not exceed {maximum}.")
            continue
        return result


def ask_options(prompt, valid_choices):
    """Read one or more unique menu choices."""
    while True:
        raw = read_input(prompt, required=True)
        choices = raw.split()
        if choices and all(choice in valid_choices for choice in choices):
            return list(dict.fromkeys(choices))
        print(f" Error: valid choices are {' '.join(valid_choices)}.")


def has_restart(model):
    """Return whether restart parameters are loaded."""
    return getattr(model, "restart_parameters", None) is not None


def require_restart(model, operation):
    """Report a restart prerequisite and return whether it is satisfied."""
    if has_restart(model):
        return True
    print(f" Error: {operation} requires a matching nep.restart file.")
    print(" Reload the tool with both nep.txt and nep.restart.")
    return False


def supports_parameter(method, parameter):
    """Check a calorine method signature for version-specific parameters."""
    return parameter in inspect.signature(method).parameters


def ask_snes_settings(operation, include_sigma_new=False):
    """Ask whether to use documented calorine SNES initialization defaults."""
    defaults = {"sigma_factor": 0.1, "sigma_floor": 1.0e-6}
    if include_sigma_new:
        defaults["sigma_new"] = 0.1 if operation == "add_species" else 0.01
    values = ", ".join(f"{key}={value:g}" for key, value in defaults.items())
    print(f" Calorine documented defaults: {values}")
    if ask_yes_no("Use these SNES initialization defaults?", default=True):
        return {}
    settings = {}
    if include_sigma_new:
        settings["sigma_new"] = ask_float(
            "Input sigma_new:", default=defaults["sigma_new"], minimum=0.0
        )
    settings["sigma_factor"] = ask_float(
        "Input sigma_factor:", default=defaults["sigma_factor"], minimum=0.0
    )
    settings["sigma_floor"] = ask_float(
        "Input sigma_floor:", default=defaults["sigma_floor"], minimum=0.0
    )
    return settings


def model_snapshot(model):
    """Return compact model fields used for change summaries."""
    return {
        "types": tuple(model.types),
        "n_neuron": model.n_neuron,
        "l_max_4b": model.l_max_4b,
        "l_max_5b": model.l_max_5b,
        "has_q_112": int(getattr(model, "has_q_112", 0)),
        "has_q_123": int(getattr(model, "has_q_123", 0)),
        "has_q_233": int(getattr(model, "has_q_233", 0)),
        "has_q_134": int(getattr(model, "has_q_134", 0)),
        "model_type": model.model_type,
        "n_parameters": model.n_parameters,
    }


def format_change(before, after):
    """Format changed model fields for the terminal and provenance log."""
    lines = []
    for key in before:
        if before[key] != after[key]:
            lines.append(f"{key}: {before[key]} -> {after[key]}")
    return lines


def show_model_info(model):
    """Display model architecture, restart, cutoff, and parameter information."""
    q_terms = " ".join(
        f"{name}={int(getattr(model, f'has_q_{name}', 0))}"
        for name in ("112", "123", "233", "134")
    )
    lines = [
        f"Version: {model.version:<12} Model type: {model.model_type}",
        f"Species ({len(model.types)}): {' '.join(model.types)}",
        f"Cutoff radial/angular: {model.radial_cutoff} / {model.angular_cutoff} Angstrom",
        f"n_max radial/angular: {model.n_max_radial} / {model.n_max_angular}",
        f"basis radial/angular: {model.n_basis_radial} / {model.n_basis_angular}",
        f"l_max 3b/4b/5b: {model.l_max_3b} / {model.l_max_4b} / {model.l_max_5b}",
        f"Higher-body switches: {q_terms}",
        f"Descriptors radial/angular: {model.n_descriptor_radial} / {model.n_descriptor_angular}",
        f"Neurons: {model.n_neuron:<12} Parameters: {model.n_parameters}",
        f"ANN parameters: {model.n_ann_parameters:<7} Descriptor parameters: {model.n_descriptor_parameters}",
        f"Restart data: {'loaded' if has_restart(model) else 'not loaded'}",
        f"ZBL: {getattr(model, 'zbl', None)}",
    ]
    if hasattr(model, "charge_mode"):
        lines.append(f"Charge mode: {model.charge_mode}")
    print_box("CURRENT NEP MODEL", lines)


def print_main_menu():
    """Print the two-column model modification menu."""
    rows = [
        ("1) Expand model capacity", "4) Remove chemical species"),
        ("2) Reduce model capacity", "5) Keep selected species"),
        ("3) Add chemical species", "6) Inspect current model"),
        ("7) Review pending changes", "8) Export model package"),
        ("0) Exit", ""),
    ]
    lines = [f"{left:<37}{right}".rstrip() for left, right in rows]
    print_box("NEP MODEL MODIFIER", lines)


def print_change_result(before, after):
    """Print fields changed by a successful operation."""
    changes = format_change(before, after)
    print(" Modification completed.")
    for change in changes:
        print(f"   {change}")
    return changes


def op_augment(model):
    """Collect and atomically apply one or more model expansions."""
    if not require_restart(model, "Expanding a model"):
        return model, None
    lines = [
        f"1) Neurons ({model.n_neuron})                  5) q_123 ({model.has_q_123})",
        f"2) l_max_4b ({model.l_max_4b})                 6) q_233 ({model.has_q_233})",
        f"3) l_max_5b ({model.l_max_5b})                 7) q_134 ({model.has_q_134})",
        f"4) q_112 ({model.has_q_112})                   8) Charge output head",
        "0) Cancel",
    ]
    print_box("EXPAND MODEL CAPACITY", lines)
    choices = ask_options("Select one or more changes, for example: 1 3", list("012345678"))
    if "0" in choices:
        print(" Operation canceled.")
        return model, None

    kwargs = {}
    if "1" in choices:
        kwargs["n_neuron"] = ask_int(
            f"Input target neuron count (current: {model.n_neuron}):",
            minimum=model.n_neuron + 1,
            maximum=120,
        )
    if "2" in choices:
        kwargs["l_max_4b"] = ask_int(
            f"Input target l_max_4b (current: {model.l_max_4b}):",
            minimum=model.l_max_4b + 1,
        )
    if "3" in choices:
        kwargs["l_max_5b"] = ask_int(
            f"Input target l_max_5b (current: {model.l_max_5b}):",
            minimum=model.l_max_5b + 1,
        )
    for choice, name in (("4", "112"), ("5", "123"), ("6", "233"), ("7", "134")):
        if choice in choices:
            if int(getattr(model, f"has_q_{name}", 0)) == 1:
                print(f" Notice: q_{name} is already enabled; no change will be made.")
            else:
                kwargs[f"has_q_{name}"] = True
    if "8" in choices:
        if model.model_type != "potential":
            print(f" Notice: a charge head cannot be added to model type '{model.model_type}'.")
        else:
            kwargs["charge_head"] = True
            if supports_parameter(model.augment, "charge_mode"):
                print(" Charge mode 1: real-space and reciprocal-space qNEP.")
                print(" Charge mode 2: reciprocal-space-only qNEP.")
                kwargs["charge_mode"] = ask_int("Input charge mode (1 or 2):", minimum=1, maximum=2)
    if not kwargs:
        print(" Operation canceled: no applicable changes were selected.")
        return model, None

    kwargs.update(ask_snes_settings("augment", include_sigma_new=True))
    print(" Requested changes:")
    for key, value in kwargs.items():
        print(f"   {key} = {value}")
    if not ask_yes_no("Apply these changes?", default=False):
        print(" Operation canceled.")
        return model, None

    before = model_snapshot(model)
    try:
        new_model = model.augment(**kwargs)
    except (TypeError, ValueError, OSError) as error:
        print(f" Error: model expansion failed: {error}.")
        return model, None
    changes = print_change_result(before, model_snapshot(new_model))
    return new_model, {"operation": "augment", "arguments": kwargs, "changes": changes, "retraining": True}


def op_prune(model):
    """Collect and atomically apply one or more model reductions."""
    if not require_restart(model, "Reducing a model"):
        return model, None
    lines = [
        f"1) Neurons ({model.n_neuron})                  5) q_123 ({model.has_q_123})",
        f"2) l_max_4b ({model.l_max_4b})                 6) q_233 ({model.has_q_233})",
        f"3) l_max_5b ({model.l_max_5b})                 7) q_134 ({model.has_q_134})",
        f"4) q_112 ({model.has_q_112})                   8) Charge output head",
        "0) Cancel",
    ]
    print_box("REDUCE MODEL CAPACITY", lines)
    choices = ask_options("Select one or more changes, for example: 1 2", list("012345678"))
    if "0" in choices:
        print(" Operation canceled.")
        return model, None

    kwargs = {}
    if "1" in choices:
        if model.n_neuron <= 1:
            print(" Notice: the neuron count cannot be reduced below one.")
        else:
            kwargs["n_neuron"] = ask_int(
                f"Input target neuron count (current: {model.n_neuron}):",
                minimum=1,
                maximum=model.n_neuron - 1,
            )
    if "2" in choices:
        if model.l_max_4b <= 0:
            print(" Notice: l_max_4b is already zero.")
        else:
            kwargs["l_max_4b"] = ask_int(
                f"Input target l_max_4b (current: {model.l_max_4b}):",
                minimum=0,
                maximum=model.l_max_4b - 1,
            )
    if "3" in choices:
        if model.l_max_5b <= 0:
            print(" Notice: l_max_5b is already zero.")
        else:
            kwargs["l_max_5b"] = ask_int(
                f"Input target l_max_5b (current: {model.l_max_5b}):",
                minimum=0,
                maximum=model.l_max_5b - 1,
            )
    for choice, name in (("4", "112"), ("5", "123"), ("6", "233"), ("7", "134")):
        if choice in choices:
            if int(getattr(model, f"has_q_{name}", 0)) == 0:
                print(f" Notice: q_{name} is already disabled; no change will be made.")
            else:
                kwargs[f"has_q_{name}"] = False
    if "8" in choices:
        if model.model_type != "potential_with_charges":
            print(" Notice: the model does not contain a charge output head.")
        else:
            kwargs["charge_head"] = True
    if not kwargs:
        print(" Operation canceled: no applicable changes were selected.")
        return model, None

    kwargs.update(ask_snes_settings("prune"))
    print(" Requested changes:")
    for key, value in kwargs.items():
        print(f"   {key} = {value}")
    if not ask_yes_no("Apply these changes?", default=False):
        print(" Operation canceled.")
        return model, None

    before = model_snapshot(model)
    try:
        new_model = model.prune(**kwargs)
    except (TypeError, ValueError, OSError) as error:
        print(f" Error: model reduction failed: {error}.")
        return model, None
    changes = print_change_result(before, model_snapshot(new_model))
    return new_model, {"operation": "prune", "arguments": kwargs, "changes": changes, "retraining": True}


def read_species(prompt):
    """Read a whitespace-separated list of unique chemical symbols."""
    while True:
        species = read_input(prompt, required=True).split()
        if len(species) != len(set(species)):
            print(" Error: duplicate species are not allowed.")
            continue
        return species


def is_typewise_cutoff(model):
    """Return whether the model stores per-species cutoff values."""
    return isinstance(model.radial_cutoff, (list, tuple))


def op_add_species(model):
    """Add species with explicit reproducibility and cutoff inputs."""
    if not require_restart(model, "Adding species"):
        return model, None
    if model.version != 4:
        print(f" Error: add_species supports NEP4 models; found version {model.version}.")
        return model, None

    print_box("ADD CHEMICAL SPECIES", [f"Current species: {' '.join(model.types)}"])
    species = read_species("Input species to add, separated by spaces:")
    existing = sorted(set(species).intersection(model.types))
    if existing:
        print(f" Error: species already present in the model: {' '.join(existing)}.")
        return model, None

    kwargs = {"species": species}
    if is_typewise_cutoff(model):
        radial_cutoffs = []
        angular_cutoffs = []
        for symbol in species:
            radial = ask_float(f"Input radial cutoff for {symbol} (Angstrom):", minimum=0.0)
            angular = ask_float(
                f"Input angular cutoff for {symbol} (Angstrom):", minimum=0.0, maximum=radial
            )
            radial_cutoffs.append(radial)
            angular_cutoffs.append(angular)
        kwargs["radial_cutoff"] = radial_cutoffs
        kwargs["angular_cutoff"] = angular_cutoffs

    kwargs["seed"] = ask_int("Input random seed for reproducible initialization:", minimum=0)
    kwargs.update(ask_snes_settings("add_species", include_sigma_new=True))
    print(" Requested changes:")
    for key, value in kwargs.items():
        print(f"   {key} = {value}")
    if not ask_yes_no("Add these species?", default=False):
        print(" Operation canceled.")
        return model, None

    before = model_snapshot(model)
    try:
        new_model = model.add_species(**kwargs)
    except (TypeError, ValueError, OSError) as error:
        print(f" Error: adding species failed: {error}.")
        return model, None
    changes = print_change_result(before, model_snapshot(new_model))
    return new_model, {"operation": "add_species", "arguments": kwargs, "changes": changes, "retraining": True}


def op_remove_species(model):
    """Remove selected species after validation and confirmation."""
    print_box("REMOVE CHEMICAL SPECIES", [f"Current species: {' '.join(model.types)}"])
    species = read_species("Input species to remove, separated by spaces:")
    missing = sorted(set(species).difference(model.types))
    if missing:
        print(f" Error: species not present in the model: {' '.join(missing)}.")
        return model, None
    if len(species) >= len(model.types):
        print(" Error: at least one species must remain in the model.")
        return model, None

    kwargs = {}
    if has_restart(model):
        kwargs.update(ask_snes_settings("remove_species"))
    print(f" Species to remove: {' '.join(species)}")
    if not ask_yes_no("Remove these species?", default=False):
        print(" Operation canceled.")
        return model, None

    before = model_snapshot(model)
    try:
        new_model = model.remove_species(species, **kwargs)
    except (TypeError, ValueError, OSError) as error:
        print(f" Error: removing species failed: {error}.")
        return model, None
    changes = print_change_result(before, model_snapshot(new_model))
    arguments = {"species": species, **kwargs}
    return new_model, {"operation": "remove_species", "arguments": arguments, "changes": changes, "retraining": False}


def op_keep_species(model):
    """Retain selected species after validation and confirmation."""
    print_box("KEEP SELECTED SPECIES", [f"Current species: {' '.join(model.types)}"])
    species = read_species("Input species to keep, separated by spaces:")
    missing = sorted(set(species).difference(model.types))
    if missing:
        print(f" Error: species not present in the model: {' '.join(missing)}.")
        return model, None

    kwargs = {}
    if has_restart(model):
        kwargs.update(ask_snes_settings("keep_species"))
    print(f" Species to keep: {' '.join(species)}")
    if not ask_yes_no("Keep only these species?", default=False):
        print(" Operation canceled.")
        return model, None

    before = model_snapshot(model)
    try:
        new_model = model.keep_species(species, **kwargs)
    except (TypeError, ValueError, OSError) as error:
        print(f" Error: keeping species failed: {error}.")
        return model, None
    changes = print_change_result(before, model_snapshot(new_model))
    arguments = {"species": species, **kwargs}
    return new_model, {"operation": "keep_species", "arguments": arguments, "changes": changes, "retraining": False}


def show_history(history, exported_count):
    """Display pending and previously exported operations."""
    if not history:
        print(" No model modifications have been made.")
        return
    print_box("MODEL MODIFICATION HISTORY", [])
    for index, entry in enumerate(history, 1):
        state = "exported" if index <= exported_count else "pending"
        print(f" {index}) {entry['operation']} ({state})")
        for key, value in entry["arguments"].items():
            print(f"      {key} = {value}")
        for change in entry["changes"]:
            print(f"      {change}")
    if any(entry["retraining"] for entry in history[exported_count:]):
        print(" Pending changes include operations that require further training.")


def safe_bundle_paths(out_dir, source_stem, include_restart):
    """Choose one suffix for a collision-free model output package."""
    base = source_stem + "_modified"
    index = 0
    while True:
        stem = base if index == 0 else f"{base}_{index}"
        paths = {
            "model": out_dir / f"{stem}.txt",
            "input": out_dir / f"{stem}.in",
            "summary": out_dir / f"{stem}.changes.txt",
        }
        if include_restart:
            paths["restart"] = out_dir / f"{stem}.restart"
        if not any(path.exists() for path in paths.values()):
            return stem, paths
        index += 1


def build_change_summary(model_path, restart_path, nep_in_path, model, history):
    """Build a reproducibility summary for the exported model package."""
    lines = [
        "GPUMDkit NEP model modification summary",
        f"Created UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Calorine version: {CALORINE_VERSION}",
        f"Source model: {model_path.resolve()}",
        f"Source restart: {restart_path.resolve() if restart_path else 'not loaded'}",
        f"Source nep.in: {nep_in_path.resolve() if nep_in_path else 'not provided'}",
        f"Final species: {' '.join(model.types)}",
        f"Final neuron count: {model.n_neuron}",
        f"Final parameter count: {model.n_parameters}",
        "",
        "Operations:",
    ]
    if not history:
        lines.append("  None (package exported without model modification)")
    for index, entry in enumerate(history, 1):
        lines.append(f"  {index}. {entry['operation']}")
        for key, value in entry["arguments"].items():
            lines.append(f"     {key} = {value}")
        for change in entry["changes"]:
            lines.append(f"     {change}")
        lines.append(f"     further training required: {'yes' if entry['retraining'] else 'optional'}")
    return "\n".join(lines) + "\n"


def export_model(model, model_path, restart_path, nep_in_path, history):
    """Export a collision-free and internally consistent model package."""
    default_dir = model_path.parent
    out_dir = Path(read_input("Output directory:", default=str(default_dir))).expanduser()
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        print(f" Error: failed to create output directory: {error}.")
        return False

    include_restart = has_restart(model)
    stem, final_paths = safe_bundle_paths(out_dir, model_path.stem, include_restart)
    architecture_params = dict(model.training_parameters)
    if nep_in_path is not None:
        try:
            old_params = read_nepfile(str(nep_in_path))
        except (OSError, ValueError, TypeError) as error:
            print(f" Error: failed to read source nep.in: {error}.")
            return False
        # Current GPUMD uses a single outer cutoff in nep.in, while calorine's
        # model representation stores the inner and outer ZBL cutoffs. Preserve
        # the source nep.in spelling instead of replacing it with a two-value
        # model representation that calorine.read_nepfile cannot read back.
        if "zbl" in old_params:
            architecture_params.pop("zbl", None)
        old_params.update(architecture_params)
        params = old_params
    else:
        print(" Notice: no source nep.in was supplied; exporting architecture fields only.")
        params = architecture_params
        if isinstance(params.get("zbl"), (list, tuple)):
            params["zbl"] = params["zbl"][-1]
        zbl_factor = getattr(model, "zbl_typewise_cutoff_factor", None)
        if zbl_factor is not None:
            params["use_typewise_cutoff_zbl"] = zbl_factor

    if "type_weight" in params and len(params["type_weight"]) != len(model.types):
        print(
            " Error: source type_weight count does not match the modified species count; "
            "update the weights explicitly before export."
        )
        return False

    print(" Output package:")
    for path in final_paths.values():
        print(f"   {path}")
    if not ask_yes_no("Write this model package?", default=False):
        print(" Operation canceled.")
        return False

    try:
        with tempfile.TemporaryDirectory(prefix=".nep_modifier_", dir=out_dir) as temp_name:
            temp_dir = Path(temp_name)
            temp_model = temp_dir / final_paths["model"].name
            model.write(str(temp_model))
            temp_paths = {"model": temp_model}

            if include_restart:
                temp_restart = temp_dir / final_paths["restart"].name
                model.write_restart(str(temp_restart))
                temp_paths["restart"] = temp_restart

            write_nepfile(params, str(temp_dir))
            generated_input = temp_dir / "nep.in"
            temp_input = temp_dir / final_paths["input"].name
            generated_input.rename(temp_input)
            temp_paths["input"] = temp_input

            # Verify that calorine can parse the file it just wrote and that
            # the architecture fields match the modified model.
            parsed_input = read_nepfile(str(temp_input))
            for key in ("version", "type", "cutoff", "n_max", "basis_size", "l_max", "neuron"):
                if key not in parsed_input:
                    raise ValueError(f"generated nep.in is missing architecture field '{key}'")
                if parsed_input[key] != params[key]:
                    raise ValueError(f"generated nep.in field '{key}' failed round-trip validation")

            temp_summary = temp_dir / final_paths["summary"].name
            temp_summary.write_text(
                build_change_summary(model_path, restart_path, nep_in_path, model, history),
                encoding="utf-8",
            )
            temp_paths["summary"] = temp_summary

            missing = [name for name, path in temp_paths.items() if not path.is_file()]
            if missing:
                raise OSError(f"temporary outputs are missing: {', '.join(missing)}")
            for name, temp_path in temp_paths.items():
                shutil.move(str(temp_path), str(final_paths[name]))
    except (OSError, TypeError, ValueError) as error:
        print(f" Error: model package export failed: {error}.")
        return False

    print(f" Model package exported with stem: {stem}")
    if any(entry["retraining"] for entry in history):
        print(" Further training is required before using the modified model in production.")
    print(" Validate the generated nep.in with the target NEP executable before training.")
    return True


def resolve_paths(cli_args):
    """Resolve model, restart, and training input paths from CLI or prompts."""
    if cli_args:
        model_path = Path(cli_args[0]).expanduser()
    else:
        model_path = Path(read_input("NEP model file:", default="nep.txt")).expanduser()

    if len(cli_args) >= 2:
        restart_value = cli_args[1]
    else:
        default_restart = model_path.with_name("nep.restart")
        restart_value = read_input("Restart file ('-' to skip):", default=str(default_restart))
    restart_path = None if restart_value == "-" else Path(restart_value).expanduser()

    if len(cli_args) >= 3:
        nep_in_value = cli_args[2]
    else:
        default_nep_in = model_path.with_name("nep.in")
        nep_in_value = read_input("Training input ('-' to skip):", default=str(default_nep_in))
    nep_in_path = None if nep_in_value == "-" else Path(nep_in_value).expanduser()
    return model_path, restart_path, nep_in_path


def load_model(model_path, restart_path, nep_in_path):
    """Validate input paths and load a calorine model."""
    if not model_path.is_file():
        print(f" Error: model file '{model_path}' does not exist.")
        return None, None, None
    if restart_path is not None and not restart_path.is_file():
        print(f" Notice: restart file '{restart_path}' was not found; loading the model without it.")
        restart_path = None
    if nep_in_path is not None and not nep_in_path.is_file():
        print(f" Notice: training input '{nep_in_path}' was not found; architecture-only export will be used.")
        nep_in_path = None

    print(" Loading model files...")
    try:
        if restart_path is None:
            model = read_model(str(model_path))
        else:
            model = read_model(str(model_path), restart_file=str(restart_path))
    except (OSError, TypeError, ValueError, IndexError) as error:
        print(f" Error: failed to load the model package: {error}.")
        return None, None, None
    print(" Model loaded successfully.")
    return model, restart_path, nep_in_path


def print_header():
    """Print program, dependency, and complete calorine citation information."""
    lines = [
        f"Calorine version: {CALORINE_VERSION}",
        "This functionality is built on calorine. Please cite:",
        "E. Lindgren, J. M. Rahm, E. Fransson, F. Eriksson, N. Österbacka, "
        "Z. Fan, and P. Erhart, \"calorine: A Python package for constructing "
        "and sampling neuroevolution potential models,\" Journal of Open "
        "Source Software 9(95), 6264 (2024).",
        "DOI: https://doi.org/10.21105/joss.06264",
        "This feature is under testing; validate outputs before further training.",
    ]
    print_box("GPUMDkit NEP MODEL MODIFIER", lines)


def main():
    """Run the interactive model modification session."""
    print_header()
    model_path, restart_path, nep_in_path = resolve_paths(CLI_ARGS)
    model, restart_path, nep_in_path = load_model(model_path, restart_path, nep_in_path)
    if model is None:
        return 1

    show_model_info(model)
    history = []
    exported_count = 0

    while True:
        print_main_menu()
        choice = read_input("Input the function number:", required=True)
        if choice == "0":
            if len(history) > exported_count:
                if not ask_yes_no("Exit without exporting pending changes?", default=False):
                    continue
            print(" Bye!")
            return 0
        if choice == "1":
            model, entry = op_augment(model)
        elif choice == "2":
            model, entry = op_prune(model)
        elif choice == "3":
            model, entry = op_add_species(model)
        elif choice == "4":
            model, entry = op_remove_species(model)
        elif choice == "5":
            model, entry = op_keep_species(model)
        elif choice == "6":
            show_model_info(model)
            entry = None
        elif choice == "7":
            show_history(history, exported_count)
            entry = None
        elif choice == "8":
            if export_model(model, model_path, restart_path, nep_in_path, history):
                exported_count = len(history)
            entry = None
        else:
            print(" Error: invalid option; please try again.")
            entry = None
        if entry is not None:
            history.append(entry)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except InputClosed:
        print("")
        print(" Input closed. Exiting without additional changes.")
        sys.exit(1)
