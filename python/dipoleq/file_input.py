"""Load the data from the cli tool *.in file
Args:
    file_path (str): The path to the *.in file
Returns:
    dict: The data from the *.in file
"""

from pathlib import Path
from typing import Any

from .input import MachineIn, PlasmaIn


def read_dotin(file_path: str | Path) -> dict[str, Any]:
    """Load the data from the cli tool *.in file
    Args:
        file_path (str): The path to the *.in file
    Returns:
        dict: The data from the *.in file
    """
    with Path(file_path).open(encoding="utf-8") as file:
        lines = file.readlines()

    path = Path(file_path)

    top_dict: dict[str, Any] = {
        "Iname": str(path.name),
        "PsiGrid": {},
        "Plasma": {},
        "Limiters": [],
        "Coils": [],
        "Shells": [],
        "Separatrices": [],
        "Measures": [],
    }
    cur_dict: dict[str, Any] = {}
    cur_coil = {}
    cur_shell = {}
    key = ""
    for liner in lines:
        line = liner.strip()
        if line.startswith("//") or line == "":
            continue
        if line.startswith("K_"):
            key = line[2:]
            cur_dict = {}
            match key:
                case "Coil":
                    cur_dict["SubCoils"] = []
                    cur_coil = cur_dict
                    top_dict["Coils"].append(cur_dict)
                case "SubCoil":
                    cur_coil["SubCoils"].append(cur_dict)
                case "Shell":
                    cur_dict["SubShells"] = []
                    cur_shell = cur_dict
                    top_dict["Shells"].append(cur_dict)
                case "SubShell":
                    cur_shell["SubShells"].append(cur_dict)
                case "Separatrix":
                    top_dict["Separatrices"].append(cur_dict)
                case "Measure":
                    top_dict["Measures"].append(cur_dict)
                case "Limiter":
                    top_dict["Limiters"].append(cur_dict)
                case "CodeControl":
                    cur_dict = top_dict
                case _:
                    cur_dict = top_dict[key]
        else:
            sub_key, value = line.split("=")
            sub_key = sub_key.strip()
            value = value.strip()
            cur_dict[sub_key] = value

    return top_dict


def input_from_dotin(filepath: str | Path) -> MachineIn:
    """Read the input file and return a dictionary of the input data
    Args:
        filepath (str): The path to the input file
    Returns:
        MachineIn: the verified input data
    """
    info = read_dotin(filepath)

    # demote plasma model data to sub dictionary...  matches new machine input
    pl = info["Plasma"]
    model = {}
    keys_to_move = [key for key in pl if key not in PlasmaIn.model_fields]
    for key in keys_to_move:
        model[key] = pl.pop(key)

    if "ModelType" in model:
        model["Type"] = model.pop("ModelType")

    pl["Model"] = model

    return MachineIn(**info)


def input_from_yaml(filepath: str | Path) -> MachineIn:
    """Read the input file and return a dictionary of the input data
    Args:
        filepath (str): The path to the input file
    Returns:
        MachineIn: the verified input data
    """
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("pyyaml is required to read yaml files") from exc

    with Path(filepath).open(encoding="utf-8") as stream:
        info = yaml.safe_load(stream)

    path = Path(filepath)

    info["Iname"] = str(path.name)
    return MachineIn(**info)
