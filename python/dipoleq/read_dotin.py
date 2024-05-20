""" Load the data from the cli tool *.in file
    Args:
        file_path (str): The path to the *.in file
    Returns:
        dict: The data from the *.in file
"""
from typing import Any


def read_dotin(file_path: str) -> dict[str, Any]:
    """Load the data from the cli tool *.in file
        Args:
            file_path (str): The path to the *.in file
        Returns:
            dict: The data from the *.in file
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    top_dict: dict[str, Any] = dict(PsiGrid={}, Plasma={},
                                    Limiters=[], Coils=[], Shells=[],
                                    Separatricies=[], Measures=[])
    cur_dict: dict[str, Any] = {}
    cur_coil = {}
    cur_shell = {}
    key = ''
    for line in lines:
        line = line.strip()
        if line.startswith('//') or line == '':
            continue
        elif line.startswith('K_'):
            key = line[2:]
            cur_dict = {}
            match key:
                case 'Coil':
                    cur_dict['SubCoils'] = []
                    cur_coil = cur_dict
                    top_dict['Coils'].append(cur_dict)
                case 'SubCoil':
                    cur_coil['SubCoils'].append(cur_dict)                    
                case 'Shell':
                    cur_dict['SubShells'] = []
                    cur_shell = cur_dict
                    top_dict['Shells'].append(cur_dict)
                case 'SubShell':
                    cur_shell['SubShells'].append(cur_dict)
                case 'Separatrix':
                    top_dict['Separatricies'].append(cur_dict)
                case 'Measure':
                    top_dict['Measures'].append(cur_dict)
                case 'Limiter':
                    top_dict['Limiters'].append(cur_dict)
                case 'CodeControl':
                    cur_dict = top_dict
                case _:
                    cur_dict = top_dict[key]
        else:
            sub_key, value = line.split('=')
            sub_key = sub_key.strip()
            value = value.strip()
            cur_dict[sub_key] = value

    return top_dict


