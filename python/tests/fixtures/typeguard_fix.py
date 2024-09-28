"""
OMAS has some bugs with typeguard, so lets exclude it by importing it first.

Yes.. this is a hack, but it works for now.
Need to get OMAS to fix problems with their regular expressions.

"""

import importlib
import sys

from typeguard import TypeguardFinder


def protect_from_typeguard(modules: list[str]) -> None:
    """Protect from typeguard"""

    import_hooks = [hook for hook in sys.meta_path if isinstance(hook, TypeguardFinder)]
    for hook in import_hooks:
        print(f"Removing {hook}")
        sys.meta_path.remove(hook)

    # import so its cached and can't be analyzed
    for module in modules:
        importlib.import_module(module)

    # add them back
    for hook in reversed(import_hooks):
        sys.meta_path.insert(0, hook)


protect_from_typeguard(["omas"])
