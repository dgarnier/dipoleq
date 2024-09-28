# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from __future__ import annotations

import importlib.metadata
import os
from pathlib import Path

import sphinx_rtd_theme

# Define the canonical URL if you are using a custom domain on Read the Docs
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True

project = "DipolEq"
copyright = "2024, Darren Garnier"
author = "Darren Garnier"
version = release = importlib.metadata.version("dipoleq")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",  # napoleon is required for autodoc to work with google-style docstrings and numpy-style docstrings
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinxcontrib.autodoc_pydantic",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

autosummary_generate = True
autodoc_docstring_signature = True
autodoc_pydantic_model_show_json = True
autodoc_pydantic_settings_show_json = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "pydata_sphinx_theme"
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# html_theme = ""
# html_static_path = ["_static"]

nitpick_ignore = []

with Path("nitpick-exceptions").open("r", encoding="utf-8") as f:
    for line in f.readlines():
        if line.strip() == "" or line.startswith("#"):
            continue
        dtype, target = line.split(None, 1)
        target = target.strip()
        nitpick_ignore.append((dtype, target))

print(nitpick_ignore)
always_document_param_types = True

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "pydantic": ("https://docs.pydantic.dev/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "h5py": ("https://docs.h5py.org/en/latest/", None),
}
