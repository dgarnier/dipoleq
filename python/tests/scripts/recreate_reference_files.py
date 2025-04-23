"""
For tests that compare exactly against reference files, it may be necessary to update
the reference files if the output is changed. This script will regenerate the
reference files used in the tests.
"""

from tests.test_omas import data_dir, make_omas_nc_file

REFERENCE_METHODS = [
    make_omas_nc_file,
]


if __name__ == "__main__":
    for reference_method in REFERENCE_METHODS:
        file_path = reference_method(data_dir / "reference")
        print(f"Wrote reference file: {file_path}")
