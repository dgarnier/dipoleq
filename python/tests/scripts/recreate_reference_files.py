from tests.test_omas import data_dir, make_omas_nc_file

REFERENCE_METHODS = [
    make_omas_nc_file,
]


if __name__ == "__main__":
    for reference_method in REFERENCE_METHODS:
        file_path = make_omas_nc_file(data_dir / "reference")
        print(f"Wrote reference file: {file_path}")
