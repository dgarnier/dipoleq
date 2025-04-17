from tests.test_omas import make_omas_nc_file, data_dir


REFERENCE_METHODS = [
    make_omas_nc_file,
]


if __name__ == '__main__':
    for reference_method in REFERENCE_METHODS:
        file_path = make_omas_nc_file(data_dir / "reference")
        print(f'Written reference file: {file_path}')
