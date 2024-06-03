"""DipolEq equilibrium solver"""

import tempfile
from argparse import ArgumentParser
from pathlib import Path

from . import Machine
from .h5togeqdsk import h5togeqdsk


def main() -> None:
    """Run the dipoleq solver from the command line using python"""

    parser = ArgumentParser(description="Solve a dipole equilibrium", prog="dipoleq")
    parser.add_argument(
        "input_file",
        metavar="input_file",
        type=str,
        help="configuration, can be .in or .yaml format",
    )
    parser.add_argument(
        "--output_file",
        "-o",
        metavar="[out].[h5|geqdsk|csv]",
        action="append",
        type=str,
        help="output file(s), can be .h5, or .geqdsk format or .csv format\n"
        "or be the same as in input stem if only suffix is given",
        default=[],
    )
    parser.add_argument(
        "--plot",
        "-p",
        action="store_true",
        default=False,
        help="Plot the equilibrium",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        default=False,
        help="Suppress INFO statements while solving",
    )
    args = parser.parse_args()

    input_file = Path(args.input_file).absolute()
    outputs = args.output_file
    if len(outputs) == 0:
        outputs = [".h5", ".geqdsk", ".csv"]
    if len(outputs) == 1 and "." not in outputs[0]:
        o = outputs[0]
        outputs = [o + ex for ex in [".h5", ".geqdsk", ".csv"]]
    else:
        absouts = []
        for output in outputs:
            if output.startswith("."):
                absouts.append(input_file.with_suffix(output))
            else:
                absouts.append(Path(output).absolute())
    odict = {output.suffix: output for output in absouts}

    with tempfile.TemporaryDirectory(prefix=".DipEqTmp") as tmpdir:
        m = Machine.from_file(input_file)
        m.LHname, m.MGname, m.RSname = "", "", ""  # no restart files
        m.solve(quiet=args.quiet)
        tmp_h5 = Path(tmpdir) / "dipEq.h5"
        m.to_hdf5(tmp_h5)
        h5togeqdsk(tmp_h5, plot=args.plot)
        if ".h5" in odict:
            tmp_h5.rename(odict[".h5"])
        if ".geqdsk" in odict:
            tmp_gfile = tmp_h5.with_suffix(".geqdsk")
            if not tmp_gfile.exists():
                for tmp_gfile in tmp_h5.parent.glob("*.geqdsk"):
                    print(tmp_gfile)
            tmp_gfile.rename(odict[".geqdsk"])
        if ".csv" in odict:
            csvstem = odict[".csv"].stem
            for tcsvf in tmp_h5.parent.glob("*.csv"):
                csvf = odict[".csv"]
                if "_" in tcsvf.stem:
                    ext = tcsvf.stem.split("_")[1]
                    csvf = csvf.with_stem(csvstem + "_" + ext)
                tcsvf.rename(csvf)


if __name__ == "__main__":
    main()
