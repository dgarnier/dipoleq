"""Nox sessions."""

import argparse
import os
import shlex
from pathlib import Path
from textwrap import dedent
from typing import Any

from nox import Session, options
from nox_uv import session

options.default_venv_backend = "uv"

package = "dipoleq"
python_versions = ["3.13", "3.12", "3.11", "3.10"]
options.sessions = (
    "pre-commit",
    "mypy",
    "tests",
    "coverage",
    "typeguard",
)


def activate_virtualenv_in_precommit_hooks(session: Session) -> None:
    """Activate virtualenv in hooks installed by pre-commit.

    This function patches git hooks installed by pre-commit to activate the
    session's virtual environment. This allows pre-commit to locate hooks in
    that environment when invoked from git.

    Args:
        session: The Session object.
    """
    assert session.bin is not None  # nosec

    # Only patch hooks containing a reference to this session's bindir. Support
    # quoting rules for Python and bash, but strip the outermost quotes so we
    # can detect paths within the bindir, like <bindir>/python.
    bindirs = [
        bindir[1:-1] if bindir[0] in "'\"" else bindir
        for bindir in (repr(session.bin), shlex.quote(session.bin))
    ]

    virtualenv = session.env.get("VIRTUAL_ENV")
    if virtualenv is None:
        return

    headers = {
        # pre-commit < 2.16.0
        "python": f"""\
            import os
            os.environ["VIRTUAL_ENV"] = {virtualenv!r}
            os.environ["PATH"] = os.pathsep.join((
                {session.bin!r},
                os.environ.get("PATH", ""),
            ))
            """,
        # pre-commit >= 2.16.0
        "bash": f"""\
            VIRTUAL_ENV={shlex.quote(virtualenv)}
            PATH={shlex.quote(session.bin)}"{os.pathsep}$PATH"
            """,
        # pre-commit >= 2.17.0 on Windows forces sh shebang
        "/bin/sh": f"""\
            VIRTUAL_ENV={shlex.quote(virtualenv)}
            PATH={shlex.quote(session.bin)}"{os.pathsep}$PATH"
            """,
    }

    hookdir = Path(".git") / "hooks"
    if not hookdir.is_dir():
        return

    for hook in hookdir.iterdir():
        if hook.name.endswith(".sample") or not hook.is_file():
            continue

        if not hook.read_bytes().startswith(b"#!"):
            continue

        text = hook.read_text()

        if not any(
            (Path("A") == Path("a") and bindir.lower() in text.lower())
            or bindir in text
            for bindir in bindirs
        ):
            continue

        lines = text.splitlines()

        for executable, header in headers.items():
            if executable in lines[0].lower():
                lines.insert(1, dedent(header))
                hook.write_text("\n".join(lines))
                break


@session(name="pre-commit", python=python_versions)
def precommit(session: Session) -> None:
    """Lint using pre-commit."""
    args = session.posargs or [
        "run",
        "--all-files",
        "--hook-stage=manual",
        "--show-diff-on-failure",
    ]
    session.install(
        "bandit",
        "darglint",
        "flake8",
        "flake8-bugbear",
        "flake8-docstrings",
        "flake8-rst-docstrings",
        "pep8-naming",
        "pre-commit",
        "pre-commit-hooks",
        "pyupgrade",
        "ruff",
        "numpy",
    )
    session.run("pre-commit", *args)
    if args and args[0] == "install":
        activate_virtualenv_in_precommit_hooks(session)


@session(
    python=python_versions,
    uv_groups=["test", "mypy"],
)
def mypy(session: Session) -> None:
    """Type-check using mypy."""
    args = session.posargs or ["python/dipoleq"]

    groups = ["test, mypy"]
    if _should_install_imas(session.python):
        groups.append("imas")

    install_cmd = ["uv", "sync", "--python", session.python]
    for group in groups:
        install_cmd.extend(["--group", group])
    session.run(*install_cmd)

    session.run("mypy", *args)
    # run mypy on noxfile.py ?  why?
    # if not session.posargs:
    #    session.run("mypy", f"--python-executable={sys.executable}", "noxfile.py")


@session(python=python_versions)
def tests(session: Session) -> None:
    """Run the test suite."""

    groups = ["test"]
    if _should_install_imas(session.python):
        groups.append("imas")

    install_cmd = ["uv", "sync", "--python", session.python]
    for group in groups:
        install_cmd.extend(["--group", group])

    session.run(*install_cmd)
    # , env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location}, external=True

    try:
        session.run("coverage", "run", "--parallel", "-m", "pytest", *session.posargs)
    finally:
        if session.interactive:
            session.notify("coverage", posargs=[])


@session(python="3.12", uv_groups=["test"], uv_no_install_project=True)
def coverage(session: Session) -> None:
    """Produce the coverage report."""
    args = session.posargs or ["report"]

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)


@session(python=python_versions)
def typeguard(session: Session) -> None:
    """Runtime type checking using Typeguard."""
    groups = ["test", "typeguard"]
    if _should_install_imas(session.python):
        groups.append("imas")

    install_cmd = ["uv", "sync", "--python", session.python]
    for group in groups:
        install_cmd.extend(["--group", group])

    session.run(*install_cmd)

    session.run("pytest", f"--typeguard-packages={package}", *session.posargs)


@session(uv_groups=["docs"], reuse_venv=True)
def docs(session: Session) -> None:
    """
    Build the docs. Pass "--serve" to serve. Pass "-b linkcheck" to check links.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--serve", action="store_true", help="Serve after building")
    parser.add_argument(
        "-b", dest="builder", default="html", help="Build target (default: html)"
    )
    args, posargs = parser.parse_known_args(session.posargs)

    if args.builder != "html" and args.serve:
        session.error("Must not specify non-HTML builder with --serve")

    # session.install("-e.[docs]", *extra_installs)
    session.chdir("docs")

    if args.builder == "linkcheck":
        session.run(
            "sphinx-build", "-b", "linkcheck", ".", "_build/linkcheck", *posargs
        )
        return

    shared_args = (
        "-n",  # nitpicky mode
        "-T",  # full tracebacks
        f"-b={args.builder}",
        ".",
        f"_build/{args.builder}",
        *posargs,
    )

    if args.serve:
        session.run("sphinx-autobuild", *shared_args)
    else:
        session.run("sphinx-build", "--keep-going", *shared_args)


def _should_install_imas(python_version: Any) -> bool:
    # IMAS-Python doesn't yet support Python 3.14
    return python_version < "3.14" and os.getenv("INSTALL_IMAS", "1") == "1"
