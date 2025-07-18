from collections.abc import Iterable
from datetime import datetime
from enum import Enum
from pathlib import Path
from sys import stderr
from typing import Annotated, Any

import flowkit as fk
import pytometry as pm
import typer
from dateutil import tz
from loguru import logger
from rich.console import Console
from rich.traceback import install
from tqdm.autonotebook import tqdm

from fcsjanitor import __version__
from fcsjanitor.janitor import take_out_the_trash

_ = install(show_locals=True)

logger.remove()

console = Console()

DEFAULT_OUTPUT_DIR = Path().cwd()


class FilterMethod(str, Enum):
    sd = "sd"
    quantile = "quantile"


class OutputFormat(str, Enum):
    anndata = "anndata"
    fcs = "fcs"


def version_callback(show_version: bool) -> None:
    """Prints the version of the package."""
    if show_version:
        console.print(f"[yellow]fcsjanitor[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


app = typer.Typer(
    name="fcsjanitor",
    help=("Utility to automatically clean junk from CyTOF FCS files according to Standard Biotools Recommendations"),
    add_completion=False,
    rich_markup_mode="markdown",
    no_args_is_help=True,
)


def unpack(iterable: Iterable[Any]):
    logger.debug(f"{iterable=}")
    for i in iterable:
        logger.debug(f"{i=}")
        if "__iter__" in dir(i):
            for j in i:
                logger.debug(f"{j=}")
                yield j
        else:
            yield i


@app.command(
    no_args_is_help=True,
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def clean_up_this_mess(
    input_files: Annotated[
        list[Path],
        typer.Option(
            "-i",
            "--input",
            help="Either a list of input files or directories in which to look for FCS files",
        ),
    ],
    output_dir: Annotated[
        Path, typer.Option("-o", "--output", help="Location to write cleaned FCS files")
    ] = DEFAULT_OUTPUT_DIR,
    method: Annotated[
        str, typer.Option("-m", "--method", help="Method to use when filtering events")
    ] = FilterMethod.sd,
    suffix: Annotated[str | None, typer.Option("-s", "--suffix", help="Suffix to append to the new files")] = None,
    output_format: Annotated[
        OutputFormat,
        typer.Option(
            "-f",
            "--format",
            help="Save the cleaned data as an FCS or AnnData object?",
        ),
    ] = OutputFormat.fcs,
    verbose: Annotated[bool, typer.Option("-v", "--verbose")] = False,
    version: Annotated[bool, typer.Option("--version", callback=version_callback)] = False,
) -> None:
    """Clean one or more CyTOF FCS files."""

    # if len(extra.args) > 0:
    #     extra_args = parse_extras(extra.args)

    _ = logger.add(f"janitor_{datetime.now(tz=tz.tzlocal()).strftime('%d-%m-%Y--%H-%M-%S')}.log", level="DEBUG")
    if verbose:
        _ = logger.add(stderr, level="DEBUG")
    else:
        _ = logger.add(stderr, level="ERROR")

    # if not isinstance(input_files, list):
    #     input_files = [input_files]
    for _ in input_files:
        if _.is_dir():
            logger.debug(f"{_} is a dir")
        elif _.is_file():
            logger.debug(f"{_} is a file")
        else:
            logger.debug(f"I don't know what {_} is")

    filelist = unpack([list(_.glob("*.fcs")) if _.is_dir() else _ for _ in input_files])

    logger.debug(f"{filelist=}")

    for i in tqdm(filelist):
        # TODO: I don't really think pytometry is needed here.  Will remove later,
        # but for now this works so it stays.
        adata = pm.io.read_fcs(i)
        pm.pp.split_signal(adata=adata, var_key="channel", option="element", data_type="cytof")
        cleaned_adata = take_out_the_trash(adata, method=method)  # type: ignore

        if suffix is not None:
            newfilename = output_dir.joinpath(f"{i.stem}{suffix}")
        else:
            newfilename = output_dir.joinpath(f"{i.stem}")
        newfilename.parent.mkdir(parents=True, exist_ok=True)

        if output_format == "fcs":
            remove_indices = adata.obs_names.difference(cleaned_adata.obs_names)
            fk_fcs = fk.Sample(i, ignore_offset_error=True)
            fk_fcs.set_flagged_events(remove_indices.astype(int).to_list())
            fk_fcs.export(newfilename.with_suffix(".fcs"), source="raw", exclude_flagged=True)
        elif output_format == "anndata":
            cleaned_adata.write(newfilename.with_suffix(".h5ad"))
