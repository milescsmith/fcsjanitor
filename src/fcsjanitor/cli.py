from enum import Enum
from pathlib import Path
from sys import stderr
from typing import Optional
from datetime import datetime

import fcsparser
import pytometry as pm
import scanpy as sc
import typer
from loguru import logger
from rich.console import Console
from rich.traceback import install
from tqdm.autonotebook import tqdm

from fcsjanitor import __version__
from fcsjanitor.janitor import export_fcs, take_out_the_trash

install(show_locals=True)

logger.remove()

console = Console()


class FilterMethod(str, Enum):
    sd = "sd"
    quantile = "quantile"


def version_callback(value: bool) -> None:
    """Prints the version of the package."""
    if value:
        console.print(f"[yellow]fcsjanitor[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


app = typer.Typer(
    name="fcsjanitor",
    help=(
        "Utility to automatically clean junk from CyTOF FCS files according to Standard Biotools Recommendations"
    ),
    add_completion=False,
    rich_markup_mode="markdown",
)


def unpack(iterable):
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
    files: list[Path] = typer.Option(
        ...,
        "-f",
        "--files",
        help="Either a list of files or directories in which to look for FCS files",
    ),
    output_dir: Optional[Path] = typer.Option(
        Path.cwd(), "-o", "--output", help="Location to write cleaned FCS files"
    ),
    method: FilterMethod = typer.Option(
        FilterMethod.sd, "-m", "--method", help="Method to use when filtering events"
    ),
    suffix: Optional[str] = typer.Option(
        "_cleaned", "-s", "--suffix", help="Suffix to append to the new files"
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Clean one or more CyTOF FCS files."""

    # if len(extra.args) > 0:
    #     extra_args = parse_extras(extra.args)

    if verbose:
        logger.add(stderr, level="DEBUG")
        logger.add(f"janitor_{datetime.now().strftime('%d-%m-%Y--%H-%M-%S')}.log", level="DEBUG")
    else:
        logger.add(stderr, level="ERROR")
        logger.add(f"janitor_{datetime.now().strftime('%d-%m-%Y--%H-%M-%S')}.log", level="ERROR")

    filelist = unpack(files)

    logger.debug(f"{filelist=}")

    for i in tqdm(filelist):
        original_fcs_meta, original_fcs_data = fcsparser.parse(i, reformat_meta=True)
        # TODO: I don't really things pytometry is needed here.  Will remove later, 
        # but for now this works.
        adata = pm.io.read_fcs(i)
        pm.pp.split_signal(
            adata=adata, var_key="channel", option="element", data_type="cytof"
        )

        take_out_the_trash(adata, method=method)  # , **extra_args)

        df = sc.get.obs_df(
            adata=adata,
            keys=adata.var_names.to_list()
            + ["Time", "Event_length", "Center", "Offset", "Width", "Residual"],
        )

        if suffix is not None:
            newfilename = i.parent.joinpath(f"{i.stem}_{suffix}{i.suffix}")
        else:
            newfilename = i

        export_fcs(
            filename=newfilename,
            df=df,
            original_fcs_channels=list(original_fcs_meta["_channel_names_"]),
        )
