from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from fcswrite import write_fcs
from loguru import logger
from tqdm.autonotebook import tqdm
from rich import print as rprint


def take_out_the_trash(
    adata: ad.AnnData,
    method: Literal["sd", "quantile"] = "sd",
    inplace: bool = True,
    **kwargs,
) -> ad.AnnData:

    if not inplace:
        adata = adata.copy()

    if method == "sd":
        filter_limits = {
            "140Ce_Beads": 1,
            "Residual": 2,
            "Center": 2,
            "Offset": 1,
            "Width": 2,
            "Event_length": 2,
            "103Rh_Viability": 2,
            "191Ir_DNA1": 2,
            "193Ir_DNA2": 1,
        }
        removed_counts = {k:0 for k in filter_limits}
    else:
        filter_limits = {
            "140Ce_Beads": {"lower": 0.000, "upper": 0.995},
            "Residual": {"lower": 0.001, "upper": 0.900},
            "Center": {"lower": 0.020, "upper": 0.980},
            "Offset": {"lower": 0.020, "upper": 0.990},
            "Width": {"lower": 0.005, "upper": 0.995},
            "Event_length": {"lower": 0.000, "upper": 0.960},
            "103Rh_Viability": {"lower": 0.000, "upper": 0.985},
            "191Ir_DNA1": {"lower": 0.010, "upper": 0.960},
            "193Ir_DNA2": {"lower": 0.015, "upper": 0.985},
        }
        removed_counts = {k:0 for k in filter_limits}

    logger.debug(f"using following as cutoffs: {filter_limits}")
    df = sc.get.obs_df(adata, keys=list(filter_limits.keys()))

    initial_number_of_events = df.shape[0]
    logger.debug(f"starting number of events: {df.shape[0]}")

    for marker in tqdm(filter_limits):
        current_count = df.shape[0]
        if method == "sd":
            num_sd = filter_limits[marker]
            upper = np.mean(df[marker]) + (num_sd * np.std(df[marker]))
            lower = np.mean(df[marker]) - (num_sd * np.std(df[marker]))
            df = df.loc[((df[marker] >= lower) & (df[marker] <= upper)), :]
        elif method == "quantile":
            df = df.loc[
                (
                    (
                        df[marker]
                        >= np.quantile(df[marker], filter_limits[marker]["lower"])
                    )
                    & (
                        df[marker]
                        <= np.quantile(df[marker], filter_limits[marker]["upper"])
                    )
                ),
                :,
            ]
        removed_counts[marker] = current_count-df.shape[0]

    logger.debug(f"removed {removed_counts}"
    )
    logger.debug(f"removed {initial_number_of_events - df.shape[0]} total")
    logger.debug(f"{100*((initial_number_of_events - df.shape[0])/initial_number_of_events):.2f}%"
    )
    rprint(f"removed [red]{initial_number_of_events - df.shape[0]}[/red] overall, or "
        f"[orange]{100*((initial_number_of_events - df.shape[0])/initial_number_of_events):.2f}%[/orange] of the total")

    
    return adata[df.index, :].copy()


def export_fcs(filename: str, df: pd.DataFrame, original_fcs_channels: list[str]) -> None:

    write_fcs(
        filename=filename,
        chn_names=list(original_fcs_channels),
        data=df.loc[:, original_fcs_channels].to_numpy(),
        compat_chn_names=False,
    )
