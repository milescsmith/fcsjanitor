from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from anndata.experimental.backed import Dataset2D
from loguru import logger
from rich import print as rprint
from tqdm.auto import tqdm


class ChannelNotFoundError(Exception):
    pass


def take_out_the_trash(
    adata: ad.AnnData,
    method: Literal["sd", "quantile"] = "sd",
    inplace: bool = True,
) -> ad.AnnData:
    if not inplace:
        adata = adata.copy()

    if "140Ce_Beads" in adata.var_names:
        bead_marker= "140Ce_Beads"
    elif "140Ce_Bead" in adata.var_names:
        bead_marker = "140Ce_Bead"
    else:
        logger.exception(
            "No channels matched the expected name for the 140Ce bead marker. It should be '140Ce_Beads' or '140Ce_Bead'"
        )
        raise ChannelNotFoundError

    if "103Rh_Viability" in adata.var_names:
        viability_marker = "103Rh_Viability"
    elif "103Rh_Live_Dead" in adata.var_names:
        viability_marker = "103Rh_Live_Dead"
    else:
        logger.exception(
            "No channels matched the expected name for the viability marker. It should be '103Rh_Viability' or '103Rh_Live_Dead'"
        )
        raise ChannelNotFoundError

    # TODO: need to peek into the adata objects and see if they have 103Rh_Viability or 103Rh_Live_Dead
    # Same for 140CE_Beads vs 140CE_Bead
    if method == "sd":
        filter_limits_sd: dict[str, int] = {
            bead_marker: 1,
            "Residual": 2,
            "Center": 2,
            "Offset": 1,
            "Width": 2,
            "Event_length": 2,
            viability_marker: 2,
            "191Ir_DNA1": 2,
            "193Ir_DNA2": 1,
        }
        filter_limits_keys: list[str] = list(filter_limits_sd.keys())
        logger.debug(f"using following as cutoffs: {filter_limits_sd}")
    else:
        filter_limits_quant: dict[str, dict[str, float]] = {
            bead_marker: {"lower": 0.000, "upper": 0.995},
            "Residual": {"lower": 0.001, "upper": 0.900},
            "Center": {"lower": 0.020, "upper": 0.980},
            "Offset": {"lower": 0.020, "upper": 0.990},
            "Width": {"lower": 0.005, "upper": 0.995},
            "Event_length": {"lower": 0.000, "upper": 0.960},
            viability_marker: {"lower": 0.000, "upper": 0.985},
            "191Ir_DNA1": {"lower": 0.010, "upper": 0.960},
            "193Ir_DNA2": {"lower": 0.015, "upper": 0.985},
        }
        filter_limits_keys = list(filter_limits_quant.keys())
        logger.debug(f"using following as cutoffs: {filter_limits_quant}")
    removed_counts = dict.fromkeys(filter_limits_keys, 0)
    df: Dataset2D | pd.DataFrame = adata.obs[filter_limits_keys]

    initial_number_of_events = df.shape[0]
    logger.debug(f"starting number of events: {initial_number_of_events}")

    for marker in tqdm(filter_limits_keys):
        current_count = df.shape[0]
        if method == "sd":
            num_sd = filter_limits_sd[marker]
            upper = np.mean(df[marker]) + (num_sd * np.std(df[marker]))
            lower = np.mean(df[marker]) - (num_sd * np.std(df[marker]))
            df = df.loc[((df[marker] >= lower) & (df[marker] <= upper)), :]
        elif method == "quantile":
            df = df.loc[
                (
                    (df[marker] >= np.quantile(df[marker], filter_limits_quant[marker]["lower"]))
                    & (df[marker] <= np.quantile(df[marker], filter_limits_quant[marker]["upper"]))
                ),
                :,
            ]
        removed_counts[marker] = current_count - df.shape[0]

    logger.debug(f"removed {removed_counts}")
    logger.debug(f"removed {initial_number_of_events - df.shape[0]} total")
    logger.debug(f"{100 * ((initial_number_of_events - df.shape[0]) / initial_number_of_events):.2f}%")
    rprint(
        f"removed [red]{initial_number_of_events - df.shape[0]}[/red] overall, or [orange]{100 * ((initial_number_of_events - df.shape[0]) / initial_number_of_events):.2f}%[/orange] of the total"
    )

    return adata[df.index, :].copy()
