from __future__ import annotations
from typing import Iterable

import pandas as pd

def read_gamd_log(log: str) -> pd.DataFrame:
    """reads an AMBER GaMD log file

    Reads an AMBER GaMD log file with standard format, i.e.
    some commented lines starting with #, with the last line
    containing the names of the columns, and then only lines
    with numbers.

    Parameters
    ----------
    log: str
        path to the AMBER GaMD log file.

    Returns
    -------
    df: pd.DataFrame
        pandas dataframe with the log values.
    """
    with open(log, "r") as f:
        prevline = None
        for line in f:
            if line.strip().startswith("#"):
                if prevline is None:
                    prevline = line
                    currline = line
                else:
                    prevline = currline
                    currline = line
            else:
                columns = currline.strip()[1:].split(",")
                columns = [c.strip() for c in columns]

    df = pd.read_csv(log, sep="\s+", skiprows=3, header=None)
    df.columns = columns
    return df

def read_gamd_logs(logs: Iterable[str]) -> pd.DataFrame:
    """reads an arbitray number of AMBER log files.

    Reads an arbitrary number of AMBER GaMD log files with
    standard format, i.e. some commented lines starting
    with #, with the last line containing the names of the
    columns, and then only lines with numbers.

    Parameters
    ----------
    logs: iterable of str
        paths to the AMBER GaMD log files.

    Returns
    -------
    df: pd.DataFrame
        pandas dataframe with the logs values.
    """
    return pd.concat([read_gamd_log(l) for l in logs], axis=0)
