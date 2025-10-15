from pathlib import Path
from typing import Dict

import geopandas as gpd
import pandas as pd


def load_chicago_shootings() -> Dict[str, object]:
    """
    Load Chicago Shootings dataset from the bundled data directory.

    This mirrors BSTPP.bstpp.main.load_Chicago_Shootings but reads files
    directly from the repository instead of using pkgutil resources.

    Returns
    -------
    dict
        events_2022: pd.DataFrame with 2022 events (columns include X, Y, T)
        events_2023: pd.DataFrame with 2023 events (columns include X, Y, T)
        covariates:  gpd.GeoDataFrame with spatial covariates
        boundaries:  gpd.GeoDataFrame with community area boundaries
    """
    # Resolve repository root from this file location: repo/utils/chicago_data.py
    # src/utils/chicago_data.py -> parents[2] == repo root
    repo_root = Path(__file__).resolve().parents[2]
    # Use the top-level data directory (no BSTPP dependency)
    # Primary location
    data_dir = repo_root / "data"
    events_2022_path = data_dir / "Chicago_2022_xyt.csv"
    events_2023_path = data_dir / "Chicago_2023_xyt.csv"
    covariates_path = data_dir / "Chicago_cov.zip"
    boundaries_path = data_dir / "Boundaries - Community Areas (current).zip"

    # Fallback to BSTPP/bstpp/data if top-level data/ not found
    if not (
        events_2022_path.exists()
        and events_2023_path.exists()
        and covariates_path.exists()
        and boundaries_path.exists()
    ):
        alt = repo_root / "BSTPP" / "bstpp" / "data"
        e2022 = alt / "Chicago_2022_xyt.csv"
        e2023 = alt / "Chicago_2023_xyt.csv"
        cov = alt / "Chicago_cov.zip"
        bnd = alt / "Boundaries - Community Areas (current).zip"
        if e2022.exists() and e2023.exists() and cov.exists() and bnd.exists():
            events_2022_path = e2022
            events_2023_path = e2023
            covariates_path = cov
            boundaries_path = bnd
        else:
            missing = [
                str(p)
                for p in [
                    events_2022_path,
                    events_2023_path,
                    covariates_path,
                    boundaries_path,
                ]
                if not p.exists()
            ]
            alt_missing = [str(p) for p in [e2022, e2023, cov, bnd] if not p.exists()]
            raise FileNotFoundError(
                f"Required data files not found. Tried {data_dir} (missing={missing}) and {alt} (missing={alt_missing})."
            )

    events_2022 = pd.read_csv(events_2022_path)
    events_2023 = pd.read_csv(events_2023_path)
    covariates = gpd.read_file(str(covariates_path))
    boundaries = gpd.read_file(str(boundaries_path))

    return {
        "events_2022": events_2022,
        "events_2023": events_2023,
        "covariates": covariates,
        "boundaries": boundaries,
    }
