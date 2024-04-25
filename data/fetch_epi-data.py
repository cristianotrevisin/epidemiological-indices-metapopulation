#!/usr/bin/env python

"""Download DPC Regioni"""

from collections import defaultdict
import pandas as pd

# GitHub raw.githubusercontent.com
_URL_TEMPLATE = (
    "https://raw.githubusercontent.com/{owner}/{repo}/{branch}/{path}".format
)

URL_PRO = _URL_TEMPLATE(
    owner="pcm-dpc",
    repo="COVID-19",
    branch="master",
    path="dati-regioni/dpc-covid19-ita-regioni.csv",
)

MA_WINDOW = "14D"


def _read_csv(p):
    """read CSV files"""

    frame = pd.read_csv(
        p,
        encoding="UTF-8",
        na_values=[""],
        keep_default_na=False,
        usecols=[
            "data",
            "totale_casi",
            "codice_regione",
            "totale_casi",
        ],
        parse_dates=["data"],
        #dtype=defaultdict(pd.CategoricalDtype) | {"totale_casi": int},
    )

    frame["data"] = frame["data"].dt.tz_localize("Europe/Rome")

    #nuts_columns = ["codice_regione"]
    #nuts_map = (
    #    frame[["codice_regione"] + nuts_columns]
    #    .dropna()
    #    .drop_duplicates()
    #    .sort_values("codice_regione")
    #    .set_index("codice_regione", verify_integrity=True)
    #)

    #for i in nuts_columns:
    #    frame[i] = frame["codice_regione"].map(nuts_map[i])

    #frame.dropna(how="any", subset=nuts_columns, inplace=True)
    frame.dropna(how="any", inplace=True)

    #return frame[["data"] + nuts_columns + ["totale_casi"]]
    return frame[["data"]+ ["codice_regione"]+ ["totale_casi"]]


def get_data():
    frame = _read_csv(URL_PRO)
    frame = frame[frame.totale_casi > 0]

    cases = frame.set_index(["data", "codice_regione"])["totale_casi"].unstack()
    cases = cases.interpolate(method="time", limit_area="inside")

    cases.index = cases.index.tz_convert(None).to_period("1D")

    
    date_range = pd.period_range(start=cases.index.min(), end=cases.index.max())
    full_index = pd.MultiIndex.from_product(
        (date_range, cases.columns), names=["date", "codice_regione"]
    )
    out = pd.DataFrame(index=full_index)
    out["cases"] = cases.unstack().swaplevel()
    out.fillna(0, inplace=True)

    assert out.unstack().stack().equals(out)

    out["new_cases"] = out["cases"].unstack().diff().stack()
    out["new_cases_ma"] = out["new_cases"].unstack().rolling(MA_WINDOW).mean().stack()

    return cases


def main():
    data = get_data()
    data.to_csv("cases_region.csv")

if __name__ == "__main__":
    main()
