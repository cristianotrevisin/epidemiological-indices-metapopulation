#!/usr/bin/env python3
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
import os



URL = "https://www.istat.it/storage/cartografia/matrici_pendolarismo/matrici_pendolarismo_2011.zip"
LCL = Path("matrici_pendolarismo_2011.zip")


def download(url, path):
    with requests.get(url, stream=True) as resp:
        resp.raise_for_status()
        with open(path, "wb") as fp:
            for chunk in resp.iter_content(chunk_size=None):
                if chunk:
                    fp.write(chunk)



def read_pendo(path):
    mpath = "MATRICE PENDOLARISMO 2011/matrix_pendo2011_10112014.txt"
    with zipfile.ZipFile(path, "r") as zc:
        with zc.open(mpath) as fp:
            zc.extract(mpath)
    print(f"Wrote {mpath}", file=sys.stderr)




def main():
    os.chdir("data/")
    print(f"Download {URL}", file=sys.stderr)
    download(URL, LCL)

    read_pendo(LCL)



if __name__ == "__main__":
    main()