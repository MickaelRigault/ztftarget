
import os
import pandas
import numpy as np
import warnings
from glob import glob

from ztfquery.filters import load_p48_filters_to_sncosmo as load_ztf_filters_to_sncosmo

PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
FILTERPATH = os.path.join(PACKAGE_PATH, "data/filters")

_NAME_ADAPT = {"p48":"ztf",
              "ioo":"lt"
             }

FILTER_PLTFORMAT = {
            # ZTF
            "ztfr":dict(marker="o", ms=7, mfc="tab:red"),
            "ztfg":dict(marker="o", ms=7, mfc="tab:green"),
            "ztfi":dict(marker="o", ms=7, mfc="tab:orange"),
            # SEDm
            "sedmu":dict(marker="s", ms=6, mfc="tab:purple"),            
            "sedmg":dict(marker="s", ms=6, mfc="darkgreen"),
            "sedmr":dict(marker="s", ms=6, mfc="darkred"),
            "sedmi":dict(marker="s", ms=6, mfc="darkorange"),            
            # LT
            "ioog":dict(marker="d", ms=6, mfc="yellowgreen"),
            "ioor":dict(marker="d", ms=6, mfc="tomato"),
            "iooi":dict(marker="d", ms=6, mfc="goldenrod"),
            "iooz":dict(marker="d", ms=6, mfc="tab:grey"), 
            
            }

    
def load_instruments_filters_to_sncosmo(instruments, **kwargs):
    """ 
    loads all the intrument filters to sncosmo.
    known instruments:
        - lt: (ioo*)
        - sedm: (sedm*)
        - ztf: (ztf*)
    
    a warning is raised if the instrument is not implemented. 

    Returns
    -------
    list
    (list of the successfully loaded instrument)
    """
    instruments = np.atleast_1d(instruments)
    loaded = {}
    for instrument in instruments:
        instrument = instrument.lower()
        if instrument in _NAME_ADAPT:
            instrument_ = _NAME_ADAPT[instrument]
        else:
            instrument_ = instrument

        try:
            filters = eval(f"load_{instrument_}_filters_to_sncosmo")(basename=instrument, **kwargs)
        except NameError:
            print(f"unknown instrument {instrument}")
            warnings.warn(f"unknown instrument {instrument}")
        else:
            loaded[instrument] = [f.name for f in filters]
            
    return loaded


def load_lt_filters_to_sncosmo(bands="*", basename="ioo"):
    """ register the given bands ['g','r', 'i' or 'z'] as basename+`band`"""
    try:
        import sncosmo
    except ImportError:
        raise ImportError("You do not have sncosmo. Please install it (pip install sncosmo)")
    
    bands_ = ["g","r","i","z"] if bands in ["all","*"] else np.atleast_1d(bands)

    bands = []
    for band_ in bands_:
        bandname = f"ioo{band_}"
        inst_path = os.path.join(FILTERPATH, f"{bandname}.dat")
        df_ = pandas.read_csv( inst_path, sep=" ").sort_values("Wavelength")
        if basename is not None and basename != "ioo":
            bandname = bandname.replace("ioo", basename)
            
        band = sncosmo.Bandpass(df_["Wavelength"].values*10, # input in nm
                                df_["transmission"].values/100,# input in %
                                name=bandname
                                )
        sncosmo.registry.register(band, force=True)
        bands.append(band)
        
    return bands

def load_sedm_filters_to_sncosmo(bands="*", basename="sedm"):
    """ register the given bands ['u', 'g','r' or 'i'] as basename+`band`"""
    try:
        import sncosmo
    except ImportError:
        raise ImportError("You do not have sncosmo. Please install it (pip install sncosmo)")
    
    bands_ = ["u","g","r","i"] if bands in ["all","*"] else np.atleast_1d(bands)

    bands = []
    for band_ in bands_:
        bandname = f"sedm{band_}"
        inst_path = os.path.join(FILTERPATH, f"{bandname}.dat")
        df_ = pandas.read_csv( inst_path, sep=" ",
                                   names=["Wavelength",
                                            "transmission"]
                              ).sort_values("Wavelength")
        if basename is not None and basename != "sedm":
            bandname = bandname.replace("sedm", basename)
            
        band = sncosmo.Bandpass(df_["Wavelength"].values, # input in nm
                                df_["transmission"].values/100,# input in %
                                name=bandname
                                )
        sncosmo.registry.register(band, force=True)
        bands.append(band)
        
    return bands


