import numpy
import sncosmo
from ztfquery import filters
filters.load_p48_filters_to_sncosmo(basename="ztf:")


def get_saltmodel():
    """ """
    dust = sncosmo.CCM89Dust()
    return sncosmo.Model("salt2", effects=[dust],
                       effect_names=['MW'],
                       effect_frames=['rest'])


