#! /usr/bin/env python
#

""" ZTF Target Lightcurve tools """

import warnings
import numpy as np

from .config import ZTFCOLOR

class ZTFLightCurve( object ):
    """ LightCurve """
    def __init__(self, dataframe, name=None):
        """ LightCurve """
        self._name = name
        self.set_data(dataframe)
        
    @classmethod
    def from_name(cls, name, source="marshal", download=True, update=False, **kwargs):
        """ """
        if source in ["marshal"]:
            from ztfquery import marshal
            lc = marshal.get_target_lightcurve(name, download=download, update=update, **kwargs)
        else:
            raise NotImplementedError(f"source '{source}' not implemented")
            
        return cls(lc, name=name)
    # ----------------- #
    #   Methods         #
    # ----------------- #    
    # -------- #
    #  SETTER  #
    # -------- #
    def set_data(self, dataframe):
        """ """
        dataframe["inst_filter"] = [d.split("+")[-1].replace('"',"").lower()
                                       for d in dataframe["instrument"]+":"+dataframe["filter"]
                                      ]
        
        if 'magpsf' in dataframe.columns:
            keys = {'filter':'inst_filter',
                    'mag':'magpsf',
                    'mag.err':'sigmamagpsf',
                    'upmag':'limmag',
                    'jdobs':'jdobs'}
        else:
            keys = {'filter':'inst_filter',
                    'mag':'mag',
                    'mag.err':'emag',
                    'upmag':'limmag',
                    'jdobs':'jdobs'}
            
        for k,v in keys.items():
            if k != v:
                dataframe[k] = dataframe.pop(v)
                
        self._data = dataframe
        
    # -------- #
    #  GETTER  #
    # -------- #    
    # -------- #
    #  PLOTTER #
    # -------- #
    def show(self, ax=None, **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates
        from astropy.time import Time        
        
        if ax is None:
            fig = mpl.figure(figsize=[5,3])
            ax = fig.add_axes([0.15,0.15,0.75,0.75])
        else:
            fig = ax.figure
        
        
        base_prop = dict(ls="None", mec="0.9", mew=0.5, ecolor="0.7")
        # DataPoints
        for filter_ in np.unique(self.data["filter"]):
            if filter_ not in  ZTFCOLOR:
                print("WARNING: Unknown instrument: %s | magnitude not shown"%filter_)
                continue
            
            jd, mag, magerr = self.data[self.data["filter"].isin([filter_]) & 
                                       ~self.data["mag"].isin([99.00])][
                                        ["jdobs","mag","mag.err"]
                                        ].values.T
        
            ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jd], 
                     mag, yerr= magerr, 
                     label="%s"%filter_, **{**base_prop,**ZTFCOLOR[filter_.replace('"',"")]})
            
        # Upper Limits       
        
        ax.invert_yaxis()  
        
        for filter_ in np.unique(self.data["filter"]):
            if filter_ not in  ZTFCOLOR:
                print("WARNING: Unknown instrument: %s | upper limits not shown"%filter_)
                continue

            jdup, upmag = self.data[self.data["filter"].isin([filter_]) & 
                                 self.data["mag"].isin([99.00])][
                            ["jdobs", "upmag"]
                        ].values.T
            ax.errorbar([Time(jd_, format="jd").datetime for jd_ in jdup], 
                                upmag, yerr=0.15, lolims=True,alpha=0.3,
                                    color=ZTFCOLOR[filter_.replace('"',"")]["mfc"], 
                            ls="None", 
                                    label="_no_legend_")
        
    
        # Data locator
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.set_ylabel("mag")
        
        return ax
    
    # ----------------- #
    #   Properties      #
    # ----------------- #    
    @property
    def name(self):
        """ target name """
        return self._name
    
    @property
    def data(self):
        """ LightCurve dataframe """
        return self._data
