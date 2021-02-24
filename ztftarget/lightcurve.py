#! /usr/bin/env python
#

""" ZTF Target Lightcurve tools """

import warnings
import numpy as np
import pandas
        
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
    def get_detection_timerange(self, filters=None, **kwargs):
        """ Get the first and latest lightcurve data point 
        
        Parameters
        ----------
        filter: [string or list of] -optional-
            Provide the filter you want to keep.
            could be:
                - None: all filters will be considered
                - filter name (or list of): only the given filter(s) will be used

        **kwargs goes to get_sorted_data() -> get_data()
        
        Returns
        -------
        DataFrame
        """
        return self.get_sorted_data(filters=filters, **kwargs).iloc[[0,-1]]
    
    
    def get_first_detection(self, n=1, filters=None, **kwargs):
        """ Get the data entry corresponding to the first n detections 

        Parameters
        ----------
        n: [int] -optional-
            number of returned entries. 1 means the first, 2 means the first two, etc.

        filter: [string or list of] -optional-
            Provide the filter you want to keep.
            could be:
                - None: all filters will be considered
                - filter name (or list of): only the given filter(s) will be used

        **kwargs goes to get_sorted_data() -> get_data()
        Returns
        -------
        DataFrame
        """
        return self.get_sorted_data(filters=filters, sort_by="jdobs", **kwargs).iloc[:n]
    
    def get_latest_detection(self, n=1, filters=None, **kwargs):
        """ Get the data entry corresponding to the latest n detections 

        Parameters
        ----------
        n: [int] -optional-
            number of returned entries. 1 means latest, 2 means the two latest, etc.

        filter: [string or list of] -optional-
            Provide the filter you want to keep.
            could be:
                - None: all filters will be considered
                - filter name (or list of): only the given filter(s) will be used
        
        Returns
        -------
        DataFrame
        """
        return self.get_sorted_data(filters=filters, sort_by="jdobs", *kwargs).iloc[-n:]

    def get_sorted_data(self, filters=None, detection=True, sort_by="jdobs", **kwargs):
        """ Get data entries sorted by the given key

        Parameters
        ----------
        filters: [string or list of] -optional-
            Provide the filter you want to keep.
            could be:
                - None: all filters will be considered
                - filter name (or list of): only the given filter(s) will be used
        
        detection: [string/bool] -optional-
            The kind of detection you want to keep:
            could be:
                - None: no selection on the detection
                - True: only detection datapoints
                - False: only non-detection datapoints.

        sort_by: [string] -optional-
            column used for sorting the returned dataframe
            
        Returns
        -------
        DataFrame
        """
        return self.get_data(filters=filters, detection=detection, **kwargs).sort_values(sort_by)

    def get_data(self, filters=None, programid=None, detection=None, **kwargs):
        """ get a filtered version of the data.
        
        Parameters
        ----------
        filter: [string or list of] -optional-
            Provide the filter you want to keep.
            could be:
                - None: all filters will be considered
                - filter name (or list of): only the given filter(s) will be used

        programid: [string or list of] -optional-
            Provide the programid you want to keep.
            could be:
                - None: all programid will be considered
                - filter name (or list of): only the given programid(s) will be used
        
        detection: [string/bool] -optional-
            The kind of detection you want to keep:
            could be:
                - None: no selection on the detection
                - True: only detection datapoints
                - False: only non-detection datapoints.
        
        **kwargs will behave as filter and programid. provide a column and accepted value(s).
        Returns
        -------
        DataFrame
        """        
        filtering = []
        
        # Filters
        if filters is not None:
            filters = np.atleast_1d(filters)
            filtering.append("filter in @filters")
        # ProgramID
        if programid is not None:
            programid = np.atleast_1d(programid)
            filtering.append("programid in @programid")

        # Any            
        for k,v in kwargs.items():
            if k not in self.data.columns:
                warnings.warn(f"unknown columns {k} - ignored.")
                continue
            if v is not None:
                v = np.atleast_1d(v)
                filtering.append(f"{k} in @v")
                
        # Detection                         
        if detection is not None:
            filtering.append("absmag < 99" if detection else "absmag == 99")

        # - output
        if len(filtering)==0:
            return self.data.copy()
        
        return self.data.query(" & ".join(filtering))
    
    def get_filtered(self, filters=None, programid=None, detection=None, **kwargs):
        """ DEPRECATED, named get_data() now
        """        
        return self.get_data(filters=filters, programid=programid, detection=detection, **kwargs)


    def get_sncosmo_table(self, filters=["ztf:r","ztf:g","ztf:i"], incl_upperlimit=True, zp=25.0,
                              as_astropy=False, **kwargs):
        """ """
        from . import utils
        
        det_data = self.get_data(filters, detection=True, **kwargs)
        
        det_fluxes = pandas.DataFrame( utils.get_fluxes(det_data, zp=zp), 
                                       index=["flux","fluxerr"],
                                       columns=det_data.index).T
        
        if incl_upperlimit:
            from .utils import get_upper_limit_fluxes            
            up_data = self.get_data(filters=["ztf:r","ztf:g","ztf:i"],detection=False)
            upperlimit = pandas.DataFrame( utils.get_upper_limit_fluxes(up_data, zp=zp),
                                              index=["flux","fluxerr"],
                                               columns=up_data.index).T
            fluxes = pandas.concat([upperlimit, det_fluxes],axis=0)
        else:
            fluxes = det_fluxes

        # End:
        alldata = self.get_data(filters, detection=None, **kwargs)
        fluxes["time"] = alldata["jdobs"].loc[fluxes.index]
        fluxes["band"] = alldata["filter"].loc[fluxes.index]
        fluxes["zp"]   = zp
        fluxes["zpsys"]= "ab"
        if as_astropy:
            from astropy import table
            fluxes = table.Table.from_pandas(fluxes)
        return fluxes


    def fit_salt(self, incl_upperlimit=True, 
             fixed=None, values=None, bounds=None,
             filterprop={},
             **kwargs):
        """ 
        Returns
        -------
        (result, fitted_model), datatable
        """
        from ztftarget import salt
        import sncosmo
    
        model = salt.get_saltmodel()
        if fixed is None:
            fixed = {}
        if values is None:
            values = {}
        
        model.set(**{**values,**fixed})
        parameters = [p_ for p_ in ['z','t0', 'x0', 'x1', 'c'] if p_ not in fixed]
    
        table_ = self.get_sncosmo_table(incl_upperlimit=incl_upperlimit, as_astropy=True, **filterprop)
        return sncosmo.fit_lc(table_, model, 
                              parameters,  # parameters of model to vary
                              bounds=bounds), table_
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
            if filter_ not in ZTFCOLOR:
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
