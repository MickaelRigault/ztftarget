""" Photometry module. Based on ztfquery.fritz and ztfquery.marshal """

import os
import warnings
import pandas
import numpy as np
from astropy import time

from ztfquery import fritz, marshal




class Lightcurve(  fritz.FritzPhotometry ):
        
    @classmethod
    def from_name(cls, name, force_dl=False, warn=True, store=True, source=["fritz","marshal"], 
                  avoid_duplicate=True,  **kwargs):
        """ """
        sources = []
        speckey = []
        prop = {**dict(force_dl=force_dl, warn=warn), **kwargs}
        
        # - fritz sources        
        if "fritz" in source:
            fritzlc = cls.from_fritz(name, store=store, **prop)
            if fritzlc is not None:
                return fritzlc

        
        # - marshal sources
        if "marshal" in source:
            marshallc = cls.from_marhal(name, **prop) # stored automatically
            if marshallc is not None:
                return marshallc
            
        return None
    
    @classmethod
    def from_fritz(cls, name, warn=True, force_dl=False, store=False, **kwargs):
        """ returns """
        try:
            return super().from_name(name, force_dl=force_dl, store=store, **kwargs)
        except OSError:
            if warn:
                warnings.warn(f"No fritz lightcurve found for {name}")
            return None

    @classmethod
    def from_marhal(cls, name,warn=True, force_dl=False, **kwargs):
        """ """
        try:
            dataframe = marshal.get_target_lightcurve(name, as_fritz=True)
        except:
            if warn:
                warnings.warn(f"No marshal lightcurve found for {name}")
            return None
        return cls(dataframe)
    
    # ============= #
    #  Methods      #
    # ============= #
    def load_intruments_to_sncosmo(self):
        """ """
        from .filters import load_instruments_filters_to_sncosmo
        instruments = self.get_instruments()
        self._sncosmo_instruments = load_instruments_filters_to_sncosmo(instruments)

    
    def set_mwebv(self, mw_ebv):
        """ Set the Milky Way extinction E(B-V)"""
        self._mwebv = float(mw_ebv)

    def set_saltresult(self, saltresult, add_residuals=True):
        """ """
        self._saltresult = saltresult
        if add_residuals:
            self._derive_saltresiduals_()
            
    def get_lcdata(self, detected=None,
                       filters='*', instruments="*",
                       time_range=None, query=None):
        """ similar to self.get_data() but apply to lcdata 
        (that may incl saltresiduals) """
        if query is None:
            query = []
        else:
            query = list(np.atleast_1d(query))
            
        # - Detected Filtering
        if detected is not None:
            if not detected:
                query.append("mag == 'NaN'") 
            else:
                query.append("mag != 'NaN'")

        # - Filters Filtering
        if filters is not None and filters not in ["*","all","any"]:
            filters = list(np.atleast_1d(filters))
            query.append(f"filter in {filters}")

        if instruments is not None and instruments not in ["*","all","any"]:
            instruments = list(np.atleast_1d(instruments))
            query.append(f"instrument_name in {instruments}")
            
        # - Time Filtering
        if time_range is not None:
            tstart, tend = time_range
            if tstart is None and tend is None:
                pass
            else:
                tindex = pandas.DatetimeIndex(time.Time(self.lcdata["mjd"], format="mjd").datetime)
                
                if tstart is None:
                    query.append(f"@tindex<'{tend}'")
                elif tend is None:
                    query.append(f"'{tstart}'<@tindex")
                else:
                    query.append(f"'{tstart}'<@tindex<'{tend}'")
                    
        # - Returns
        if len(query)==0:
            return self.lcdata.copy()
        return self.lcdata.query(" and ".join(query))
        
    def get_sncosmo_table(self, filters=None,
                              instruments="*",
                              incl_upperlimit=True, zp=25.0,
                              as_astropy=False, dropna=True, **kwargs):
        """ get the lightcurve data in a format sncosmo understands for fitting the lightcurve.
        (used bu fit_salt) """
        from . import utils

        if filters is None:
            filters = self.sncosmo_filters
            warnings.warn(f"filters are not specified. Only these (known by sncosmo) will be used {filters}")
            
        
        lcprop = {**dict(filters=filters, instruments=instruments),**kwargs}
        
        det_data = self.get_lcdata(detected=True, **lcprop)
        
        det_fluxes = pandas.DataFrame( utils.get_fluxes(det_data, zp=zp, magkey="mag",
                                                            magerrkey="magerr"), 
                                       index=["flux","fluxerr"],
                                       columns=det_data.index).T
        
        if incl_upperlimit:
            up_data = self.get_lcdata(detected=False, **lcprop)
            upperlimit = pandas.DataFrame( utils.get_upper_limit_fluxes(up_data, zp=zp,
                                                                        upmagkey="limiting_mag"),
                                              index=["flux","fluxerr"],
                                               columns=up_data.index).T
            fluxes = pandas.concat([upperlimit, det_fluxes],axis=0)
        else:
            fluxes = det_fluxes

        # End:
        alldata = self.get_lcdata(detected=None, **lcprop)
        fluxes["time"] = alldata["mjd"].loc[fluxes.index]
        fluxes["band"] = alldata["filter"].loc[fluxes.index]
        fluxes["zp"]   = zp
        fluxes["zpsys"]= "ab"

        if dropna:
            fluxes = fluxes.dropna()
        if as_astropy:
            from astropy import table
            fluxes = table.Table.from_pandas(fluxes)
        return fluxes

    def fit_salt(self, incl_upperlimit=True, filterprop={},
                     filters=["ztfr","ztfg","ztfi"],
                 force_color=None, force_x1=None,
                 force_t0=None, force_redshift=None,
                 fixed=None, values=None, bounds=None,
                 get_object=True,
                 set_it=True,
                 **kwargs):
        """ 
        Parameters
        ----------

        // data options

        incl_upperlimit: [bool] -optional-
            Shall the upper limit be used when fitting salt ?

        filters: [None, str or list of str] -optional-
            list of filters to be used. 
            None (or '*') means no filter selection (so all).
            Careful, filters should be known by sncosmo.

        filterprop: [dict or None] -optional-
            kwargs option for get_sncosmo_table() and then get_lcdata()
            
        
        // fit options

        force_redshift, force_color, force_x1, force_t0: [float or None] -optional-
            Force the redshift, color, stretch (x1) or time of peak (t0) values.
            None, means not forced.
            These options overwrite the fixed dictionary entries
            = force_t0 is in mjd =

        fixed: [dict or None] -optional-
            fix values of the model. Format: e.g. fixed = {'x1':2, 'c':0.1}

        values: [dict or None] -optional-
            provide initial guess values, Format: e.g. values = {'x1':2, 'c':0.1}

        bounds: [dict or None] -optional-
            provide boundaries for each values. Format: e.g. bounds = {'x1':[-3,3], 'c':[-1,3]}

        // returns options

        get_object: [bool] -optional-
            should this method returns the saltresult object or the raw sncomos output ?

        set_it: [bool] -optional-
            if the fit is successful should it be set to this instance ? 

        **kwargs goes to sncosmo.fit_lc()
            
        Returns
        -------
        (result, fitted_model), datatable
        """
        from . import salt
        import sncosmo

        if self.mwebv is None:
            warnings.warn("No MW E(B-V) used. (hence, assumed 0)")
            
        model = salt.get_saltmodel(mwebv=self.mwebv)
        if fixed is None:
            fixed = {}
        if values is None:
            values = {}
            
        if force_redshift is not None:
            z = float(force_redshift)
            fixed["z"] = z

        if force_color is not None:
            fixed["c"] = float(force_color)
            
        if force_x1:
            fixed["x1"] = float(force_x1)
            
        if force_t0:
            fixed["t0"] = float(force_t0)
            
        model.set(**{**values,**fixed})
        parameters = [p_ for p_ in ['z','t0', 'x0', 'x1', 'c'] if p_ not in fixed]
    
        table_ = self.get_sncosmo_table(filters=filters,
                                        incl_upperlimit=incl_upperlimit, as_astropy=True,
                                            **filterprop)
        print(table_)
        #
        # - Fit LightCurve
        (result, fitted_model)= sncosmo.fit_lc(table_, model, 
                                               parameters,  # parameters of model to vary
                                               bounds=bounds, **kwargs)
        saltres = salt.SALTResult(result=result, model=fitted_model, data=table_)
        #
        # - Output
        if set_it:
            self.set_saltresult(saltres)

        if not get_object:
            return (result, fitted_model), table_

        return saltres
    
    # --------- #
    #  Data     #
    # --------- #
    def show(self, filters=None,
                 axes=None, incl_salt=True, influx=True,
                 show_model=True, show_header=False, 
                 ulength=0.1, ualpha=0.1, ylimmag=21.8,
                 fig=None,as_phase=False, lcprop={},
                 **kwargs):
        """ 
        
        Parameters
        ----------
        lcprop: [dict] -optional-
            used as kwargs for get_lcdata()
            
        """
        from .filters import FILTER_PLTFORMAT
        
        if not incl_salt or not self.has_saltresult():
            return super().show(**kwargs)

        # So this include salt.
        from matplotlib import dates as mdates
        from astropy.time import Time


        if filters is None:
            filters = self.sncosmo_filters
        # - Data
        #  data first to known  who many bands
        lightcurves = self.get_lcdata(filters=filters, **lcprop)
        
        data_filters = lightcurves["filter"].unique()
        sncosmo_filters = self.sncosmo_filters # filters from data & known by sncosmo
        #  In case some where cutted out because of selection in get_lcdata():
        #  limit to the one actually in the current lightcurve data
        sncosmo_filters = np.asarray(sncosmo_filters)[np.in1d(sncosmo_filters, data_filters)]
        nfilters = len(sncosmo_filters)
        #
        # - Axes
        if axes is not None:
            if len(np.atleast_1d(axes))==1:
                ax = np.atleast_1d(axes)[0]
                bottom_ax = ax
                axres = {"ztfg":None, "ztfr":None, "ztfi":None}
            elif len(axes) == 4:            
                resg, resr, resi, ax = axes
                axres = {"ztfg":resg, "ztfr":resr, "ztfi":resi}
                bottom_ax = axres["ztfg"]
            else:
                raise ValueError("axes if given must be a single axes or a list of 4 axes [resg, resr, resi, main]")
            fig = ax.figure
            
        else:
            if fig is None:
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=[10,4+nfilters])
                
            ratio = 4/(4+nfilters)
            print(ratio)
            left, bottom, width, heigth, resheigth_tot = 0.15, 0.1, 0.75, ratio, 0.9-ratio
            vspan_tot = resheigth_tot*0.1 # 10%
            resheigth = resheigth_tot*0.9/nfilters # 10%
            vspan = (vspan_tot / (nfilters-1)) if nfilters >1 else vspan_tot
            extra_vspan = 0.01
            # per line
            
            ax = fig.add_axes([left, bottom+resheigth_tot+extra_vspan, width, heigth])
            axres = {f: fig.add_axes([left, bottom+i*(resheigth+vspan), width, resheigth])
                         for i,f in enumerate(sncosmo_filters)}
            
            bottom_ax = axres[sncosmo_filters[0]]

        allaxes = list(axres.values()) + [ax]
        # - Axes
        #

        # 
        bands = np.unique(lightcurves["filter"])
        modeltime, modelbands = self.saltresult.get_lightcurve(bands, as_phase=as_phase,
                                                                as_dataframe=False,
                                                              influx=influx)

        if not as_phase:
            modeltime=Time(modeltime, format="mjd").datetime
            
        # - Data
        #

        #
        # - Properties
        base_prop = dict(ls="None", mec="0.9", mew=0.5, ecolor="0.7")
        lineprop = dict(color="0.7", zorder=1, lw=0.5)
        # - Properties
        #
        
        #
        # - Plots
        for band_ in bands:
            if band_ not in FILTER_PLTFORMAT:
                warnings.warn(f"Plotting issue: unknown instrument: {band_} | magnitude not shown")
                continue
            
            bdata = lightcurves[lightcurves["filter"]==band_]
            if not as_phase:
                datatime = Time(bdata["mjd"], format="mjd").datetime
            else:
                datatime = bdata["phase"]

            # = In Magnitude                
            if not influx:
                flag_notdet = bdata["detection"].isna()
                y, dy = bdata["mag"], bdata["magerr"]
                # detected
                ax.errorbar(datatime[~flag_notdet],
                         y[~flag_notdet],  yerr= dy[~flag_notdet], 
                         label=band_, 
                         **{**base_prop,**FILTER_PLTFORMAT[band_]}
                       )
                                    
            # = In Flux
            else:
                y, dy = bdata["flux"], bdata["fluxerr"]
                ax.errorbar(datatime,
                         y,  yerr= dy, 
                         label=band_, 
                         **{**base_prop,**FILTER_PLTFORMAT[band_]}
                       )
                
            # - Residual in sigma
            if axres[band_] is not None:
                color = FILTER_PLTFORMAT[band_]["mfc"]
                axres[band_].plot(datatime, 
                                bdata["fluxres"]/bdata["fluxerr"],
                                    marker=FILTER_PLTFORMAT[band_]["marker"], ls="None", 
                                ms=FILTER_PLTFORMAT[band_]["ms"]/1.2, 
                                mfc=color,
                                mec="0.5"
                           )
                axres[band_].set_ylabel(band_, fontsize="medium", color=color)
            # = Models
            if show_model:
                ax.plot(modeltime, modelbands[band_], color=FILTER_PLTFORMAT[band_]["mfc"])

        if not influx:
            ax.invert_yaxis()

            # = upperlimit
            for band_ in bands:
                if band_ not in FILTER_PLTFORMAT:
                    warnings.warn(f"WARNING: Unknown instrument: {band_} | magnitude not shown")
                    continue
            
                bdata = lightcurves[lightcurves["filter"]==band_]
                if not as_phase:
                    datatime = Time(bdata["mjd"], format="mjd").datetime
                else:
                    datatime = bdata["phase"]

                flag_notdet = bdata["detection"].isna()
                upmag = bdata["limiting_mag"]
                ax.errorbar(datatime[flag_notdet], upmag[flag_notdet],
                                 yerr=ulength, lolims=True, alpha=ualpha,
                                 color=FILTER_PLTFORMAT[band_]["mfc"], 
                                 ls="None",  label="_no_legend_")

        # - Plots        
        #
        
        for k_, ax_  in axres.items():
            if ax_ is None:
                continue
            if k_ not in bands:
                ax_.text(0.5,0.5, f"No {k_} data", va="center", ha="center",
                             transform=ax_.transAxes, 
                        color=FILTER_PLTFORMAT[k_]["mfc"])
                ax_.set_yticks([])
                ax_.set_xticks([])
            else:
                ax_.set_xlim(*ax.get_xlim())
                ax_.set_ylim(-8,8)
                ax_.axhline(0, **lineprop)
                if k_ in FILTER_PLTFORMAT:
                    color = FILTER_PLTFORMAT[k_]["mfc"]
                else:
                    color = "0.7"
                ax_.axhspan(-2,2, color=color, zorder=2, alpha=0.05)
                ax_.axhspan(-5,5, color=color, zorder=2, alpha=0.05)            

            clearwhich = ["left","right","top"] # "bottom"
            [ax_.spines[which].set_visible(False) for which in clearwhich]
            ax_.tick_params(axis="y", labelsize="small", 
                           labelcolor="0.7", color="0.7")
        # Upper Limits       

        #ax.invert_yaxis()  

        # Data locator
        if not as_phase:
            locator = mdates.AutoDateLocator()
            formatter = mdates.ConciseDateFormatter(locator)
            bottom_ax.xaxis.set_major_locator(locator)
            bottom_ax.xaxis.set_major_formatter(formatter)
        else:
            bottom_ax.set_xlabel("phase [days]")

        if influx:
            ax.set_ylabel("flux", fontsize="large")
            ax.axhline(0, **lineprop)
        else:
            ax.set_ylabel("mag", fontsize="large")
            ax.set_ylim(ylimmag)

        [ax_.set_xlim(*ax.get_xlim()) for ax_ in axres.values() if ax_ is not None]
        [ax_.xaxis.set_ticklabels([]) for ax_ in allaxes if ax_ != bottom_ax and ax_ is not None]

        if show_header:
            ax.set_title(self.name, loc="left", fontsize="medium")
            s_ = self.saltresult.get_parameters()
            label = f"x1={s_['x1']:.2f}±{s_['x1_err']:.2f}"
            label+= f" | c={s_['c']:.2f}±{s_['c_err']:.2f}"
            ax.text(1,1, label, va="bottom", ha="right", fontsize="small", color="k", 
                        transform=ax.transAxes)

        #
        # - Align
        
        #
        #
        return ax, axres

    
    # ============= #
    #  Internal     #
    # ============= #    
    def _build_lcdata_(self, keys=["ra","dec", "mjd", "filter",
                                   "mag","magerr","magsys","limiting_mag",
                                    "instrument_name"]):
        """ """
        from . import utils
        lcdata = self.get_data()[keys].copy()
        lcdata.loc[lcdata.index,"instrument_name"] = lcdata["instrument_name"].str.lower()
        detected_flux = lcdata[~lcdata["mag"].isna()]
        flux,fluxerr = utils.mag_to_flux(*detected_flux[["mag","magerr"]].values.T)
        fdf = pandas.DataFrame(np.asarray([flux,fluxerr, flux/fluxerr]).T, 
                               index=detected_flux.index, 
                               columns=["flux","fluxerr", "detection"])
        
        
        self._lcdata = pandas.merge(lcdata, fdf, left_index=True, right_index=True, how="outer")

    def _derive_saltresiduals_(self):
        """ """
        
        self.lcdata["phase"] = self.saltresult.get_phase(self.lcdata["mjd"])
        self.lcdata["fluxres"] = np.NaN # Nan by default, filling what I can.

        lcdata_ = self.lcdata[self.lcdata["filter"].isin(self.sncosmo_filters)]
        bands = lcdata_[["mjd","filter"]]
        model = self.saltresult.get_lightcurve(bands["filter"], jd=bands["mjd"],
                                                   as_phase=True, as_dataframe=True)
        for filter_ in bands["filter"].unique():
            b_index = self.lcdata.index[lcdata_["filter"]==filter_]
            self.lcdata.loc[b_index, "fluxres"] =  lcdata_.loc[b_index, "flux"]-model.loc[b_index][filter_]

    # ============= #
    #  Internal     #
    # ============= #    
    @property
    def lcdata(self):
        """ """
        if not hasattr(self, "_lcdata"):
            self._build_lcdata_()
            
        return self._lcdata

    @property
    def sncosmo_instruments(self):
        """ """
        if not hasattr(self, "_sncosmo_instruments"):
            self.load_intruments_to_sncosmo()
            
        return list(self._sncosmo_instruments.keys())

    @property
    def sncosmo_filters(self):
        """ """
        if not hasattr(self, "_sncosmo_instruments"):
            self.load_intruments_to_sncosmo()

        return list(np.concatenate( list(self._sncosmo_instruments.values()) ))
    
    @property
    def saltresult(self):
        """ """
        if not hasattr(self, "_saltresult"):
            return None
        return self._saltresult
    
    def has_saltresult(self):
        """ """
        return self.saltresult is not None

    @property
    def mwebv(self):
        """ milky way extinction in the target direction """
        if not hasattr(self, "_mwebv"):
            return None
        return self._mwebv
