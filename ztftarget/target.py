#! /usr/bin/env python
#

""" ZTF Target based on fritz """

import pandas
import numpy as np
import  warnings

from ztfquery import fritz
from astropy import time

class FritzTarget( object ):
    def __init__(self, name=None, 
                 lightcurve=None, spectra=None, 
                 source=None, alerts=None):
        """ """
        if name is not None:
            self.set_name(name)
            
        if lightcurve is not None:
            self.set_lightcurve(lightcurve)
            
        if spectra is not None:
            self.set_spectra(spectra)

        if source is not None:
            self.set_source(source)

        if alerts is not None:
            self.set_alerts(alerts)

    @classmethod
    def from_name(cls, targetname, warn=False, store=True, force_dl=False, 
                  incl_source=True, incl_lightcurve=True, incl_spectra=True, incl_alerts=True, 
                  **kwargs):
        """ """
        prop = {**dict(store=store, force_dl=force_dl), **kwargs}
        if incl_lightcurve:
            lightcurve = fritz.FritzPhotometry.from_name(targetname, **prop)
        else:
            lightcurve=None
            
        if incl_spectra:
            spectra = fritz.FritzSpectrum.from_name(targetname, warn=warn, **prop)
        else:
            spectra=None            
            
        if incl_source:
            source = fritz.FritzSource.from_name(targetname, **prop)
        else:
            source=None
            
        if incl_alerts:
            alerts = fritz.FritzAlerts.from_name(targetname, **prop)
        else:
            alerts=None
        
        return cls(name=targetname, 
                   lightcurve=lightcurve, spectra=spectra, 
                    source=source, alerts=alerts)
    
    def store(self, **kwargs):
        """ calls the individual lightcurve and spectra store methods. """
        
        if self.has_lightcurve():
            self.lightcurve.store(**kwargs)
            
        if self.has_source():
            self.source.store(**kwargs)
            
        if self.has_spectra():
            self._call_down_spectra_("store", isfunc=True)

        if self.has_alerts():
            self.alerts.store(**kwargs)

    # ============= #
    #   Methods     #
    # ============= #
    def set_name(self, targetname):
        """ """
        self._name = targetname
        
    def set_lightcurve(self, lightcurve):
        """ """
        self._lightcurve = lightcurve
        
    def set_spectra(self, spectra):
        """ """
        self._spectra = np.atleast_1d(spectra)
        
    def set_source(self, source):
        """ """
        self._source = source
        
    def set_alerts(self, alerts):
        """ """
        self._alerts = alerts
        
    # - derived
    def set_saltresult(self, saltresult, add_residuals=True):
        """ """
        self._saltresult = saltresult
        if add_residuals:
            self._derive_saltresiduals_()
            
    def set_redshift(self, redshift, redshifterr=None):
        """ """
        self._redshift = float(redshift)
        self._redshifterr = float(redshifterr) if redshifterr is not None else redshifterr
    
    # ------- #
    # GETTER  #
    # ------- #
    #
    # - Redshift
    def get_redshift(self, loaded=True, full=False, **kwargs):
        """ Get the object redshift. 
        
        Parameters
        ----------
        loaded: [bool] -optional-
            get the redshift stored at self.redshift (set it using self.set_redshift())
            If you did not use set_redshift, self.redshift is taken from self.source is any.

        full: [bool] -optional-
            Do you want the full fritz details concerning the redshift ?

        - kwargs goes to self.source.get_redshift()

        Returns
        -------
        float (or None)
        """
        if full or len(kwargs)>0:
            if loaded:
                warnings.warn("full or kwargs given, loaded set to False")
            loaded = False

        if loaded:
            return self.redshift
            
        if not self.has_source():
            raise AttributeError("No source set.")
        return self.source.get_redshift(full=full, **kwargs)
    
    #
    # - Classification
    def get_classification(self, full=False, **kwargs):
        """ """
        if not self.has_source():
            raise AttributeError("No source set.")
            
        return self.source.get_classification(full=full, **kwargs)
        
    #
    # - Coordinates
    def get_coordinates(self, full=False, perband=False, filters="*", clearna=True, usestat='mean', **kwargs):
        """ shortcut to self.lightcurve.get_coordinates()
        (see also self.alerts.get_coordinates() or self.source.get_coordinates())
        
        Parameters
        ----------
        clearna: [bool] -optional-
            If the returned coordinates is not a simple Serie, shall this remove the NaN rows ?
        """
        if not self.has_lightcurve():
            raise AttributeError("No lightcurve set, maybe see self.source.get_coordiantes()")
            
        cdata = self.lightcurve.get_coordinates(full=full, perband=perband, filters=filters, 
                                                method=usestat, **kwargs)
        if clearna and type(cdata) == pandas.DataFrame:
            return cdata[~cdata.isna().all(axis=1)]
        return cdata
            
    #
    # - ccd location        
    def get_ccdpos(self, full=False, perband=False, groupby='field', usestat='mean', **kwargs):
        """ shortcut to self.alerts.get_ccdpos()
        
        Parameters
        ----------
        clearna: [bool] -optional-
            If the returned coordinates is not a simple Serie, shall this remove the NaN rows ?
        """
        if not self.has_alerts():
            raise AttributeError("No alerts set.")
            
        return self.alerts.get_ccdpos(full=full, perband=perband, groupby=groupby, 
                                           usestat=usestat, **kwargs)
    
    #
    # - Lightcurves
    def get_lcdata(self, detected=None, filters='*', time_range=None, query=None):
        """ similar to self.lightcurve.get_data() but apply to lcdata (that may incl saltresiduals)"""
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
            query.append("filter == @filters")
            
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
    
    def get_sncosmo_table(self, filters=["ztfr","ztfg","ztfi"],
                              incl_upperlimit=True, zp=25.0,
                              as_astropy=False, **kwargs):
        """ get the lightcurve data in a format sncosmo understands for fitting the lightcurve.
        (used bu fit_salt) """
        from . import utils
        
        det_data = self.get_lcdata(filters=filters, detected=True, **kwargs)
        
        det_fluxes = pandas.DataFrame( utils.get_fluxes(det_data, zp=zp, magkey="mag",
                                                            magerrkey="magerr"), 
                                       index=["flux","fluxerr"],
                                       columns=det_data.index).T
        
        if incl_upperlimit:
            up_data = self.get_lcdata(filters=filters,detected=False)
            upperlimit = pandas.DataFrame( utils.get_upper_limit_fluxes(up_data, zp=zp,
                                                                        upmagkey="limiting_mag"),
                                              index=["flux","fluxerr"],
                                               columns=up_data.index).T
            fluxes = pandas.concat([upperlimit, det_fluxes],axis=0)
        else:
            fluxes = det_fluxes

        # End:
        alldata = self.get_lcdata(filters=filters, detected=None, **kwargs)
        fluxes["time"] = alldata["mjd"].loc[fluxes.index]
        fluxes["band"] = alldata["filter"].loc[fluxes.index]
        fluxes["zp"]   = zp
        fluxes["zpsys"]= "ab"
        if as_astropy:
            from astropy import table
            fluxes = table.Table.from_pandas(fluxes)
        return fluxes

    # ----------- #
    #  Extra      #
    # ----------- #
    def fit_salt(self, incl_upperlimit=True, filterprop={},
                 force_color=None, force_x1=None,
                 force_t0=None, force_redshift=None,
                 fix_redshift=True,
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

        filterprop: [dict or None] -optional-
            kwargs option for get_sncosmo_table() and then get_lcdata()
            
        // fit options

        force_redshift, force_color, force_x1, force_t0: [float or None] -optional-
            Force the redshift, color, stretch (x1) or time of peak (t0) values.
            None, means not forced.
            These options overwrite the fixed dictionary entries
            = force_t0 is in mjd =

        fix_redshift: [bool] -optional-
            shortcut to force_redshift to the target redshift value (see self.redshift).

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
                
        model = salt.get_saltmodel()
        if fixed is None:
            fixed = {}
        if values is None:
            values = {}
            
        if fix_redshift or force_redshift is not None:
            z = self.get_redshift() if force_redshift is None else float(force_redshift)
            fixed["z"] = z

        if force_color is not None:
            fixed["c"] = float(force_color)
            
        if force_x1:
            fixed["x1"] = float(force_x1)
            
        if force_t0:
            fixed["t0"] = float(force_t0)
            
        model.set(**{**values,**fixed})
        parameters = [p_ for p_ in ['z','t0', 'x0', 'x1', 'c'] if p_ not in fixed]
    
        table_ = self.get_sncosmo_table(incl_upperlimit=incl_upperlimit, as_astropy=True, **filterprop)
        
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
        
    # ----------- #
    #  Plot       #
    # ----------- #
    def view_on_fritz(self):
        """ open the source fritz page on the browers.
        = calls self.source.view_on_fritz() =
        """
        if not self.has_source():
            raise AttributeError("No source loaded.")
        return self.source.view_on_fritz()
    
    def show_lightcurve(self, axes=None, incl_salt=True, influx=True,
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
        ZTFCOLOR = { # ZTF
            "ztfr":dict(marker="o",ms=7,  mfc="C3"),
            "ztfg":dict(marker="o",ms=7,  mfc="C2"),
            "ztfi":dict(marker="o",ms=7, mfc="C1")
                }

        
        if not self.has_lightcurve():
            raise AttributeError("No lightcurve set.")
            
        if not incl_salt or not self.has_saltresult():
            return self.lightcurve.show(**kwargs)
        
        
        import matplotlib.pyplot as mpl
        from matplotlib import dates as mdates
        from astropy.time import Time
    
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
                fig = mpl.figure(figsize=[7,5])
        
            left, bottom, width, heigth, resheigth = 0.15,0.1,0.75,0.55, 0.07
            vspan, extra_vspan=0.02, 0
            ax = fig.add_axes([left, bottom+3*(resheigth+vspan)+extra_vspan, width, heigth])
            axres = {'ztfg': fig.add_axes([left, bottom+0*(resheigth+vspan), width, resheigth]),
                     'ztfr': fig.add_axes([left, bottom+1*(resheigth+vspan), width, resheigth]),
                     'ztfi': fig.add_axes([left, bottom+2*(resheigth+vspan), width, resheigth])}
            
        
            bottom_ax = axres["ztfg"]

        allaxes = [axres["ztfg"], axres["ztfr"], axres["ztfi"], ax]
        # - Axes
        #

        # 
        # - Data
        lightcurves = self.get_lcdata(**lcprop)
        bands = np.unique(lightcurves["filter"])
        modeltime, modelbands = self.saltresult.get_lightcurve(bands, as_phase=as_phase, as_dataframe=False,
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
            if band_ not in ZTFCOLOR:
                warnings.warn(f"WARNING: Unknown instrument: {band_} | magnitude not shown")
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
                         **{**base_prop,**ZTFCOLOR[band_]}
                       )
                                    
            # = In Flux
            else:
                y, dy = bdata["flux"], bdata["fluxerr"]
                ax.errorbar(datatime,
                         y,  yerr= dy, 
                         label=band_, 
                         **{**base_prop,**ZTFCOLOR[band_]}
                       )
                
            # - Residual in sigma
            if axres[band_] is not None:
                axres[band_].plot(datatime, 
                                bdata["fluxres"]/bdata["fluxerr"],
                                    marker="o", ls="None", 
                                ms=ZTFCOLOR[band_]["ms"]/2, 
                                mfc=ZTFCOLOR[band_]["mfc"],
                                mec="0.5"
                           )
            # = Models
            if show_model:
                ax.plot(modeltime, modelbands[band_], color=ZTFCOLOR[band_]["mfc"])

        if not influx:
            ax.invert_yaxis()

            # = upperlimit
            for band_ in bands:
                if band_ not in ZTFCOLOR:
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
                                 color=ZTFCOLOR[band_]["mfc"], 
                                 ls="None",  label="_no_legend_")

        # - Plots        
        #
        
        for k_, ax_  in axres.items():
            if ax_ is None:
                continue
            if k_ not in bands:
                ax_.text(0.5,0.5, f"No {k_} data", va="center", ha="center", transform=ax_.transAxes, 
                        color=ZTFCOLOR[k_]["mfc"])
                ax_.set_yticks([])
                ax_.set_xticks([])
            else:
                ax_.set_xlim(*ax.get_xlim())
                ax_.set_ylim(-8,8)
                ax_.axhline(0, **lineprop)
                ax_.axhspan(-2,2, color=ZTFCOLOR[k_]["mfc"], zorder=2, alpha=0.05)
                ax_.axhspan(-5,5, color=ZTFCOLOR[k_]["mfc"], zorder=2, alpha=0.05)            

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
            ax.set_ylabel("flux")
            ax.axhline(0, **lineprop)
        else:
            ax.set_ylabel("mag")
            ax.set_ylim(ylimmag)

        [ax_.set_xlim(*ax.get_xlim()) for ax_ in axres.values() if ax_ is not None]
        [ax_.xaxis.set_ticklabels([]) for ax_ in allaxes if ax_ != bottom_ax and ax_ is not None]

        if show_header:
            ax.set_title(targetname, loc="left", fontsize="medium")
        
            s_ = self.saltresults.get_target_parameters(targetname)
            label = f"x1={s_['x1']:.2f}±{s_['x1_err']:.2f}"
            label+= f" | c={s_['c']:.2f}±{s_['c_err']:.2f}"
            ax.text(1,1, label, va="bottom", ha="right", fontsize="small", color="0.7", 
                        transform=ax.transAxes)

        #
        # - Align
        
        #
        #
        return ax, axres
        


    def show_spectra(self, ax=None, normed=True, 
                     add_restline=[], proprest ={},
                     add_obsline=[], propobs={}, legend=True):
        """ show the target spectra
        """
        import matplotlib.pyplot as mpl
        if not self.has_spectra():
            raise AttributeError("No Spectra loaded.")
        from .utils import read_lines
        if ax is None:
            fig = mpl.figure(figsize=[7,4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        for i,spec_ in enumerate(self.spectra):
            specdata_ = self.specdata.iloc[i]
            time_ = time.Time(specdata_["mjd"], format="mjd").datetime
            day = f"{time_.year:04d}/{time_.month:02d}/{time_.day:02d}"
            label = f"{specdata_['obj_id']} | {day}"
            if "phase" in specdata_:
                label += f" | phase={specdata_.phase:+.1f}"
            
            spec_.show(ax=ax, color=f"C{i}", normed=normed, label=label)

        if add_restline is not None and len(add_restline)>0:
            prop = {**dict(ls="-", color="0.7", zorder=1), **proprest}
            _ = [ax.axvline(l, **prop) for l in read_lines(add_restline)]
            
        if add_obsline is not None and len(add_obsline)>0:
            z = self.get_redshift()
            prop = {**dict(ls="--", color="0.7", zorder=1), **propobs}
            _ = [ax.axvline(l*(1+z), **prop) for l in read_lines(add_obsline)]

        if legend:
            ax.legend(loc="best", frameon=False, fontsize="small")
            
        return ax
    # ============= #
    #  Internal     #
    # ============= #
    def _call_down_spectra_(self, what, isfunc, *args, **kwargs):
        """ """
        if not isfunc:
            return [getattr(s_,what) for s_ in self.spectra]
        return [getattr(s_,what)(*args, **kwargs) for s_ in self.spectra]
    
    def _build_lcdata_(self, keys=["ra","dec", "mjd", "filter", "mag","magerr","magsys","limiting_mag"]):
        """ """
        from . import utils
        lcdata = self.lightcurve.get_data()[keys]
        detected_flux = lcdata[~lcdata["mag"].isna()]
        flux,fluxerr = utils.mag_to_flux(*detected_flux[["mag","magerr"]].values.T)
        fdf = pandas.DataFrame(np.asarray([flux,fluxerr, flux/fluxerr]).T, 
                               index=detected_flux.index, 
                               columns=["flux","fluxerr", "detection"])
        
        
        self._lcdata = pandas.merge(lcdata, fdf, left_index=True, right_index=True, how="outer")

        
    def _build_specdata_(self, keys=["obj_id", "mjd"]):
        """ """
        self._specdata = pandas.DataFrame([self._call_down_spectra_(k, isfunc=False) for k in keys],
                 index=keys).T

    def _derive_saltresiduals_(self):
        """ """
        if self.has_lightcurve():
            self.lcdata["phase"] = self.saltresult.get_phase(self.lcdata["mjd"])

            bands = self.lcdata[["mjd","filter"]]
            model = self.saltresult.get_lightcurve(bands["filter"], jd=bands["mjd"], as_phase=True, as_dataframe=True)
            for filter_ in bands["filter"].unique():
                b_index = self.lcdata.index[self.lcdata["filter"]==filter_]
                self.lcdata.loc[b_index, "fluxres"] = self.lcdata.loc[b_index, "flux"]-model.loc[b_index][filter_]
        if self.has_spectra():
            self.specdata["phase"] = self.saltresult.get_phase(self.specdata["mjd"])


    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def name(self):
        """ """
        return self._name
    
    @property
    def lightcurve(self):
        """ """
        if not hasattr(self,"_lightcurve"):
            return None
        return self._lightcurve
    
    def has_lightcurve(self):
        """ """
        return self.lightcurve is not None
    
    @property
    def lcdata(self):
        """ """
        if not hasattr(self, "_lcdata"):
            if not self.has_lightcurve():
                return None
            self._build_lcdata_()
            
        return self._lcdata
    
    @property
    def source(self):
        """ """
        if not hasattr(self,"_source"):
            return None
        return self._source
    
    def has_source(self):
        """ """
        return self.source is not None
    
    @property
    def spectra(self):
        """ """
        if not hasattr(self,"_spectra"):
            return None
        return self._spectra
    
    def has_spectra(self):
        """ """
        return self.spectra is not None
    
    @property
    def specdata(self):
        """ """
        if not hasattr(self, "_specdata"):
            if not self.has_spectra():
                return None
            self._build_specdata_()
            
        return self._specdata
        
    
    @property
    def nspectra(self):
        """ """
        if not self.has_spectra():
            return 0
        return len(self.spectra)
    
    @property
    def alerts(self):
        """ """
        if not hasattr(self,"_alerts"):
            return None
        return self._alerts
    
    def has_alerts(self):
        """ """
        return self.alerts is not None
    
    # - Derived
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
    def redshift(self):
        """ """
        if not hasattr(self,"_redshift"):
            if self.has_source():
                self.set_redshift( self.get_redshift(loaded=False, full=False) )
            else:
                return None
        return self._redshift

    @property
    def redshifterr(self):
        """ """
        if not hasattr(self,"_redshifterr"):
            return None
        return self._redshifterr
