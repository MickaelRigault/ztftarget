#! /usr/bin/env python
#

""" ZTF Target based on fritz """

import pandas
import numpy as np
import warnings

from ztfquery import fritz
from astropy import time

from . import spectroscopy
from . import photometry

class Target( object ):
    def __init__(self, name=None, 
                 lightcurve=None, spectra=None, 
                 source=None, alerts=None,
                     load_mwebv=True):
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

        if load_mwebv:
            self.load_mwebv()
            
    @classmethod
    def from_name(cls, targetname, warn=False, store=True, force_dl=False, 
                  incl_source=True, incl_lightcurve=True, incl_spectra=True, incl_alerts=True,
                  spec_source=["fritz", "sedm"], lc_source=["fritz","marshal"],
                  load_mwebv=True,
                  **kwargs):
        """ """
        prop = {**dict(store=store, force_dl=force_dl), **kwargs}
        if incl_lightcurve:
            lightcurve = photometry.Lightcurve.from_name(targetname, warn=warn, source=lc_source, **prop)
        else:
            lightcurve=None
            load_mwebv = False
            
        if incl_spectra and spec_source is not None and len(spec_source)>0:
            spectra = spectroscopy.Spectrum.from_name(targetname, warn=warn,
                                                          source=spec_source,**prop)
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
                    source=source, alerts=alerts, load_mwebv=load_mwebv)
    
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
        """ shortcut to self.lightcurve.set_saltresult() """
        self.lightcurve.set_saltresult(saltresult, add_residuals=add_residuals)
            
    def set_snidresult(self, snidresult, which):
        """ """
        print("DEPRECATED")
        if not hasattr(self, "_snidresult") or self._snidresult is None:
            self._snidresult = {}
        self._snidresult[which] = snidresult
            
    def set_redshift(self, redshift, redshifterr=None):
        """ """
        self._redshift = float(redshift)
        self._redshifterr = float(redshifterr) if redshifterr is not None else redshifterr

    def set_mwebv(self, mw_ebv, pass_down=True):
        """ Set the Milky Way extinction E(B-V)"""
        self._mwebv = float(mw_ebv)        
        if pass_down and self.has_lightcurve():
            self.lightcurve.set_mwebv(mw_ebv)
            
    # ------- #
    # LOAD    #
    # ------- #
    def load_mwebv(self, sourcemap="sfd"):
        """ """
        from .utils import get_mwebv
        mwebv = get_mwebv(*self.get_coordinates(), sourcemap=sourcemap)
        self.set_mwebv(mwebv)

    # ------- #
    # GETTER  #
    # ------- #
    #
    # - Extinction
    def get_mwebv(self, source=None):
        """ Get the Milky Way extinction in the target direction """
        if source is None or source in ["stored", "self"]:
            return self.mwebv
        from .utils import get_mwebv
        return get_mwebv(*self.get_coordinates(), sourcemap=source)
    
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
    def get_coordinates(self, full=False, perband=False, filters="*", clearna=True,
                            usestat='mean', **kwargs):
        """ shortcut to self.lightcurve.get_coordinates()
        (see also self.alerts.get_coordinates() or self.source.get_coordinates())
        
        Parameters
        ----------
        clearna: [bool] -optional-
            If the returned coordinates is not a simple Serie, shall this remove the NaN rows ?
        """
        if not self.has_lightcurve():
            warnings.warn("No lightcurve set, maybe see self.source.get_coordiantes()")
            return None
            
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
    def get_lcdata(self, detected=None, filters='*', time_range=None, query=None,  **kwargs):
        """ shortcut to  self.lightcurve.get_lcdata()  """
        if not self.has_lightcurve():
            raise AttributeError("no lightcurve set. See set_lightcurve().")
        
        return self.lightcurve.get_lcdata(detected=None, filters='*', time_range=None, query=None, **kwargs)
    
    # - Spectra
    def get_spectrum_phase(self, which=None):
        """ """
        if not self.has_spectra():
            raise AttributeError("No spectra loaded.")
        if not "phase" in self.specdata:
            raise AttributeError("No Phase loaded in specdata. Most likely fit_salt() did not run")
        
        if which is None:
            return np.asarray(self.specdata.phase, dtype="float")
        
        return self.specdata.iloc[which].phase
    
    def get_snidresult(self, which=None, safeout=True):
        """ """
        if not self.has_spectra():
            raise AttributeError("No spectra loaded.")
        if which is None:
            return self._call_down_spectra_("snidresult", isfunc=False)        

        if not safeout and not self.spectra[which].has_snidresult():
            raise AttributeError(f"fit_snid has not been ran for {which}")
                
        return self.spectra[which].snidresult
    
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
        if fix_redshift:
            force_redshift = self.get_redshift()
            
        output = self.lightcurve.fit_salt(incl_upperlimit=incl_upperlimit, filterprop=filterprop,
                                            force_color=force_color, force_x1=force_x1,
                                            force_t0=force_t0, force_redshift=force_redshift,
                                            fixed=fixed, values=values, bounds=bounds,
                                            get_object=get_object,
                                            set_it=set_it,
                                            **kwargs)
        if set_it:
            self._derive_saltresiduals_()

    def fit_snid(self, which="*", delta_phase=5, delta_redshift=None,
                     use_redshift=False, use_phase=True, squeeze=True, **kwargs):
        """ """
        if not self.has_spectra():
            raise AttributeError("No spectra loaded")

        if which in ["all","*"]:
            which = np.arange(self.nspectra)

        snidres = []
        for which_ in np.atleast_1d(which):        

            snid_prop = dict(delta_phase=delta_phase, delta_redshift=delta_redshift)
            if use_redshift:
                snid_prop["redshift"] = self.get_redshift()
            
            if use_phase:
                if not self.has_saltresult():
                    warnings.warn("No salt fitted, no phase available (which use_phase=True)")
                else:
                    snid_prop["phase"] = self.get_spectrum_phase(which_)

            snidres.append(self.spectra[which].fit_snid(**{**snid_prop,**kwargs}))

        return np.squeeze(snidres) if squeeze else snidres



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
    
    def show_lightcurve(self, **kwargs):
        """ """
        if not self.has_lightcurve():
            raise AttributeError("No lightcurve set. see self.set_lightcurve()")

        return self.lightcurve.show(**kwargs)
    
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
            label = f"{specdata_['instrument'].lower()} | {day}"
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

    def _build_specdata_(self, keys=["obj_id","instrument","mjd"]):
        """ """
        self._specdata = pandas.DataFrame([self._call_down_spectra_(k, isfunc=False)
                                               for k in keys],
                                            index=keys).T
        if self.has_saltresult():
            self._derive_saltresiduals_(["spectra"])
            
    def _derive_saltresiduals_(self, which=["lightcurve", "spectra"]):
        """ """
        if "lightcurve" in which and self.has_lightcurve():
            self.lightcurve._derive_saltresiduals_()
        if "spectra" in which and self.has_spectra():
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
        if self.has_lightcurve():
            return self.lightcurve.lcdata
        return None
    
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
        if self.has_lightcurve():
            return self.lightcurve.saltresult
        return None
    
    def has_saltresult(self):
        """ """
        return self.saltresult is not None
    
    # redshift
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
    
    # Extinction
    @property
    def mwebv(self):
        """ milky way extinction in the target direction """
        if not hasattr(self, "_mwebv") or self._mwebv is None:
            coord = self.get_coordinates()
            if coord is not None:
                self.load_mwebv()
            else:
                self._mwebv = None
                
        return self._mwebv
