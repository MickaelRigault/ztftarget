""" Spectroscopy module. Based on ztfquery.fritz and ztfquery.sedm """

import os
import numpy as np
import warnings

from ztfquery import fritz, sedm


class Spectrum(  fritz.FritzSpectrum ):

    @classmethod
    def from_name(cls, name, force_dl=False, warn=True, store=True, source=["fritz","sedm"], 
                  avoid_duplicate=True,  **kwargs):
        """ """
        sources = []
        speckey = []
        prop = {**dict(force_dl=force_dl, warn=warn), **kwargs}
        
        # - fritz sources        
        if "fritz" in source:
            fritzspec = cls.from_fritz(name, store=store, **prop)
            if fritzspec is not None:
                sources += [s for s in fritzspec if (not avoid_duplicate or s.filekey not in speckey)] 
                speckey = [s.filekey for s in sources]
            
        # - pharos sources
        if "sedm" in source:
            pharosspec = cls.from_pharos(name, **prop) # stored automatically
            if pharosspec is not None:
                sources += [s for s in np.atleast_1d(pharosspec) if (not avoid_duplicate or s.filekey not in speckey)] 
                speckey = [s.filekey for s in sources]

        return sources
    
    @classmethod
    def from_fritz(cls, name, warn=True, force_dl=False, store=False, **kwargs):
        """ """
        return super().from_name(name, warn=warn, force_dl=force_dl, store=store, **kwargs)
    
    @classmethod
    def from_pharos(cls, name, warn=True, force_dl=False, show_progress=False, **kwargs):
        """ """
        squery = sedm.SEDMQuery()
        specfiles = squery.get_target_spectra(name, extension=".fits", force_dl=force_dl,
                                                  show_progress=show_progress, **kwargs)
        
        if specfiles is None or len(specfiles)==0:
            warnings.warn(f"no spectra on pharos for {name}")
            return None
        
        if len(specfiles)>1 and warn:
            warnings.warn(f"{name} has several spectra, list of FritzSpectrum returned")
            
        specs = []
        for specfile_ in specfiles:
            spec = cls.read_fits(specfile_)
            spec._fritz_it_(warn=warn, instrument="SEDM")
            specs.append(spec)
            
        return spec

    # ============= #
    #  Methods      #
    # ============= #
    def set_snidresult(self, snidresult):
        """ """
        self._snidresult = snidresult

    def fit_snid(self, phase=None, redshift=None, delta_phase=5, delta_redshift=None,
                     lbda_range=[4000,8000], set_it=True,
                     verbose=False, quiet=True, **kwargs):
        """ """
        from . import snid
        snid_prop = dict(quiet=quiet, lbda_range=lbda_range, verbose=verbose)
        
        #
        # - Phase
        if phase is not None:
            snid_prop["phase_range"]=[phase-delta_phase, phase+delta_phase]
        #
        # - redshift            
        if redshift is not None:
            snid_prop["forcez"] = redshift
            if delta_redshift is not None:
                snid_prop["redshift_range"] = [redshift-delta_redshift, redshift+delta_redshift]
                

        # - Running SNID
        snidf = snid.SNID()
        outfile = snidf.run(self.filename, **{**snid_prop,**kwargs})
        if outfile is None:
            warnings.warn("SNID fit failed. Nothing returned")
            return None

        snidres = snid.SNIDReader.from_filename(outfile)
        if set_it:
            self.set_snidresult(snidres)
            
        return snidres

    # ============= #
    #  Internal     #
    # ============= #    
    def _fritz_it_(self, warn=True, instrument="unknown"):
        """ Internal method used when a non-fritz spectra is loaded. """
        from astropy import time
        # not using self.header to avoid recursivity
        try:
            self._header.loc["OBSERVED_AT"] = time.Time(self._header.loc["MJD_OBS"],
                                                            format="mjd").isot
            self._header.loc["OBJID"] = self._header.loc["OBJNAME"].replace('"','')
        except:
            if warn:
                warnings.warn("_fritzit_: Cannot set OBSERVED_AT and OBJ_ID from header")
            
        self._header.loc["INSTNAME"] = instrument
        if self.filename is not None:
            self._header.loc["OFNAME"] = os.path.basename(self.filename).split(".")[0]
        elif warn:
            warnings.warn("filename is None. cannot set OFNAME")


    # ============= #
    #  Properties   #
    # ============= #    
    @property
    def snidresult(self):
        """ """
        if not hasattr(self, "_snidresult"):
            return None
        return self._snidresult
    
    def has_snidresult(self):
        """ """
        return self.snidresult is not None
