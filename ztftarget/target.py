#! /usr/bin/env python
#

""" ZTF Target Lightcurve tools """


import numpy as np
import pandas

from . import lightcurve
from . import spectrum
from . import io

class ZTFTarget( object ):
    """ ZTF Target object """
    def __init__(self, name=None):
        """ Name of the ZTF target """
        self._name = name
        
    @classmethod
    def from_name(cls, name, lc_source="marshal", spectral_source="marshal", setup_source="marshal"):
        """ Load ZTF data given the name.
        
        Parameters
        ----------
        name: [string]
            Name of the ZTF target

        lc_source: [string] -optional-
            The source of the lightcurve
            available:
                - marshal

        spectral_source: [string] -optional-
            The source of the spectra
            available:
                - marshal

        setup_source: [string] -optional-
            Set the target information (redshift, coordinates and classification) from the given source
            available:
                - marshal
            
        Returns
        -------
        ZTFTarget
        """
        this = cls(name)
        if lc_source is not None:
            this.set_lightcurve( lightcurve.ZTFLightCurve.from_name(this.name, source=lc_source) )
            
        if spectral_source is not None:
            this.set_spectra( spectrum.ZTFSpectrum.from_name(this.name, source=spectral_source) )
            
        if setup_source is not None:
            this.setup(source=setup_source)
            
        return this
    
    # ----------------- #
    #   Methods         #
    # ----------------- #
    # -------- #
    #  LOADER  #
    # -------- #
    def setup(self, source=None, z=None, zerr=None, radec=None, classification=None):
        """ This sets the target information (redshift, coordinates, classification).
        They could come from a source

        Parameters
        ----------
        source: [string] -optional-
            Set the target information (redshift, coordinates and classification) from this source.
            = If not None, all the other parameters are ignored =
            available:
                - marshal

        // Only if source=None //
        z, zerr: [float] -optional-
            target redshift and its error. Any could be None.

        radec: [2d array] -optional-
            RA and Dec in degree (like radec=[210.7736, 54.2736])

        classification: [string] -optional-
            Classification type of the transient, like 'SNIa'
        
        Returns
        -------
        Void

        Example:
        --------
        self.setup(source='marshal')
        """
        if source is not None:
            if source in ["marshal"]:
                if not self.name in io.MARSHALQUERY.target_sources["name"].values:
                    warnings.warn(f"Cannot set target information from the marshal, {self.name} is unknown from your local data ; see io.MARSHALQUERY")
                    return None
                z, zerr = io.MARSHALQUERY.get_target_redshift(self.name).values[0], None
                radec = io.MARSHALQUERY.get_target_coordinates(self.name).values[0]
                classification = io.MARSHALQUERY.get_target_classification(self.name).values[0]
            else:
                raise NotImplementedError(f"source {source} has not been implemented.")
            
        self.set_redshift(z, zerr=zerr)
        self.set_coordinates( *radec )
        self.set_classification( classification )
        
    # -------- #
    #  SETTER  #
    # -------- #
    def set_redshift(self, z, zerr=None):
        """ Redshift of the target """
        self._z = z
        self._zerr = zerr

    def set_coordinates(self, ra, dec):
        """ RA and Dec of the target """
        self._ra, self._dec  = ra, dec

    def set_classification(self, classification):
        """ classification (type) of the target (e.g. SNIa) """
        self._classification = classification
        
    def set_lightcurve(self, lightcurve):
        """ Provide the target lightcurve 
        
        Parameters
        ----------
        lightcurve: [dataframe or ZTFLightCurve]
            If you provide a dataframe, it will be converted into a ZTFLightCurve.
            
        Returns
        -------
        Void
        """
        if type(lightcurve) == pandas.DataFrame:
            lightcurve = ZTFLightCurve(lightcurve, self.name)
        self._lc = lightcurve
        
    def set_spectra(self, spectra):
        """ """
        self._spectra = spectra
        
    # -------- #
    #  GETTER  #
    # -------- #
    # -------- #
    #  Other   #
    # -------- #
    def view_on_marshal(self):
        """ opens your browser at the corresponding target marshal page"""
        if self.name is None:
            raise AttributeError("self.name is not set. Cannot launch target Marshal page.")
        
        import webbrowser
        return webbrowser.open(f'http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name={self.name}', new=2)
        
    # -------- #
    #  PLOTTER #
    # -------- #
    def show(self):
        """ """
        import matplotlib.pyplot as mpl
        from astropy import time
        
        fig = mpl.figure(figsize=[9,3])
        axlc = fig.add_axes([0.1,0.175,0.35,0.7])
        axspec = fig.add_axes([0.525,0.175,0.4,0.7])
    
        self.lightcurve.show(ax=axlc)

        color_spec = "grey"
        self.spectra.show(ax=axspec, color=color_spec)
        axlc.axvline(time.Time(self.spectra.header.get("JD"), format="jd").datetime, 
                         ls="-", lw=1, color=color_spec, ymin=0.95,ymax=1, 
                         ms=3, marker="o", markevery=[True, False])
        axspec.set_ylabel("Flux []")
        axspec.set_xlabel(r"Wavelength [$\AA$]")
        fig.text(0.02,0.98, self.name, color=mpl.cm.binary(0.7),
                     va="top", ha="left")
        
    # ----------------- #
    #   Properties      #
    # ----------------- #    
    @property
    def name(self):
        """ ZTF target name """
        return self._name
    
    @property
    def radec(self):
        """ ZTF target RA, Dec """
        return self._ra, self._dec

    @property
    def z(self):
        """ redshift (and error) of the target """
        if not hasattr(self,"_z"):
            self._z = None
        return self._z

    @property
    def zerr(self):
        """ redshift (and error) of the target """
        if not hasattr(self,"_zerr"):
            self._zerr = None
        return self._zerr

    @property
    def classification(self):
        """ target classification """
        if not hasattr(self,"_classification"):
            self._classification = None
        return self._classification

    
    @property
    def lightcurve(self):
        """ ZTF target lightcurve """
        return self._lc

    def has_lightcurve():
        """ Test if the lightcurve has been set. """
        return hasattr(self,"_lc") and self._lc is not None
    
    @property
    def spectra(self):
        """ """
        return self._spectra

    def has_spectra():
        """ Test if the lightcurve has been set. """
        return hasattr(self,"_spectra") and self._spectra is not None
    
    def has_multiple_spectra():
        """ Test if the lightcurve has been set. """
        return self.has_spectra() and np.atleast_1d(self.spectra)>1
