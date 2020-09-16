#! /usr/bin/env python
#

""" ZTF Spectra """

import numpy as np
from pyifu import spectroscopy

def read_sedmspec(ascii_data):
    """ Returns header and data """
    from pyifu import get_spectrum
    from astropy.io import fits
    header = fits.Header()
    for d in ascii_data:
        if not d.startswith("#"): continue
        if "SNID" in d: continue
        header.set(*d.replace("# ","").split(": "))

    spec_data = np.asarray([d.split() for d in ascii_data if not d.startswith("#")], dtype="float").T
    return {"header":header, "lbda":spec_data[0], "flux":spec_data[1], "variance":spec_data[2]}



class ZTFSpectrum( spectroscopy.Spectrum ):
    
    def __init__(self, *args, name=None, **kwargs):
        """ """
        super().__init__(*args, **kwargs)
        self._name = name
        
    @classmethod
    def from_name(cls, name, source="marshal", only_sedm=True):
        """ """
        if source in ["marshal"]:
            from ztfquery import marshal
            specdata = [read_sedmspec(spec_) for spec_ in marshal.get_target_spectra(name, only_sedm=only_sedm).values()]
            # - Single            
            if len(specdata)==1:
                this = ZTFSpectrum(None, name=name)
                this.create(data=specdata[0]["flux"], variance=specdata[0]["variance"], 
                            header=specdata[0]["header"], lbda=specdata[0]["lbda"], 
                            logwave=None)
                return this
            # - List of
            these = []
            for i in range( len(specdata) ):
                this = ZTFSpectrum(None, name=name)
                this.create(data=specdata[i]["flux"], variance=specdata[i]["variance"], 
                            header=specdata[i]["header"], lbda=specdata[i]["lbda"], 
                            logwave=None)
                these.append(this)
            return these
                
        else:
            raise NotImplementedError(f"source '{source}' not implemented")
            

    # ----------------- #
    #   Properties      #
    # ----------------- #    
    @property
    def name(self):
        """ target name """
        return self._name
