import numpy as np


def get_dustmap(sourcemap, useweb=False):
    """ get the dustmap (from the dustmaps package) of the given source.

    Parameters
    ---------
    sourcemap: [string]
        origin of the MW extinction information. 
        currently implemented: planck, sfd
        
    useweb: [bool] -optional-
        shall this query from the web 
        = only implemented for sfd, ignored otherwise =
        
    Returns
    -------
    dustmaps.Dustmap
    """
    if sourcemap.lower() == "sfd":
        from dustmaps import sfd
        return sfd.SFDQuery() if not useweb else sfd.SFDWebQuery()
    
    if sourcemap.lower() == "planck":
        from dustmaps import planck
        return planck.PlanckQuery()
    
    raise NotImplementedError(f"Only Planck and SFD maps implemented. {sourcemap} given.")

def get_mwebv(ra, dec, sourcemap="sfd"):
    """ get the Milky Way E(B-V) at the given coordinates
    
    Parameters
    ----------
    ra, dec: [float or list of]
        RA and Dec (in deg) ; could be lists

    sourcemap: [string] -optional-
        origin of the MW extinction information. 
        (see get_dustmap) e.g.:
        - planck or sfd

    Returns
    -------
    float (or list following the ra, dec input)
    """    
    from astropy.coordinates import SkyCoord
    coords = SkyCoord(ra, dec, unit="deg")
    return get_dustmap(sourcemap)(coords)


    

def mag_to_flux(mag, magerr=None, units="zp", zp=25.0, wavelength=None):
    """converts magnitude into flux
    Parameters
    ----------
    mag: [float or array]
        AB magnitude(s)
    magerr: [float or array] -optional-
        magnitude error if any
    units: [string] -optional-
        Unit system in which to return the flux:
        - 'zp':   units base on zero point as required for sncosmo fits
        - 'phys': physical units [erg/s/cm^2/A)
    zp: [float or array] -optional-
        zero point of for flux; required if units == 'zp'
    wavelength: [float or array] -optional-
        central wavelength of the photometric filter.
        In Angstrom; required if units == 'phys'
    Returns
    -------
    - float or array (if magerr is None)
    - float or array, float or array (if magerr provided)
    """
    if units not in ["zp", "phys"]:
        raise ValueError("units must be 'zp' or 'phys'")
    elif units == "zp":
        if zp is None:
            raise ValueError("zp must be float or array if units == 'zp'")
        flux = 10 ** (-(mag - zp) / 2.5)
    else:
        if wavelength is None:
            raise ValueError("wavelength must be float or array if units == 'phys'")
        flux = 10 ** (-(mag + 2.406) / 2.5) / wavelength ** 2

    if magerr is None:
        return flux

    dflux = np.abs(flux * (-magerr / 2.5 * np.log(10)))  # df/f = dcount/count
    return flux, dflux


def flux_to_mag(flux, dflux, units="zp", zp=25.0, wavelength=None):
    """ Converts fluxes (erg/s/cm2/A) into AB magnitudes """
    if units not in ["zp", "phys"]:
        raise ValueError("units must be 'zp' or 'phys'")
    elif units == "zp":
        if zp is None:
            raise ValueError("zp must be float or array if units == 'zp'")
        wavelength = 1.0
    else:
        if wavelength is None:
            raise ValueError("wavelength must be float or array if units == 'phys'")
        zp = -2.406

    if dflux is None:
        return -2.5 * np.log10(flux * wavelength ** 2) + zp, None

    err = -2.5 / np.log(10) * dflux / flux

    return -2.5 * np.log10(flux * wavelength ** 2) + zp, np.abs(err)


def get_upper_limit_fluxes(dataframe, units="zp", zp=25.0, upmagkey="upmag"):
    """ """
    if units != "zp":
        raise NotImplementedError("only units=zp implemented.")
        
    error = (
        mag_to_flux(
            np.asarray(dataframe[upmagkey]),
            None,
            units=units,
            zp=zp,
        ) / 5
    )
    flux = np.random.normal(loc=0, scale=error)
    return flux, error


def get_fluxes(dataframe, units="zp", zp=25.0, magkey="mag", magerrkey="mag.err"):
    """ """
    if units != "zp":
        raise NotImplementedError("only units=zp implemented.")
        
    return mag_to_flux(dataframe[magkey], dataframe[magerrkey], units=units, zp=zp)        



def read_lines(lines):
    """ """
    import warnings
    DICTLINES = {"ha":6562, "hb":4861, "hg":4340, "hd":4101}
    linesout = []
    for l in np.atleast_1d(lines).tolist():
        if type(l) is str or type(l) in [np.str_,np.str]:
            if l in DICTLINES:
                linesout.append(DICTLINES[l])
            else:
                try:
                    linesout.append(float(l))
                except:
                    warnings.warn(f"cannot parse line {l} ; ignored")
                    continue
        else:
            linesout.append(float(l))

    return linesout
