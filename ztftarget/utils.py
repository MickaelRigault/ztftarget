import numpy as np


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


def get_upper_limit_fluxes(dataframe, units="zp", zp=25.0):
    """ """
    if units != "zp":
        raise NotImplementedError("only units=zp implemented.")
        
    error = (
        mag_to_flux(
            np.asarray(dataframe["upmag"]),
            None,
            units=units,
            zp=zp,
        ) / 5
    )
    flux = np.random.normal(loc=0, scale=error)
    return flux, error


def get_fluxes(dataframe, units="zp", zp=25.0):
    """ """
    if units != "zp":
        raise NotImplementedError("only units=zp implemented.")
        
    return mag_to_flux(dataframe["mag"], dataframe["mag.err"], units=units, zp=zp)        
