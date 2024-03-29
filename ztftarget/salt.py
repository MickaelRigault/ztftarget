import numpy as np
import sncosmo
import pandas

from ztfquery import filters
filters.load_p48_filters_to_sncosmo(basename="ztf")


def get_saltmodel( mwebv=None):
    """ """
    import sncosmo
    dust = sncosmo.F99Dust()
    model =  sncosmo.Model("salt2", effects=[dust],
                       effect_names=['mw'],
                       effect_frames=['obs'])
    if mwebv is not None:
        model.set(mwebv=mwebv)
        
    return model

def sncosmoresult_to_pandas(result):
    """ """
    error = pandas.Series( dict(result.get("errors") ), name="error")
    values = pandas.Series( result.get("parameters"),
                                index=result.get("param_names"),
                                name="value")

    cov = pandas.DataFrame( result.get("covariance"), 
                                index=error.index, columns="cov_"+error.index)

    fit_res = pandas.concat( [values,error, cov], axis=1)
    fit_meta = pandas.Series( {k:result.get(k) for k in ["success", "ncall", "chisq", "ndof"]} )
    fit_meta["chi2dof"] = fit_meta["chisq"]/fit_meta["ndof"]
    return fit_res, fit_meta

class SALTResult( object ):

    def __init__(self, result=None, model=None, data=None):
        """ """
        if result is not None:
            self.set_result(result)

        if model is not None:
            self.set_model(model)

        if data is not None:
            self.set_data(data)
            
    @classmethod
    def from_fit(cls, table, values={}, fixed={}, bounds=None):
        """"""
        model = get_saltmodel()
        model.set(**{**values,**fixed})
        parameters = [p_ for p_ in ['z','t0', 'x0', 'x1', 'c'] if p_ not in fixed]
        result, fitted_model = sncosmo.fit_lc(table, model, 
                                              parameters,  # parameters of model to vary
                                              bounds=bounds)
        return cls(result=result, fitted_model=fitted_model, data=table)
    
    # ============= #
    #  Methods      #
    # ============= #
    # ------- #
    # SETTER  # 
    # ------- #
    def set_result(self, result):
        """ """
        self._result = result
        
    def set_model(self, model):
        """ """
        self._model = model
        
    def set_data(self, data):
        """ """
        from astropy.table import Table
        if type(data) == Table:
            data = data.to_pandas()
            
        self._data = data
        
    # ------- #
    # GETTER  # 
    # ------- #
    def to_pandas(self):
        """ """
        return sncosmoresult_to_pandas(self.sncosmo_result)
    
    def get_model(self):
        """ Set the model to the parameters and returns it. """
        self.model.set(**self.get_parameters(inclerrors=False, as_serie=False))
        return self.model
    
    def get_phase(self, mjd):
        """ """
        return mjd - self.get_parameters()["t0"]
    
    def get_parameters(self, inclerrors=True, as_serie=True):
        """ """
        param_ = {k:v for k,v in zip(self.sncosmo_result.param_names, self.sncosmo_result.parameters)}
        if inclerrors:
            for k_,v_ in self.sncosmo_result.errors.items():
                param_[f"{k_}_err"] = v_
                
        if as_serie:
            return pandas.Series(param_)
        return param_

    def get_lightcurve(self, bands, jd=None, timerange=[-20,50], bins=70,
                             zp=25, zpsys="ab", squeeze=False, as_phase=False,
                             influx=True, as_dataframe=True):
        """ returns the model lightcurve. 
        
        Parameters
        ----------
        bands: [string (or list of)]
            Name of the band for which you want a lightcurve model.
        """
        model = self.get_model()
        t0 = model.get("t0")
        if jd is None:
            jd = np.linspace(t0+timerange[0], t0+timerange[1], bins)
            
            
        if influx:
            modelfunc = model.bandflux
            modelprop = dict(zp=zp, zpsys=zpsys)
        else:
            modelfunc = model.bandmag
            modelprop = dict(magsys=zpsys)
            
        time_ =  jd if not as_phase else jd-t0
        if type(bands) == str and squeeze:
            bandflux = modelfunc(band=bands, time=jd, **modelprop)
            if not as_dataframe:
                return time_,bandflux
            return pandas.DataFrame({bands:bandflux,"time":time_})
        
        bandflux = {b_:modelfunc(band=b_, time=jd, **modelprop) for b_ in np.atleast_1d(bands)}
        if not as_dataframe:
            return time_, bandflux
        bandflux["time"] = time_
        return pandas.DataFrame(bandflux)
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def result(self):
        """ fitted_result and fitted_meta """
        return self.to_pandas()[0]

    @property
    def result_meta(self):
        """ """
        return self.to_pandas()[1]
    
    @property
    def sncosmo_result(self):
        """ """
        if not hasattr(self, "_result"):
            return None
        return self._result
    
    @property
    def model(self):
        """ """
        if not hasattr(self, "_model"):
            return None
        return self._model
    
    @property
    def data(self):
        """ """
        if not hasattr(self, "_data"):
            return None
        return self._data
