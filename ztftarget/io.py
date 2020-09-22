#! /usr/bin/env python
#

""" Managing the IO for the module """


def load_marshalquery(program="*"):
    """ """
    from ztfquery import marshal
    global MARSHALQUERY
    MARSHALQUERY = marshal.MarshalAccess.load_local("*")
    
    
def load_ztf_obslogs(start="2018-03-01", end=None, **kwargs):
    """ """
    from ztfquery import skyvision
    global ZTF_OBSLOGS
    ZTF_OBSLOGS = skyvision.CompletedLog.from_daterange(start=start, end=end, **kwargs )
    
