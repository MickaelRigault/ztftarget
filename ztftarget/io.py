#! /usr/bin/env python
#

""" Managing the IO for the module """

import warnings
from ztfquery import marshal

MARSHALQUERY = marshal.MarshalAccess.load_local("*")

def update_marshal_sources(program="*"):
    """ update the sources and store them locally ; see marshal.MarshalAccess.store() 
    = This function is slow. =
    """
    MARSHALQUERY.load_target_sources(program)
    MARSHALQUERY.store()
    
