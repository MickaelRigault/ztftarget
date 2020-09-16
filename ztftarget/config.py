#! /usr/bin/env python
#

""" Configuration options """


import matplotlib.pyplot as mpl

ZTFCOLOR = { # ZTF
        "ztf:r":dict(marker="o",ms=7,  mfc="C3"),
        "ztf:g":dict(marker="o",ms=7,  mfc="C2"),
        "ztf:i":dict(marker="o",ms=7, mfc="C1"),
        # Swift
        "uvot:B":   dict(marker="s",  ms=5, mfc="C0"),
        "uvot:u":   dict(marker="s",  ms=5, mfc=mpl.cm.Blues(0.7)),
        "uvot:uvm2":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.6)),
        "uvot:uvm2":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.8)),
        "uvot:uvm1":dict(marker="s", ms=5, mfc=mpl.cm.Purples(0.4)),
        "uvot:V":   dict(marker="s", ms=5, mfc=mpl.cm.Greens(0.9)),
        # 
        "ioo:u":   dict(marker="d", ms=6,mfc=mpl.cm.Blues(0.6)),
        "ioo:g":   dict(marker="d", ms=6,mfc=mpl.cm.Greens(0.6)),
        "ioo:r":   dict(marker="d", ms=6,mfc=mpl.cm.Reds(0.7)),
        "ioo:i":   dict(marker="d",ms=6, mfc=mpl.cm.Oranges(0.6)),
        "ioo:z":   dict(marker="d", ms=6,mfc=mpl.cm.binary(0.8))
        }
