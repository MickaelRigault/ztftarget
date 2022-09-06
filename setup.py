#! /usr/bin/env python
#

DESCRIPTION = "ztftarget: ZTF tools assocated to target data"
LONG_DESCRIPTION = """ ztftarget: ZTF tools assocated to target data """

DISTNAME = "ztftarget"
AUTHOR = "Mickael Rigault"
MAINTAINER = "Mickael Rigault"
MAINTAINER_EMAIL = "m.rigault@ipnl.in2p3.fr"
URL = "https://github.com/MickaelRigault/ztftarget/"
LICENSE = "Apache 2.0"
DOWNLOAD_URL = "https://github.com/MickaelRigault/ztftarget/tarball/0.4"
VERSION = "0.4.0"

try:
    from setuptools import setup, find_packages

    _has_setuptools = True
except ImportError:
    from distutils.core import setup

    _has_setuptools = False


def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import ztfquery
    except ImportError:
        install_requires.append("ztfquery")

    try:
        import pandas
    except ImportError:
        install_requires.append("pandas")

    return install_requires


if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = ["ztftarget"]

    setup(
        name=DISTNAME,
        author=AUTHOR,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        scripts=[],
        packages=packages,
        include_package_data=True,
        package_data={'ztftarget': ['data/*','data/filters/*' ]},
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "License :: OSI Approved :: Apache License",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
    )
