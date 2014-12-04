#-------------------------------------------------------------------------------
# Name:		setup.py
# Purpose:      Installation of OpenWind
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	19.11.2014
# Last modified:04.12.2014 13:12
# Copyright:    (c) NERSC
# License:      GNU General Public License
#-------------------------------------------------------------------------------
import sys
import os

NAME                = 'openwind'
MAINTAINER          = "Nansat and OpenWind Developers"
MAINTAINER_EMAIL    = "nansat-dev@googlegroups.com"
DESCRIPTION         = "A python package for estimating high resolution wind from SAR images"
LONG_DESCRIPTION    = "A python package for estimating high resolution wind from SAR images"
MAJOR               = 0
MINOR               = 7
MICRO               = 0
# Make sure to tag releases as "v%d.%d.%d" %(MAJOR, MINOR, MICRO)
URL                 = "https://github.com/nansencenter/openwind"
DOWNLOAD_URL        = "https://github.com/nansencenter/openwind/archive/v%d.%d.%d.tar.gz" %(MAJOR, MINOR, MICRO)
LICENSE             = "GNU General Public License"
CLASSIFIERS         = '***'  # filter(None, CLASSIFIERS.split('\n'))
AUTHOR              = ("Morten W. Hansen, Knut-Frode Dagestad, Asuka Yamakawa")
AUTHOR_EMAIL        = "nansat-dev@googlegroups.com"
PLATFORMS           = ["UNKNOWN"]
ISRELEASED          = False
VERSION             = '%d.%d-dev.%d' % (MAJOR, MINOR, MICRO) 
REQS                = []

import_error_msg = "OpenWind v%d.%d.%d requires %s, which should be installed separately"

# Check if required packages are installed
try:
    import nansat
except ImportError:
    raise ImportError(import_error_msg %(MAJOR,MINOR,MICRO,'nansat'))

#----------------------------------------------------------------------------#
#                               Install package
#----------------------------------------------------------------------------#
from setuptools import setup, find_packages
from setuptools.command.install_scripts import install_scripts
from distutils import log
from distutils.errors import CCompilerError, DistutilsExecError,\
    DistutilsPlatformError

# Windows batch file handling
# from https://matthew-brett.github.io/pydagogue/installing_scripts.html
BAT_TEMPLATE = \
r"""@echo off
set mypath=%~dp0
set pyscript="%mypath%{FNAME}"
set /p line1=<%pyscript%
if "%line1:~0,2%" == "#!" (goto :goodstart)
echo First line of %pyscript% does not start with "#!"
exit /b 1
:goodstart
set py_exe=%line1:~2%
call %py_exe% %pyscript% %*
"""
class my_install_scripts(install_scripts):
    def run(self):
        install_scripts.run(self)
        if not os.name == "nt":
            return
        for filepath in self.get_outputs():
            # If we can find an executable name in the #! top line of the script
            # file, make .bat wrapper for script.
            with open(filepath, 'rt') as fobj:
                first_line = fobj.readline()
            if not (first_line.startswith('#!') and
                    'python' in first_line.lower()):
                log.info("No #!python executable found, skipping .bat "
                            "wrapper")
                continue
            pth, fname = os.path.split(filepath)
            froot, _ = os.path.splitext(fname)
            bat_file = os.path.join(pth, froot + '.bat')
            bat_contents = BAT_TEMPLATE.replace('{FNAME}', fname)
            log.info("Making %s wrapper for %s" % (bat_file, filepath))
            if self.dry_run:
                continue
            with open(bat_file, 'wt') as fobj:
                fobj.write(bat_contents)

# the following is adapted from simplejson's setup.py
if sys.platform == 'win32' and sys.version_info > (2, 6):
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    # It can also raise ValueError http://bugs.python.org/issue7511
    ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError,
                  IOError, ValueError)
else:
    ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)

def run_setup():
    kw = dict()

    setup(
        name=NAME,
        version=VERSION,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        platforms=PLATFORMS,
        packages=find_packages(),
        package_data={NAME:['tests/plots/*.png']},
        install_requires=REQS,
        **kw
        )

run_setup()

import shutil
shutil.rmtree('dist')
shutil.rmtree('build')
shutil.rmtree('openwind.egg-info')
