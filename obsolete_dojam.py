#! /usr/bin/env python
import os,sys

#allow the script to be invoked from any dir, but use the script's parent as $PHYCAS_ROOT
scriptPar = os.path.split(os.path.abspath(sys.argv[0]))[0]
os.chdir(scriptPar)

env = os.environ
warnings = []
if "PYTHON_ROOT" not in env:
    sys.exit("""Please define PYTHON_ROOT environmental variable as path to Python24 directory
For example (on Mac):
    if using bash: export PYTHON_ROOT=/Library/Frameworks/Python.framework/Versions/2.4
    if using tcsh: setenv PYTHON_ROOT /Library/Frameworks/Python.framework/Versions/2.4
""")
print "PYTHON_ROOT =", env["PYTHON_ROOT"]

if "PYTHON_VERSION" not in env:
        sys.exit("""Please define PYTHON_VERSION environmental variable as Python major version
For example:
    if using bash: export PYTHON_VERSION=2.4
    if using tcsh: setenv PYTHON_VERSION 2.4
""")
print "PYTHON_VERSION =", env["PYTHON_VERSION"]

if False:
    if "BOOST_ROOT" not in env:
        sys.exit("""Please define BOOST_ROOT environmental variable as path to boost directory
    For example:
      if using bash: export BOOST_ROOT=$HOME/boost_1_33_1
      if using tcsh: setenv BOOST_ROOT $HOME/boost_1_33_1
    """)
    print "BOOST_ROOT =", env["BOOST_ROOT"]

if "PHYCAS_ROOT" not in env:
    env["PHYCAS_ROOT"] = os.path.split(scriptPar)[0]
    s = '"PHYCAS_ROOT" not in enviroment. The default: "%s"\n is being used.' % env["PHYCAS_ROOT"]
    warnings.append(s)
else:
    print "PHYCAS_ROOT =", env["PHYCAS_ROOT"]

if "BUILD" not in env:
    env["BUILD"] = "debug"
    s = '"BUILD" not in enviroment, tho default is "debug" (set to "BUILD=release" before invoking this script to override this)'
    warnings.append(s)
else:
    print "BUILD =", env["BUILD"]

if "TOOLS" not in env:
    if sys.platform == "darwin":
        env["TOOLS"] = "darwin"
    elif sys.platform == "win32":
        env["TOOLS"] = "vc-7_1"
    else:
        env["TOOLS"] = "gcc"
    s = '"TOOLS" not in enviroment, the default for your platform is "%s"\n(see the Boost Build system documentation for the list of recognized values)' % env["TOOLS"]
    warnings.append(s)
else:
    print "TOOLS =", env["TOOLS"]

for w in warnings:
    print "Warning: ", w

if not os.system('bjam -q') == 0:
    sys.exit('Build failed')

for w in warnings:
    print "Warning: ", w
