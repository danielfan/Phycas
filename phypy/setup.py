# To create a Windows installer...
#   python setup.py bdist_wininst --title="Phycas 1.0"
# This will create two directories: build and dist
# The installer will be placed in dist, build can be deleted

from distutils.core import setup

setup(name='phypy',
      version='1.0',
      description='Phycas: Python Software for Phylogenetic Analysis',
      author='Mark T. Holder, Paul O. Lewis and David L. Swofford',
      author_email='phycas@phypy.org',
      url='http://www.phypy.org/',
      packages=['phypy',
                'phypy.Conversions',
                'phypy.DataMatrix',
                'phypy.Likelihood',
                'phypy.ProbDist',
                'phypy.Phycas',
                'phypy.Phylogeny',
                'phypy.ReadNexus'],
      package_data={'phypy.Conversions': ['*.pyd'],
                    'phypy.DataMatrix': ['*.pyd'],
                    'phypy.Likelihood': ['*.pyd'],
                    'phypy.ProbDist': ['*.pyd'],
                    'phypy.Phylogeny': ['*.pyd'],
                    'phypy.ReadNexus': ['*.pyd'],
                    'phypy': ['Tests/doctestall.py',
                              'Tests/cleanall.bat',
                              'Tests/doctestall.bat',
                              'Tests/runall.bat',
                              'Examples/Data/green.nex',
                              'Examples/Data/nyldna4.nex',
                              'Examples/ExplorePrior/__init__.py',
                              'Examples/ExplorePrior/ExplorePrior.py',
                              'Examples/ExplorePrior/reference_output/nodata.nex.p',
                              'Examples/ExplorePrior/reference_output/nodata.nex.t',
                              'Examples/FixedParams/__init__.py',
                              'Examples/FixedParams/FixedParams.py',
                              'Examples/FixedParams/reference_output/params.p',
                              'Examples/FixedParams/reference_output/trees.t',
                              'Examples/GelfandGhosh/__init__.py',
                              'Examples/GelfandGhosh/GelfandGhosh.py',
                              'Examples/GelfandGhosh/reference_output/analHKY.nex.p',
                              'Examples/GelfandGhosh/reference_output/analHKY.nex.t',
                              'Examples/GelfandGhosh/reference_output/analHKYflex.nex.p',
                              'Examples/GelfandGhosh/reference_output/analHKYflex.nex.t',
                              'Examples/GelfandGhosh/reference_output/analHKYg.nex.p',
                              'Examples/GelfandGhosh/reference_output/analHKYg.nex.t',
                              'Examples/GelfandGhosh/reference_output/ggout.txt',
                              'Examples/LikelihoodTest/__init__.py',
                              'Examples/LikelihoodTest/LikelihoodTest.py',
                              'Examples/LikelihoodTest/reference_output/check.nex',
                              'Examples/LikelihoodTest/reference_output/simulated.nex',
                              'Examples/Polytomies/__init__.py',
                              'Examples/Polytomies/Polytomies.py',
                              'Examples/Polytomies/reference_output/analHKY.nex.p',
                              'Examples/Polytomies/reference_output/analHKY.nex.t',
                              'Examples/Polytomies/reference_output/simHKY.nex',
                              'Examples/Simulator/__init__.py',
                              'Examples/Simulator/Simulator.py',
                              'Examples/Simulator/reference_output/simulated.nex',
                              'Examples/*.py',
                              ],
                    }
     )

#from setuptools import setup, find_packages
#setup(
#    name = "phypy",
#    version = "1.0",
#    packages = find_packages(),
#    #scripts = ['say_hello.py'],
#
#    # Project uses reStructuredText, so ensure that the docutils get
#    # installed or upgraded on the target machine
#    #install_requires = ['docutils>=0.3'],
#
#    #package_data = {
#    #    # If any package contains *.txt or *.rst files, include them:
#    #    '': ['*.txt', '*.rst'],
#    #    # And include any *.msg files found in the 'hello' package, too:
#    #    'hello': ['*.msg'],
#    #}
#
#    # metadata for upload to PyPI
#    author = "Mark T. Holder, Paul O. Lewis, David L. Swofford",
#    author_email = "phycoids@phycas.org",
#    description = "PhyPy: Python Software for Phylogenetic Analysis",
#    license = "PSF",
#    keywords = "phypy phycas Bayesian phylogenetics phylogeny likelihood",
#    url = "http://www.phypy.org/",   # project home page, if any
#
#    # could also include long_description, download_url, classifiers, etc.
#)
