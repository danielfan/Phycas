from setuptools import setup
#from distutils.core import setup

setup(name='phypy',
      version='1.0',
      description='Python Software for Phylogenetic Analysis',
      author='Paul O. Lewis',
      author_email='paul.lewis@uconn.edu',
      url='http://www.phypy.org/',
      packages=['pyphy','pyphy.Conversions','pyphy.DataMatrix','pyphy.Likelihood','pyphy.ProbDist','pyphy.Phylogeny','pyphy.ReadNexus'],
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
#    keywords = "pyphy phypy phycas Bayesian phylogenetics phylogeny likelihood",
#    url = "http://www.phypy.org/",   # project home page, if any
#
#    # could also include long_description, download_url, classifiers, etc.
#)
