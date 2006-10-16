# To create a Windows installer...
#   python setup.py bdist_wininst --title="Phycas 1.0"
# This will create two directories: build and dist
# The installer will be placed in dist, build can be deleted

from distutils.core import setup

data_all = [
          'Tests/doctestall.py',
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
          'Examples/*.py'
          ]

data_windows_only = [
                  'Tests/cleanall.bat',
                  'Tests/doctestall.bat',
                  'Tests/runall.bat'
                  ]
data_windows_only.extend(data_all)

windows_package_data = {
                    'phypy.Conversions': ['*.pyd','*.dll'],
                    'phypy.DataMatrix': ['*.pyd'],
                    'phypy.Likelihood': ['*.pyd'],
                    'phypy.ProbDist': ['*.pyd'],
                    'phypy.Phylogeny': ['*.pyd'],
                    'phypy.ReadNexus': ['*.pyd']
                    }
windows_package_data.update({'phypy': data_windows_only})

phycas_description = """\
Phycas and the PhyPy library:
 	
Phycas is a Python application for carrying out phylogenetic analyses.
The PhyPy library is a C++ and Python library used by Phycas, but which
can be used to create new applications or to extend the functionality
currently built into Phycas.
"""

setup(name='Phycas',
      version='1.0',
      description='Phycas: Python Software for Phylogenetic Analysis',
      author='Phycas Development Team',
      author_email='phycas@phypy.org',
      url='http://www.phypy.org/',
      license='GNU General Public License (GPL)',
      scripts=['win_shortcuts.py'],
      package_dir = {'': 'phypy'},
      packages=['phypy',
                'phypy.Conversions',
                'phypy.DataMatrix',
                'phypy.Likelihood',
                'phypy.ProbDist',
                'phypy.Phycas',
                'phypy.Phylogeny',
                'phypy.ReadNexus'],
      package_data=windows_package_data,
      long_description=phycas_description,
      platforms=['Linux', 'MacOS X', 'Windows'],
      keywords=['phylogeny', 'phylogenetics', 'MCMC', 'Bayesian', 'bioinformatics'],
      classifiers=[
                  'Development Status :: 3 - Alpha',
                  'Environment :: Console',
                  'Environment :: Win32 (MS Windows)',
                  'Environment :: MacOS X',
                  'Intended Audience :: End Users/Desktop',
                  'Intended Audience :: Education',
                  'Intended Audience :: Science/Research',
                  'License :: OSI Approved :: GNU General Public License (GPL)',
                  'Natural Language :: English',
                  'Operating System :: MacOS :: MacOS X',
                  'Operating System :: Microsoft :: Windows',
                  'Operating System :: POSIX',
                  'Programming Language :: Python',
                  'Programming Language :: C++',
                  'Topic :: Scientific/Engineering :: Bio-Informatics'
                  ],      
     )
