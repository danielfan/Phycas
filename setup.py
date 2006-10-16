# To create a Windows installer...
#   python setup.py bdist_wininst --title="Phycas 1.0"
# This will create two directories: build and dist
# The installer will be placed in dist, build can be deleted

from distutils.core import setup

data_all = ['phypy/Tests/doctestall.py',
          'phypy/Examples/Data/green.nex',
          'phypy/Examples/Data/nyldna4.nex',
          'phypy/Examples/ExplorePrior/__init__.py',
          'phypy/Examples/ExplorePrior/ExplorePrior.py',
          'phypy/Examples/ExplorePrior/reference_output/nodata.nex.p',
          'phypy/Examples/ExplorePrior/reference_output/nodata.nex.t',
          'phypy/Examples/FixedParams/__init__.py',
          'phypy/Examples/FixedParams/FixedParams.py',
          'phypy/Examples/FixedParams/reference_output/params.p',
          'phypy/Examples/FixedParams/reference_output/trees.t',
          'phypy/Examples/GelfandGhosh/__init__.py',
          'phypy/Examples/GelfandGhosh/GelfandGhosh.py',
          'phypy/Examples/GelfandGhosh/reference_output/analHKY.nex.p',
          'phypy/Examples/GelfandGhosh/reference_output/analHKY.nex.t',
          'phypy/Examples/GelfandGhosh/reference_output/analHKYflex.nex.p',
          'phypy/Examples/GelfandGhosh/reference_output/analHKYflex.nex.t',
          'phypy/Examples/GelfandGhosh/reference_output/analHKYg.nex.p',
          'phypy/Examples/GelfandGhosh/reference_output/analHKYg.nex.t',
          'phypy/Examples/GelfandGhosh/reference_output/ggout.txt',
          'phypy/Examples/LikelihoodTest/__init__.py',
          'phypy/Examples/LikelihoodTest/LikelihoodTest.py',
          'phypy/Examples/LikelihoodTest/reference_output/check.nex',
          'phypy/Examples/LikelihoodTest/reference_output/simulated.nex',
          'phypy/Examples/Polytomies/__init__.py',
          'phypy/Examples/Polytomies/Polytomies.py',
          'phypy/Examples/Polytomies/reference_output/analHKY.nex.p',
          'phypy/Examples/Polytomies/reference_output/analHKY.nex.t',
          'phypy/Examples/Polytomies/reference_output/simHKY.nex',
          'phypy/Examples/Simulator/__init__.py',
          'phypy/Examples/Simulator/Simulator.py',
          'phypy/Examples/Simulator/reference_output/simulated.nex',
          'phypy/Examples/*.py'
          ]

data_windows_only = ['phypy/Tests/cleanall.bat',
                  'phypy/Tests/doctestall.bat',
                  'phypy/Tests/runall.bat'
                  ]
data_windows_only.extend(data_all)

windows_package_data = {'phypy.Conversions': ['*.pyd','*.dll'],
                    'phypy.DataMatrix': ['*.pyd'],
                    'phypy.Likelihood': ['*.pyd'],
                    'phypy.ProbDist': ['*.pyd'],
                    'phypy.Phylogeny': ['*.pyd'],
                    'phypy.ReadNexus': ['*.pyd']}
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
      package_dir = {'': 'phypy'}
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
