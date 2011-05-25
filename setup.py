import sys

# Get version from phycasver.txt file
phycas_major, phycas_minor, phycas_bugfix = tuple(open('phycasver.txt','r').read().split('.'))

from distutils.core import setup, Extension
import distutils.sysconfig
import distutils
import os

# using MANIFEST now to specify all this
#data_all = [
#            'src/*',
#            'src/thirdparty/*',
#            'src/thirdparty/dcdflib/doc/*',
#            'src/thirdparty/dcdflib/src/*',
#            'src/thirdparty/num_util/*',
#            'Tests/Data/*.nex',
#            'Tests/ExplorePrior/*.py',
#            'Tests/ExplorePrior/reference_output/*',
#            'Tests/SplitTest/*.py',
#            'Tests/SplitTest/reference_output/*',
#            'Tests/PDFTree/*.py',
#            'Tests/PDFTree/reference_output/*',
#            'Tests/GTRTest/*.py',
#            'Tests/GTRTest/reference_output/*',
#            'Tests/SumP/*.py',
#            'Tests/SumP/*.p',
#            'Tests/SumP/reference_output/*',
#            'Tests/SumT/*.py',
#            'Tests/SumT/*.t',
#            'Tests/SumT/reference_output/*',
#            'Tests/LikelihoodTest/*.py',
#            'Tests/LikelihoodTest/reference_output/*',
#            'Tests/LikelihoodTest/acceptable_diff/check.nex',          
#            'Tests/FixedParams/*.py',
#            'Tests/FixedParams/reference_output/*',
#            'Tests/Underflow/*.py',
#           'Tests/Underflow/*.tre',
#            'Tests/Underflow/reference_output/*',
#            'Tests/CodonTest/*.py',
#            'Tests/CodonTest/reference_output/*',
#            'Examples/Paradox/*',
#            'Examples/Steppingstone/*',
#            'Examples/Simplest/*',
#            'Utilities/*.py',
#            'TreeViewer/*.py',
#            'PDFGen/*.py',
#            'PDFGen/AFM/*.afm',
#            ]

#windows_package_data = {
#            'phycas.Conversions': ['*.pyd','*.dll'],
#            'phycas.DataMatrix': ['*.pyd'],
#            'phycas.Likelihood': ['*.pyd'],
#            'phycas.ProbDist': ['*.pyd'],
#            'phycas.Phylogeny': ['*.pyd'],
#            'phycas.ReadNexus': ['*.pyd'],
#            }
#windows_package_data.update({'phycas': data_all})

isWin = sys.platform == 'win32'

sharedObjSuffix = [distutils.sysconfig.get_config_var("SO")]
dynamicLibSuffix = isWin and ['.dll'] or [".dylib"]
conversionsDataFiles = sharedObjSuffix + dynamicLibSuffix

#test_data = [
#          'Tests/cleanall.py',
#          'Tests/doctestall.py',
#          'Tests/runall.py',
#          ]
#_pack_data = {
#        'phycas': data_all,
#        }
phycas_description = """\
Phycas:
    
Phycas is a Python application for carrying out phylogenetic analyses.
It is also a C++ and Python library that can be used to create new
applications or to extend the current functionality.
"""

phycas_full_version = phycas_major+'.'+phycas_minor+'.'+phycas_bugfix

setupArgs = {
    'name':'Phycas',
    'version':phycas_full_version,
    'description':'Phycas: Python Software for Phylogenetic Analysis',
    'author':'Phycas Development Team',
    'author_email':'developers@phycas.org',
    'url':'http://www.phycas.org/',
    'license':'GNU General Public License (GPL)',
    'package_dir':{'': '.'},
    'packages':[
        'phycas',
        'phycas.Conversions',
        'phycas.DataMatrix',
        'phycas.Likelihood',
        'phycas.PDFGen',
        'phycas.ProbDist',
        'phycas.Phycas',
        'phycas.Phylogeny',
        'phycas.ReadNexus',
        'phycas.TreeViewer',
        'phycas.Utilities',
        ],
    'long_description':phycas_description,
    'platforms':['Linux', 'MacOS X', 'Windows'],
    'keywords':['phylogeny', 'phylogenetics', 'MCMC', 'Bayesian', 'bioinformatics'],
    'classifiers':[
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
     }

#if isWin:
#    setupArgs.update({
#        'package_data' : windows_package_data,
#        })
#else:
#    #setupArgs.update({
#    #    'package_data' : _pack_data,
#    #    })
#    parent_dir = os.path.split(sys.argv[0])[0]
#    boost_target_dir = 'boost_1_42_0'
#    include_dirs = [
#        "/usr/include",
#        "/usr/local/include/"+boost_target_dir, 
#        os.path.abspath(parent_dir),
#        os.path.abspath(os.path.join(parent_dir, "phycas", "src")),
#        ]
#    libraries=["boost_python"]
#    library_dirs=['/usr/local/lib']

setup(**setupArgs)


