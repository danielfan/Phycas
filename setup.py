# To create a Windows installer...
#   python setup.py bdist_wininst --title="Phycas 1.0"
# This will create two directories: build and dist
# The installer will be placed in dist, build can be deleted

import sys

fakeCompilingExtensions = False # deprecated
compilingExt = sys.platform == "darwin"

if fakeCompilingExtensions:
    import distutils.ccompiler
    distutils.ccompiler.compiler_class['fake'] = ('ccompiler', 'FakeCompiler', "Dummy for a compiler that clones already built sources")
    import fakecompiler 
    distutils.ccompiler.FakeCompiler = fakecompiler.FakeCompiler
                               
from distutils.core import setup, Extension
import distutils.sysconfig
import distutils
import os

data_all = [
          'Tests/cleanall.py',
          'Tests/doctestall.py',
          'Tests/runall.py',
          'Tests/Data/green.nex',
          'Tests/Data/nyldna4.nex',
          'Tests/ExplorePrior/__init__.py',
          'Tests/ExplorePrior/ExplorePrior.py',
          'Tests/ExplorePrior/reference_output/nodata.nex.p',
          'Tests/ExplorePrior/reference_output/nodata.nex.t',
          'Tests/FixedParams/__init__.py',
          'Tests/FixedParams/FixedParams.py',
          'Tests/FixedParams/reference_output/params.p',
          'Tests/FixedParams/reference_output/trees.t',
          'Tests/GelfandGhosh/__init__.py',
          'Tests/GelfandGhosh/GelfandGhosh.py',
          'Tests/GelfandGhosh/reference_output/analHKY.nex.p',
          'Tests/GelfandGhosh/reference_output/analHKY.nex.t',
          'Tests/GelfandGhosh/reference_output/analHKYflex.nex.p',
          'Tests/GelfandGhosh/reference_output/analHKYflex.nex.t',
          'Tests/GelfandGhosh/reference_output/analHKYg.nex.p',
          'Tests/GelfandGhosh/reference_output/analHKYg.nex.t',
          'Tests/GelfandGhosh/reference_output/ggout.txt',
          'Tests/LikelihoodTest/__init__.py',
          'Tests/LikelihoodTest/LikelihoodTest.py',
          'Tests/LikelihoodTest/reference_output/check.nex',
          'Tests/LikelihoodTest/reference_output/simulated.nex',
          'Tests/Polytomies/__init__.py',
          'Tests/Polytomies/Polytomies.py',
          'Tests/Polytomies/reference_output/analHKY.nex.p',
          'Tests/Polytomies/reference_output/analHKY.nex.t',
          'Tests/Polytomies/reference_output/simHKY.nex',
          'Tests/Simulator/__init__.py',
          'Tests/Simulator/Simulator.py',
          'Tests/Simulator/reference_output/simulated.nex',
          'Tests/*.py',
          'Examples/Paradox/Paradox.py',
          'Examples/Paradox/ShoupLewis.nex'
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

isWin = sys.platform == 'win32'

sharedObjSuffix = [distutils.sysconfig.get_config_var("SO")]
dynamicLibSuffix = isWin and ['.dll'] or [".dylib"]
conversionsDataFiles = sharedObjSuffix + dynamicLibSuffix

test_data = [
          'Tests/cleanall.py',
          'Tests/doctestall.py',
          'Tests/runall.py'
          ]
_pack_data = {
        'phypy': data_all,
        }
phycas_description = """\
Phycas and the PhyPy library:
 	
Phycas is a Python application for carrying out phylogenetic analyses.
The PhyPy library is a C++ and Python library used by Phycas, but which
can be used to create new applications or to extend the functionality
currently built into Phycas.
"""

setupArgs = {
    'name':'Phycas',
    'version':'0.11.0',
    'description':'Phycas: Python Software for Phylogenetic Analysis',
    'author':'Phycas Development Team',
    'author_email':'phycas@phypy.org',
    'url':'http://www.phypy.org/',
    'license':'GNU General Public License (GPL)',
        'package_dir':{'': 'phypy'},
        'packages':[
         '',
         'phypy',
         'phypy.Conversions',
         'phypy.DataMatrix',
         'phypy.Likelihood',
         'phypy.ProbDist',
         'phypy.Phycas',
         'phypy.Phylogeny',
         'phypy.ReadNexus',
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

if isWin:
    setupArgs.update({
        'package_data' : windows_package_data,
        'scripts':['win_shortcuts.py'],
        })
else:
    setupArgs.update({
        'package_data' : _pack_data,
        })
    parent_dir = os.path.split(sys.argv[0])[0]
    include_dirs = [
        "/usr/include",
        "/usr/local/include/boost-1_33_1", 
        os.path.abspath(parent_dir),
        os.path.abspath(os.path.join(parent_dir, "phypy", "src")),
        ]
    libraries=["boost_python"]
    library_dirs=['/usr/local/lib']

if compilingExt:
    if fakeCompilingExtensions:
        em = [Extension('Conversions._Conversions', ['phypy/src/conversions_pymod.cpp'], 
                    
                    extra_link_args = [
                        #'--colocate-lib=libboost_python.dylib',
                        '--built-under=phypy',
                        '--path-from-package=phypy/Conversions', ]),
            Extension('DataMatrix._DataMatrixBase', ['phypy/src/data_matrix_pymod.cpp'], 
                    extra_link_args = [
                        '--built-under=phypy',
                        '--path-from-package=phypy/DataMatrix', ]),
            Extension('Likelihood._LikelihoodBase', ['phypy/src/likelihood_pymod.cpp'], 
                    extra_link_args = [
                        '--built-under=phypy',
                        '--path-from-package=phypy/Likelihood', ]),
            Extension('Phylogeny._Phylogeny', ['phypy/src/phylogeny_pymod.cpp'], 
                    extra_link_args = [
                        '--built-under=phypy',
                        '--path-from-package=phypy/Phylogeny', ]),
            Extension('ProbDist._ProbDist', ['phypy/src/probdist_pymod.cpp'], 
                    extra_link_args = [
                        '--built-under=phypy',
                        '--path-from-package=phypy/ProbDist', ]),
            Extension('ReadNexus._ReadNexus', ['phypy/src/read_nexus_pymod.cpp'], 
                    extra_link_args = [
                        '--built-under=phypy',
                        '--path-from-package=phypy/ReadNexus', ]),
            ]
    else:
        import subprocess
        arch = subprocess.Popen(["arch"], stdout=subprocess.PIPE).communicate()[0].strip()
        argDict = {
                "library_dirs": library_dirs,
                "libraries": libraries,
                "include_dirs": include_dirs,
                "depends": [],
                "extra_compile_args": [
                        "-include", 
                        "phypy/src/phypy_config.h", 
                        "-arch", 
                        arch,
                        "-ftemplate-depth-256",
                        "-finline-functions",
                        "-Wno-inline",
                        "-fPIC",
                        "-Wall",
                        "-Wno-four-char-constants",
                        "-Wno-unknown-pragmas",
                        "-F/Library/Frameworks",
                        ],
                "extra_link_args": [
                        "-arch", 
                        arch,
                        "-Wl,-x",
                        "-multiply_defined",
                        "suppress",
                        "-twolevel_namespace",
                        "-F/Library/Frameworks",
                        "-framework",
                        "Python",
                        ],
                "define_macros":[
                        ('NDEBUG', '1'),
                        ("BOOST_PYTHON_DYNAMIC_LIB", "1"),
                        ],
                 }
        
        em = [
              Extension('Conversions._Conversions', 
                    [
                    "phypy/src/conversions_pymod.cpp",
                    "phypy/src/boost_assertion_failed.cpp",
                    "phypy/src/basic_tree.cpp",
                    "phypy/src/basic_tree_node.cpp",
                    "phypy/src/phypy_string.cpp"
                    ], **argDict),
              Extension("DataMatrix._DataMatrixBase", 
                    [
                    "phypy/src/data_matrix_pymod.cpp",
                    "phypy/src/boost_assertion_failed.cpp",
                    "phypy/src/cipres/CipresDataMatrixHelper.cpp",
                    ],**argDict),
              Extension('Likelihood._LikelihoodBase', 
                    [
                    'phypy/src/basic_lot.cpp',
                    'phypy/src/basic_tree.cpp',
                    'phypy/src/basic_tree_node.cpp',
                    'phypy/src/boost_assertion_failed.cpp',
                    'phypy/src/bush_move.cpp',
                    'phypy/src/cipres/CipresDataMatrixHelper.cpp',
                    'phypy/src/cond_likelihood_storage.cpp',
                    'phypy/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phypy/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phypy/src/edge_move.cpp',
                    'phypy/src/internal_data.cpp',
                    'phypy/src/larget_simon_move.cpp',
                    'phypy/src/likelihood_loops.cpp',
                    'phypy/src/likelihood_models.cpp',
                    'phypy/src/likelihood_pymod.cpp',
                    'phypy/src/model_pymod.cpp',
                    'phypy/src/linalg.cpp',
                    'phypy/src/mcmc_chain_manager.cpp',
                    'phypy/src/mcmc_param.cpp',
                    'phypy/src/mcmc_flexcat_param.cpp',
                    'phypy/src/mcmc_updater.cpp',
                    'phypy/src/ncat_move.cpp',
                    'phypy/src/phypy_string.cpp',
                    'phypy/src/q_matrix.cpp',
                    'phypy/src/sim_data.cpp',
                    'phypy/src/slice_sampler.cpp',
                    'phypy/src/tip_data.cpp',
                    'phypy/src/topo_prior_calculator.cpp',
                    'phypy/src/tree_likelihood.cpp',
                    'phypy/src/tree_manip.cpp',
                    'phypy/src/underflow_policy.cpp',
                    'phypy/src/updater_pymod.cpp',
                    ], **argDict),
              Extension('Phylogeny._Phylogeny', 
                    [
                    'phypy/src/phylogeny_pymod.cpp',
                    'phypy/src/boost_assertion_failed.cpp',
                    'phypy/src/basic_tree.cpp',
                    'phypy/src/basic_tree_node.cpp',
                    'phypy/src/tree_manip.cpp',
                    'phypy/src/basic_lot.cpp',
                    'phypy/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phypy/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phypy/src/basic_cdf.cpp',
                    'phypy/src/phypy_string.cpp', 
                    ],**argDict),
              Extension('ProbDist._ProbDist', 
                    [
                    'phypy/src/probdist_pymod.cpp',
                    'phypy/src/boost_assertion_failed.cpp',
                    'phypy/src/phypy_string.cpp',
                    'phypy/src/basic_lot.cpp',
                    'phypy/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phypy/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phypy/src/probability_distribution.cpp',
                    'phypy/src/slice_sampler.cpp',
                    ],**argDict),
              Extension('ReadNexus._ReadNexus', 
                    [
                    'phypy/src/read_nexus_pymod.cpp',
                    'phypy/src/boost_assertion_failed.cpp',
                    'phypy/src/cipres/CipresDataMatrixHelper.cpp',
                    'phypy/src/cipres/cipres_nexus_reader.cpp',
                    'phypy/src/ncl/nxs_basic_manager.cpp',
                    'phypy/src/ncl/nxs_block.cpp',
                    'phypy/src/ncl/nxs_exception.cpp',
                    'phypy/src/ncl/nxs_reader.cpp',
                    'phypy/src/ncl/nxs_token.cpp',
                    'phypy/src/ncl/characters/nxs_char_listener.cpp',
                    'phypy/src/ncl/characters/nxs_characters_block.cpp',
                    'phypy/src/ncl/characters/nxs_characters_block_commands.cpp',
                    'phypy/src/ncl/characters/nxs_characters_manager.cpp',
                    'phypy/src/ncl/command/nxs_auto_command.cpp',
                    'phypy/src/ncl/command/nxs_choice_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_command.cpp',
                    'phypy/src/ncl/command/nxs_command_manager.cpp',
                    'phypy/src/ncl/command/nxs_command_output.cpp',
                    'phypy/src/ncl/command/nxs_file_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_mixed_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_primitive_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_restricted_string_cmd_param.cpp',
                    'phypy/src/ncl/command/nxs_set_cmd_param.cpp',
                    'phypy/src/ncl/misc/arg_stream.cpp',
                    'phypy/src/ncl/misc/nxs_data_type.cpp',
                    'phypy/src/ncl/misc/eliminated_index_slider.cpp',
                    'phypy/src/ncl/misc/nxs_discrete_matrix.cpp',
                    'phypy/src/ncl/misc/nxs_file_path.cpp',
                    'phypy/src/ncl/misc/nxs_index_set.cpp',
                    'phypy/src/ncl/misc/nxs_test.cpp',
                    'phypy/src/ncl/misc/string_extensions.cpp',
                    'phypy/src/ncl/output/nxs_output_stream_wrapper.cpp',
                    'phypy/src/ncl/output/generic_output_classes.cpp',
                    'phypy/src/ncl/output/nxs_console_out_stream_impl.cpp',
                    'phypy/src/ncl/output/nxs_input.cpp',
                    'phypy/src/ncl/output/nxs_table.cpp',
                    'phypy/src/ncl/output/nxs_table_cell.cpp',
                    'phypy/src/ncl/output/nxs_typist.cpp',
                    'phypy/src/ncl/output/nxs_user_query.cpp',
                    'phypy/src/ncl/output/nxs_xml_socket_output_stream.cpp',
                    'phypy/src/ncl/taxa/base_taxa_manager.cpp',
                    'phypy/src/ncl/taxa/nxs_alternative_taxa_block.cpp',
                    'phypy/src/ncl/taxa/nxs_taxa_block.cpp',
                    'phypy/src/ncl/taxa/nxs_taxa_listener.cpp',
                    'phypy/src/ncl/taxa/nxs_taxa_manager.cpp',
                    'phypy/src/ncl/trees/newick_verifier.cpp',
                    'phypy/src/ncl/trees/nxs_tree_listener.cpp',
                    'phypy/src/ncl/trees/nxs_trees_block.cpp',
                    'phypy/src/ncl/trees/nxs_trees_manager.cpp',
                    'phypy/src/oldphycas/characters_manager.cpp',
                    'phypy/src/oldphycas/const_site_info.cpp',
                    'phypy/src/oldphycas/multiline_bound.cpp',
                    'phypy/src/probability_distribution.cpp',
                    'phypy/src/oldphycas/distribution_command_param.cpp',
                    'phypy/src/oldphycas/distribution_description.cpp',
                    'phypy/src/basic_lot.cpp',
                    'phypy/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phypy/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phypy/src/oldphycas/taxa_manager.cpp',
                    'phypy/src/oldphycas/draw_context.cpp',
                    'phypy/src/oldphycas/draw_tree.cpp',
                    'phypy/src/oldphycas/tree_manip.cpp',
                    'phypy/src/oldphycas/tree_node.cpp',
                    'phypy/src/oldphycas/tree.cpp',
                    'phypy/src/oldphycas/trees_manager.cpp',
                    'phypy/src/oldphycas/split.cpp',
                    ],**argDict),
              ]
    setupArgs.update({
        'ext_package':'phypy',
        'ext_modules':em,
        })

setup(**setupArgs)


# What follows is some args to Extension for the Conversions extension from an
#   aborted attempt to get Extensions built using g++ directly (instead of bjam driving
#   g++)
#      ['phypy/src/conversions_pymod.cpp',
#                     'phypy/src/boost_assertion_failed.cpp',
#                     'phypy/src/basic_tree.cpp',
#                     'phypy/src/basic_tree_node.cpp',
#                     'phypy/src/phypy_string.cpp'], 
#                     include_dirs= [
#                         os.environ.get('BOOST_ROOT'), 
#                         os.environ.get('PHYCAS_ROOT'),
#                         os.path.expandvars('$PHYCAS_ROOT/phypy/src')
#                         ],
#                     define_macros=[
#                         ('BOOST_PYTHON_DYNAMIC_LIB', '1')
#                         ],
                   
