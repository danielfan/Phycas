import sys

# If build_number_from_svn_info = True, gleans svn revision from subprocess call to "svn info -r HEAD"
# if False, need to set svn revision manually in svn_revision
build_number_from_svn_info = False

# the following setting is only used if build_number_from_svn_info is False, or
# regular expression search of svn output fails to find pattern 'Revision: (\d+)'
svn_revision = 760

phycas_major = '0'
phycas_minor = '11'

boost_target_dir = 'boost_1_34_0'

fakeCompilingExtensions = True # deprecated
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
          'PDFGen/AFM/*.afm',
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
          'Tests/FixedParams/reference_output/fixed.p',
          'Tests/FixedParams/reference_output/fixed.t',
          'Tests/FixedTopology/__init__.py',
          'Tests/FixedTopology/FixedTopology.py',
          'Tests/FixedTopology/reference_output/fixdtree.p',
          'Tests/FixedTopology/reference_output/fixdtree.t',
          'Tests/FixedTopology/reference_output/simulated.nex',
          'Tests/GTRTest/__init__.py',
          'Tests/GTRTest/GTRTest.py',
          'Tests/GTRTest/reference_output/gtr_test.p',
          'Tests/GTRTest/reference_output/gtr_test.t',
          #'Tests/GelfandGhosh/__init__.py',
          #'Tests/GelfandGhosh/GelfandGhosh.py',
          #'Tests/GelfandGhosh/reference_output/analHKY.nex.p',
          #'Tests/GelfandGhosh/reference_output/analHKY.nex.t',
          #'Tests/GelfandGhosh/reference_output/analHKYflex.nex.p',
          #'Tests/GelfandGhosh/reference_output/analHKYflex.nex.t',
          #'Tests/GelfandGhosh/reference_output/analHKYg.nex.p',
          #'Tests/GelfandGhosh/reference_output/analHKYg.nex.t',
          #'Tests/GelfandGhosh/reference_output/ggout.txt',
          #'Tests/LikelihoodTest/__init__.py',
          #'Tests/LikelihoodTest/LikelihoodTest.py',
          #'Tests/LikelihoodTest/reference_output/check.nex',
          #'Tests/LikelihoodTest/reference_output/simulated.nex',
          #'Tests/Polytomies/__init__.py',
          #'Tests/Polytomies/Polytomies.py',
          #'Tests/Polytomies/reference_output/HKYpolytomy.p',
          #'Tests/Polytomies/reference_output/HKYpolytomy.t',
          #'Tests/Polytomies/reference_output/simHKY.nex',
          'Tests/Simulator/__init__.py',
          'Tests/Simulator/Simulator.py',
          'Tests/Simulator/reference_output/simulated.nex',
          'Tests/SplitTest/__init__.py',
          'Tests/SplitTest/SplitTest.py',
          'Tests/SplitTest/reference_output/out.txt',
          'Tests/PDFTree/__init__.py',
          'Tests/PDFTree/PDFTree.py',
          'Tests/PDFTree/reference_output/test.pdf',
          'Tests/SumT/__init__.py',
          'Tests/SumT/SumT.py',
          'Tests/SumT/test.t',
          'Tests/SumT/reference_output/trees.tre',
          'Tests/SumT/reference_output/trees.pdf',
          'Tests/SumT/reference_output/splits.pdf',
          'Tests/SumT/reference_output/logfile.txt',
          'Tests/*.py',
          'Examples/Paradox/Paradox.py',
          'Examples/Paradox/ShoupLewis.nex',
          'Examples/Simplest/Simplest.py',
          'Examples/Simplest/green.nex'
          ]

#data_windows_only = [
#                  'Tests/cleanall.bat',
#                  'Tests/doctestall.bat',
#                  'Tests/runall.bat'
#                  ]
#data_windows_only.extend(data_all)

windows_package_data = {
                    'phycas.Conversions': ['*.pyd','*.dll'],
                    'phycas.DataMatrix': ['*.pyd'],
                    'phycas.Likelihood': ['*.pyd'],
                    'phycas.ProbDist': ['*.pyd'],
                    'phycas.Phylogeny': ['*.pyd'],
                    'phycas.ReadNexus': ['*.pyd'],
                    'phycas.PDFGen': ['AFM/*.afm']
                    }
#windows_package_data.update({'phycas': data_windows_only})
windows_package_data.update({'phycas': data_all})

isWin = sys.platform == 'win32'

sharedObjSuffix = [distutils.sysconfig.get_config_var("SO")]
dynamicLibSuffix = isWin and ['.dll'] or [".dylib"]
conversionsDataFiles = sharedObjSuffix + dynamicLibSuffix

test_data = [
          'Tests/cleanall.py',
          'Tests/doctestall.py',
          'Tests/runall.py',
          ]
_pack_data = {
        'phycas': data_all,
        }
phycas_description = """\
Phycas:
 	
Phycas is a Python application for carrying out phylogenetic analyses.
It is also a C++ and Python library that can be used to create new
applications or to extend the current functionality.
"""

if build_number_from_svn_info:
    import subprocess, re
    svnheadinfo = subprocess.Popen('svn info -r HEAD', shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
    re_match = re.search('Revision: (\d+)', svnheadinfo)
    if re_match:
        svninfo = subprocess.Popen('svn status --non-interactive', shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
        if len(svninfo) == 0:
            # we're up to date
            svn_revision = re_match.group(1)
        else:
            # we're not up to date, once we commit (which hopefully will happen before
            # any other developer commits!), revision number will increase by 1
            next_revision_number = 1 + int(re_match.group(1))
            svn_revision = str(next_revision_number)
            
phycas_full_version = phycas_major+'.'+phycas_minor+'.'+str(svn_revision)
        
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
     'phycas.wxPhycas',
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
        "/usr/local/include/"+boost_target_dir, 
        os.path.abspath(parent_dir),
        os.path.abspath(os.path.join(parent_dir, "phycas", "src")),
        ]
    libraries=["boost_python"]
    library_dirs=['/usr/local/lib']

if compilingExt:
    if fakeCompilingExtensions:
    	phycasRoot = os.path.dirname(sys.argv[0])
    	builtUnder = "--built-under=" + phycasRoot
        em = [Extension('Conversions._Conversions', ['phycas/src/conversions_pymod.cpp'], 
                    extra_link_args = [
                        #'--colocate-lib=libboost_python.dylib',
                        builtUnder,
                        '--path-from-package=phycas/Conversions', ]),
            Extension('DataMatrix._DataMatrixBase', ['phycas/src/data_matrix_pymod.cpp'], 
                    extra_link_args = [
                        builtUnder,
                        '--path-from-package=phycas/DataMatrix', ]),
            Extension('Likelihood._LikelihoodBase', ['phycas/src/likelihood_pymod.cpp'], 
                    extra_link_args = [
                        builtUnder,
                        '--path-from-package=phycas/Likelihood', ]),
            Extension('Phylogeny._Phylogeny', ['phycas/src/phylogeny_pymod.cpp'], 
                    extra_link_args = [
                        builtUnder,
                        '--path-from-package=phycas/Phylogeny', ]),
            Extension('ProbDist._ProbDist', ['phycas/src/probdist_pymod.cpp'], 
                    extra_link_args = [
                        builtUnder,
                        '--path-from-package=phycas/ProbDist', ]),
            Extension('ReadNexus._ReadNexus', ['phycas/src/read_nexus_pymod.cpp'], 
                    extra_link_args = [
                        builtUnder,
                        '--path-from-package=phycas/ReadNexus', ]),
            ]
    else:
        import subprocess
        arch = subprocess.Popen(["arch"], shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
        argDict = {
                "library_dirs": library_dirs,
                "libraries": libraries,
                "include_dirs": include_dirs,
                "depends": [],
                "extra_compile_args": [
                        "-include", 
                        "phycas/src/phycas_config.h", 
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
                    "phycas/src/conversions_pymod.cpp",
                    "phycas/src/boost_assertion_failed.cpp",
                    "phycas/src/basic_tree.cpp",
                    "phycas/src/basic_tree_node.cpp",
                    "phycas/src/phycas_string.cpp"
                    ], **argDict),
              Extension("DataMatrix._DataMatrixBase", 
                    [
                    "phycas/src/data_matrix_pymod.cpp",
                    "phycas/src/boost_assertion_failed.cpp",
                    "phycas/src/cipres/CipresDataMatrixHelper.cpp",
                    ],**argDict),
              Extension('Likelihood._LikelihoodBase', 
                    [
                    'phycas/src/basic_lot.cpp',
                    'phycas/src/basic_tree.cpp',
                    'phycas/src/basic_tree_node.cpp',
                    'phycas/src/boost_assertion_failed.cpp',
                    'phycas/src/bush_move.cpp',
                    'phycas/src/cipres/CipresDataMatrixHelper.cpp',
                    'phycas/src/cond_likelihood_storage.cpp',
                    'phycas/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phycas/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phycas/src/edge_move.cpp',
                    'phycas/src/internal_data.cpp',
                    'phycas/src/larget_simon_move.cpp',
                    'phycas/src/likelihood_loops.cpp',
                    'phycas/src/likelihood_models.cpp',
                    'phycas/src/likelihood_pymod.cpp',
                    'phycas/src/model_pymod.cpp',
                    'phycas/src/linalg.cpp',
                    'phycas/src/mcmc_chain_manager.cpp',
                    'phycas/src/mcmc_param.cpp',
                    'phycas/src/mcmc_flexcat_param.cpp',
                    'phycas/src/mcmc_updater.cpp',
                    'phycas/src/ncat_move.cpp',
                    'phycas/src/phycas_string.cpp',
                    'phycas/src/q_matrix.cpp',
                    'phycas/src/sim_data.cpp',
                    'phycas/src/slice_sampler.cpp',
                    'phycas/src/tip_data.cpp',
                    'phycas/src/topo_prior_calculator.cpp',
                    'phycas/src/tree_likelihood.cpp',
                    'phycas/src/tree_manip.cpp',
                    'phycas/src/underflow_policy.cpp',
                    'phycas/src/updater_pymod.cpp',
                    ], **argDict),
              Extension('Phylogeny._Phylogeny', 
                    [
                    'phycas/src/phylogeny_pymod.cpp',
                    'phycas/src/boost_assertion_failed.cpp',
                    'phycas/src/basic_tree.cpp',
                    'phycas/src/basic_tree_node.cpp',
                    'phycas/src/tree_manip.cpp',
                    'phycas/src/basic_lot.cpp',
                    'phycas/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phycas/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phycas/src/basic_cdf.cpp',
                    'phycas/src/phycas_string.cpp', 
                    ],**argDict),
              Extension('ProbDist._ProbDist', 
                    [
                    'phycas/src/probdist_pymod.cpp',
                    'phycas/src/boost_assertion_failed.cpp',
                    'phycas/src/phycas_string.cpp',
                    'phycas/src/basic_lot.cpp',
                    'phycas/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phycas/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phycas/src/probability_distribution.cpp',
                    'phycas/src/slice_sampler.cpp',
                    ],**argDict),
              Extension('ReadNexus._ReadNexus', 
                    [
                    'phycas/src/read_nexus_pymod.cpp',
                    'phycas/src/boost_assertion_failed.cpp',
                    'phycas/src/cipres/CipresDataMatrixHelper.cpp',
                    'phycas/src/cipres/cipres_nexus_reader.cpp',
                    'phycas/src/ncl/nxs_basic_manager.cpp',
                    'phycas/src/ncl/nxs_block.cpp',
                    'phycas/src/ncl/nxs_exception.cpp',
                    'phycas/src/ncl/nxs_reader.cpp',
                    'phycas/src/ncl/nxs_token.cpp',
                    'phycas/src/ncl/characters/nxs_char_listener.cpp',
                    'phycas/src/ncl/characters/nxs_characters_block.cpp',
                    'phycas/src/ncl/characters/nxs_characters_block_commands.cpp',
                    'phycas/src/ncl/characters/nxs_characters_manager.cpp',
                    'phycas/src/ncl/command/nxs_auto_command.cpp',
                    'phycas/src/ncl/command/nxs_choice_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_command.cpp',
                    'phycas/src/ncl/command/nxs_command_manager.cpp',
                    'phycas/src/ncl/command/nxs_command_output.cpp',
                    'phycas/src/ncl/command/nxs_file_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_mixed_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_primitive_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_restricted_string_cmd_param.cpp',
                    'phycas/src/ncl/command/nxs_set_cmd_param.cpp',
                    'phycas/src/ncl/misc/arg_stream.cpp',
                    'phycas/src/ncl/misc/nxs_data_type.cpp',
                    'phycas/src/ncl/misc/eliminated_index_slider.cpp',
                    'phycas/src/ncl/misc/nxs_discrete_matrix.cpp',
                    'phycas/src/ncl/misc/nxs_file_path.cpp',
                    'phycas/src/ncl/misc/nxs_index_set.cpp',
                    'phycas/src/ncl/misc/nxs_test.cpp',
                    'phycas/src/ncl/misc/string_extensions.cpp',
                    'phycas/src/ncl/output/nxs_output_stream_wrapper.cpp',
                    'phycas/src/ncl/output/generic_output_classes.cpp',
                    'phycas/src/ncl/output/nxs_console_out_stream_impl.cpp',
                    'phycas/src/ncl/output/nxs_input.cpp',
                    'phycas/src/ncl/output/nxs_table.cpp',
                    'phycas/src/ncl/output/nxs_table_cell.cpp',
                    'phycas/src/ncl/output/nxs_typist.cpp',
                    'phycas/src/ncl/output/nxs_user_query.cpp',
                    'phycas/src/ncl/output/nxs_xml_socket_output_stream.cpp',
                    'phycas/src/ncl/taxa/base_taxa_manager.cpp',
                    'phycas/src/ncl/taxa/nxs_alternative_taxa_block.cpp',
                    'phycas/src/ncl/taxa/nxs_taxa_block.cpp',
                    'phycas/src/ncl/taxa/nxs_taxa_listener.cpp',
                    'phycas/src/ncl/taxa/nxs_taxa_manager.cpp',
                    'phycas/src/ncl/trees/newick_verifier.cpp',
                    'phycas/src/ncl/trees/nxs_tree_listener.cpp',
                    'phycas/src/ncl/trees/nxs_trees_block.cpp',
                    'phycas/src/ncl/trees/nxs_trees_manager.cpp',
                    'phycas/src/oldphycas/characters_manager.cpp',
                    'phycas/src/oldphycas/const_site_info.cpp',
                    'phycas/src/oldphycas/multiline_bound.cpp',
                    'phycas/src/probability_distribution.cpp',
                    'phycas/src/oldphycas/distribution_command_param.cpp',
                    'phycas/src/oldphycas/distribution_description.cpp',
                    'phycas/src/basic_lot.cpp',
                    'phycas/src/thirdparty/dcdflib/src/dcdflib.cpp',
                    'phycas/src/thirdparty/dcdflib/src/ipmpar.cpp',
                    'phycas/src/oldphycas/taxa_manager.cpp',
                    'phycas/src/oldphycas/draw_context.cpp',
                    'phycas/src/oldphycas/draw_tree.cpp',
                    'phycas/src/oldphycas/tree_manip.cpp',
                    'phycas/src/oldphycas/tree_node.cpp',
                    'phycas/src/oldphycas/tree.cpp',
                    'phycas/src/oldphycas/trees_manager.cpp',
                    'phycas/src/oldphycas/split.cpp',
                    ],**argDict),
              ]
    setupArgs.update({
        'ext_package':'phycas',
        'ext_modules':em,
        })

setup(**setupArgs)


# What follows is some args to Extension for the Conversions extension from an
#   aborted attempt to get Extensions built using g++ directly (instead of bjam driving
#   g++)
#      ['phycas/src/conversions_pymod.cpp',
#                     'phycas/src/boost_assertion_failed.cpp',
#                     'phycas/src/basic_tree.cpp',
#                     'phycas/src/basic_tree_node.cpp',
#                     'phycas/src/phycas_string.cpp'], 
#                     include_dirs= [
#                         os.environ.get('BOOST_ROOT'), 
#                         os.environ.get('PHYCAS_ROOT'),
#                         os.path.expandvars('$PHYCAS_ROOT/phycas/src')
#                         ],
#                     define_macros=[
#                         ('BOOST_PYTHON_DYNAMIC_LIB', '1')
#                         ],
                   
