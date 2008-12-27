# Assumes these environmental variables have been defined
# PHYCAS_ROOT is where the phycas/src folder (containing the C++ sources for building phycas) resides (e.g. $(HOME)/pdev/trunk)
# NCL_ROOT is where the ncl folder (containing the source files for building NCL) resides (e.g. $(HOME)/ndev/branches/v2.1)
# BOOST_ROOT is where the boost folder resides (e.g. boost_1_37_0)

TARGETS = profiletest treetest internaldatatest
all: $(TARGETS)

clean: 
	rm -f *.o
	rm -f  $(TARGETS)

INCLUDE_DIRS =                   \
   -I$(PHYCAS_ROOT)              \
   -I$(PHYCAS_ROOT)/phycas/src   \
   -I$(NCL_ROOT)                 \
   -I$(BOOST_ROOT)
FORCE_INCLUDE = -include $(PHYCAS_ROOT)/profiling/test_force_incl.hpp
PHYCAS_THIRDPARTY_SRC = $(PHYCAS_ROOT)//phycas/src/thirdparty/dcdflib/src/
PHYCAS_SRC = $(PHYCAS_ROOT)/phycas/src
NCL_SRC = $(NCL_ROOT)/ncl

# Standard Make variables
VPATH = $(PHYCAS_SRC):$(PHYCAS_THIRDPARTY_SRC):$(NCL_SRC)
CXXFLAGS = $(INCLUDE_DIRS) $(FORCE_INCLUDE)

# Rules for compiling the profile target
PROFILETEST_OBJS = profiletest.o nxs_file_path.o phycas_nexus_reader.o nxsreader.o nxstoken.o nxsblock.o \
					nxspublicblocks.o nxsassumptionsblock.o nxscharactersblock.o nxstaxablock.o nxsexception.o \
					nxssetreader.o nxsstring.o nxstreesblock.o nxsunalignedblock.o nxscxxdiscretematrix.o nxsdatablock.o \
					nxsdistancesblock.o mcmc_updater.o slice_sampler.o tree_likelihood.o likelihood_models.o larget_simon_move.o \
					mcmc_chain_manager.o probability_distribution.o univents.o internal_data.o tip_data.o univent_prob_mgr.o \
					basic_tree.o basic_tree_node.o edgelen_master_param.o likelihood_loops.o split.o square_matrix.o \
					basic_lot.o basic_cdf.o dcdflib.o ipmpar.o underflow_policy.o flex_rate_param.o flex_prob_param.o \
					pinvar_param.o mapping_move.o tree_manip.o hyperprior_param.o mcmc_param.o state_freq_param.o kappa_param.o \
					jc_model.o hky_model.o gtr_model.o codon_model.o q_matrix.o omega_param.o sim_data.o gtr_rate_param.o \
					discrete_gamma_shape_param.o linalg.o cond_likelihood_storage.o mcmc_flexcat_param.o
profiletest: test_force_incl.hpp $(PROFILETEST_OBJS)
	$(CXX) $(CXXFLAGS) -o profiletest $(PROFILETEST_OBJS)

# Rules for compiling the internaldatatest target
INTERNALDATATEST_OBJS = internaldatatest.o  internal_data.o univents.o
internaldatatest: test_force_incl.hpp $(INTERNALDATATEST_OBJS)
	$(CXX) $(CXXFLAGS) -o internaldatatest $(INTERNALDATATEST_OBJS)

# Rules for compiling the treetest target
TREETEST_OBJS = treetest.o  basic_cdf.o basic_lot.o basic_tree.o basic_tree_node.o tree_manip.o \
				split.o probability_distribution.o dcdflib.o ipmpar.o
treetest: test_force_incl.hpp $(TREETEST_OBJS)
	$(CXX) $(CXXFLAGS) -o treetest $(TREETEST_OBJS)

