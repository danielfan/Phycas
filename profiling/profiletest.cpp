#define SIMPLEST_CASE

#include <algorithm>
#include <iostream>
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/larget_simon_move.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/phycas_nexus_reader.hpp"
#include "memchk.hpp"

using namespace phycas;

typedef boost::shared_ptr<LargetSimonMove> LargetSimonMoveShPtr; 
typedef boost::shared_ptr<MCMCChainManager> MCMCChainManagerShPtr; 
typedef boost::shared_ptr<TreeLikelihood> TreeLikelihoodShPtr; 
typedef boost::shared_ptr<NxsCXXDiscreteMatrix> NxsCXXDiscreteMatrixShPtr;
	
void run()
	{
	LotShPtr lot;
	ProbDistShPtr edgelen_prior;

#if !defined(SIMPLEST_CASE)
	ProbDistShPtr kappa_prior;
	ProbDistShPtr state_freq_param_prior;
	ProbDistShPtr shape_prior;
	ProbDistShPtr edgelen_hyperprior;
#endif

	ModelShPtr model;
	//HKYShPtr hky;
	//TreeLikelihoodShPtr likelihood;
	TreeShPtr t;
	MCMCChainManagerShPtr mcmc;
	LargetSimonMoveShPtr lsmove;
	
	// Read the data matrix
	PhycasNexusReaderShPtr reader(new PhycasNexusReader(-1));
	std::string file_path = "green.nex";
	reader->ReadFilePath(file_path);
	
	std::vector<std::string> taxon_names = reader->GetTaxLabels();
	unsigned ntax = (unsigned)taxon_names.size();
	unsigned nchar = reader->GetNChar();
	std::cout << "Data matrix " << file_path;
	std::cout << "\n  ntax  = " << ntax;
	std::cout << "\n  nchar = " << nchar;
	std::cout << "\n  taxon labels:\n\t";
	std::copy(taxon_names.begin(), taxon_names.end(), std::ostream_iterator<std::string>(std::cout, "\n\t"));
	std::cout << std::endl;

	// Get access to the embedded data
	bool convert_gaps_to_missing = true;
	NxsCXXDiscreteMatrixShPtr data_matrix = NxsCXXDiscreteMatrixShPtr(GetLastDiscreteMatrix(*reader, convert_gaps_to_missing));
	PHYCAS_ASSERT(data_matrix != NULL);
	PHYCAS_ASSERT(ntax == data_matrix->getNTax());
	PHYCAS_ASSERT(nchar == data_matrix->getNChar());

	// Create a random number generator
	lot.reset(new Lot());
	lot->SetSeed(13579);
	
	// Create priors
#if !defined(SIMPLEST_CASE)
	kappa_prior.reset(new ExponentialDistribution(1.0));
	state_freq_param_prior.reset(new ExponentialDistribution(1.0));
	shape_prior.reset(new ExponentialDistribution(1.0));
	edgelen_hyperprior.reset(new InverseGammaDistribution(2.1, 0.909));
#endif
	edgelen_prior.reset(new ExponentialDistribution(2.0));

	// Replace random number generator in each prior with lot
#if !defined(SIMPLEST_CASE)
	kappa_prior->SetLot(lot.get());
	state_freq_param_prior->SetLot(lot.get());
	shape_prior->SetLot(lot.get());
	edgelen_hyperprior->SetLot(lot.get());
#endif
	edgelen_prior->SetLot(lot.get());
	
	// Create a model
	unsigned ncat = 3;
	{
#if defined(SIMPLEST_CASE)
		JCShPtr jc(new JC());
		jc->setInternalEdgeLenPrior(edgelen_prior);
		jc->setExternalEdgeLenPrior(edgelen_prior);
		model = jc;
#else
		HKYShPtr hky(new HKY());
		hky->setKappa(4.0);
		hky->setKappaPrior(kappa_prior);
		hky->setNucleotideFreqs(0.1, 0.2, 0.3, 0.4);
		hky->setStateFreqParamPrior(state_freq_param_prior);
		hky->setNGammaRates(ncat);
		hky->setShape(0.5);
		hky->setDiscreteGammaShapePrior(shape_prior);
		hky->setPriorOnShapeInverse(false);
		hky->setNotPinvarModel();
		hky->setInternalEdgeLenPrior(edgelen_prior);
		hky->setExternalEdgeLenPrior(edgelen_prior);
		hky->setEdgeLenHyperPrior(edgelen_hyperprior);
		// model has use count 0 and hky has use count 1 here
		model = hky;
		// both model and hky have use count 2 here
		//std::cerr << "hky use count = " << hky.use_count() << std::endl;
		//std::cerr << "model use count = " << model.use_count() << std::endl;
#endif
	}
	// model has use count 1 here
	unsigned nrates = ncat;
	unsigned nstates = 4;
	
	// Create a TreeLikelihood object for the model and copy data from data_matrix
	TreeLikelihoodShPtr likelihood(new TreeLikelihood(model));
	//likelihood.reset(new TreeLikelihood(model));
	// likelihood has use count 1 here
	// model has use count 3 here
	//	1 from original new
	//	1 from copy to model data member in TreeLikelihood constructor
	//  1 from copy to univentProbMgr data member in TreeLikelihood constructor
	likelihood->copyDataFromDiscreteMatrix(*data_matrix);
	unsigned npatterns = likelihood->getNPatterns();
	std::cout << "Number of data patterns = " << npatterns << std::endl;
	std::cout << "Patterns found:" << std::endl;
	std::cout << likelihood->listPatterns(false) << std::endl;

	// Build a random tree
	t.reset(new Tree());
	TreeManip tree_manip(t);
	std::cerr << "\n***** edgelen_prior use count before equiprobTree = " << edgelen_prior.use_count() << " (expecting 3)" << std::endl;
	tree_manip.equiprobTree(ntax, lot, edgelen_prior, edgelen_prior);
	std::cerr << "\n***** edgelen_prior use count after equiprobTree = " << edgelen_prior.use_count() << " (expecting 3)" << std::endl;
	std::cout << "Starting tree = " << t->MakeNewick() << std::endl;

	// Prepare the tree
	likelihood->prepareForLikelihood(t);
	
	double lnL = likelihood->calcLnL(t);
	std::cout << "lnL of starting tree = " << lnL << std::endl;

	std::cerr << "\n***** edgelen_prior use count before creating mcmc = " << edgelen_prior.use_count() << " (expecting 3)" << std::endl;
	// model has use count 3 here
	
	// Create chain, add parameters and moves
	mcmc.reset(new MCMCChainManager());
	std::cerr << "\n***** lot use count before addMCMCUpdaters = " << lot.use_count() << " (expecting 1)" << std::endl;
	std::cerr << "\n***** edgelen_prior use count after creating mcmc = " << edgelen_prior.use_count() << " (expecting 3)" << std::endl;
	mcmc->addMCMCUpdaters(model, t, likelihood, lot, 0, 1);	
	mcmc->debugUpdaterReport("after addMCMCUpdaters");
	
	std::cerr << "\n***** edgelen_prior use count after addMCMCUpdaters = " << edgelen_prior.use_count() << " (expecting 4)" << std::endl;
	std::cerr << "\n***** lot use count after addMCMCUpdaters = " << lot.use_count() << " (expecting 2)" << std::endl;

	// likelihood has use count 9 here (1 from before plus 8 updaters (1 kappa + 4 base freqs +  1 shape + 1 edgelen master + 1 hyper) held by mcmc)
	// model has use count 11 here (3 from before plus 8 updaters held by mcmc)
	
	//   9 undeleted memory elements introduced here
	// -----------------------------------------------------------
	//  12 bytes remaining from allocation at profiletest.cpp (40) lot
	//  48 bytes remaining from allocation at profiletest.cpp (44) kappa_prior
	//  48 bytes remaining from allocation at profiletest.cpp (45) state_freq_param_prior
	//  48 bytes remaining from allocation at profiletest.cpp (46) shape_prior
	//  48 bytes remaining from allocation at profiletest.cpp (47) edgelen_prior
	//  48 bytes remaining from allocation at profiletest.cpp (48) edgelen_hyperprior
	// 308 bytes remaining from allocation at profiletest.cpp (59) hky
	// 352 bytes remaining from allocation at profiletest.cpp (78) likelihood
	// 116 bytes remaining from allocation at profiletest.cpp (86) t

	lsmove.reset(new LargetSimonMove());
	// lsmove now has use count = 1 
	lsmove->setName("Larget-Simon move");
    lsmove->setWeight(100);
    lsmove->setTree(t);
    lsmove->setModel(model);
    lsmove->setTreeLikelihood(likelihood);
    lsmove->setLot(lot);
    lsmove->setLambda(0.2);
	mcmc->addMove(lsmove);
	// lsmove now has use count = 2 because a second shared pointer has been added to MCMCChainManager's moves vector
	mcmc->finalize();
	mcmc->debugUpdaterReport("after finalize");	
	// lsmove now has use count = 3 because a third shared pointer has been added to MCMCChainManager's all_updaters vector
	// likelihood now has use count = 10: 1 (new) + 8 held by updaters and 1 held by lsmove

	std::cerr << "\n***** lot use count after finalize = " << lot.use_count() << " (expecting 3)" << std::endl;
	
	// Compute starting log-likelihood
	mcmc->refreshLastLnLike();
	std::cout << "Starting log-likelihood = " << mcmc->getLastLnLike() << std::endl;

	// Run chain for 1 cycle
	for (unsigned cycle = 0; cycle < 1; ++cycle)
		{
		const MCMCUpdaterVect & all_updaters = mcmc->getAllUpdaters();
		for (MCMCUpdaterVect::const_iterator p = all_updaters.begin(); p != all_updaters.end(); ++p)
			{
			std::cout << "Updating " << (*p)->getName() << "..." << std::endl;
			unsigned w = (*p)->getWeight();
			for (unsigned rep = 0; rep < w; ++rep)
				{
				(*p)->update();
				}
			}
		}
	std::cout << "Log-likelihood after 1 cycle = " << mcmc->getLastLnLike() << std::endl;
	
	// Build 10 equiprobable trees and compute the log-likelihood for each
	for (unsigned k = 0; k < 10; ++k)
		{
		tree_manip.equiprobTree(ntax, lot, edgelen_prior, edgelen_prior);
		likelihood->prepareForLikelihood(t);
		double lnL = likelihood->calcLnL(t);
		std::cout << "lnL of starting tree = " << lnL << std::endl;
		}
		
	//mcmc->releaseUpdaters();
	//likelihood->releaseModel();
	model->releaseUpdaters();
	
	// 	************ before exiting run function ***********
	// 	lot                    use count = 3
	// 	edgelen_prior          use count = 2 vs 4
	// 	likelihood             use count = 3
	// 	t                      use count = 6
	// 	model                  use count = 4 vs 5
	// 	mcmc                   use count = 1
	// 	lsmove                 use count = 3 vs. 1

	std::cerr << "\n\n************ before exiting run function ***********" << std::endl;
	std::cerr << "lot                    use count = " << lot.use_count() << std::endl;
	std::cerr << "edgelen_prior          use count = " << edgelen_prior.use_count() << std::endl;
#if !defined(SIMPLEST_CASE)
	std::cerr << "kappa_prior            use count = " << kappa_prior.use_count() << std::endl;
	std::cerr << "state_freq_param_prior use count = " << state_freq_param_prior.use_count() << std::endl;
	std::cerr << "shape_prior            use count = " << shape_prior.use_count() << std::endl;
	std::cerr << "edgelen_hyperprior     use count = " << edgelen_hyperprior.use_count() << std::endl;
#endif
	//std::cerr << "hky                    use count = " << hky.use_count() << std::endl;
	std::cerr << "likelihood             use count = " << likelihood.use_count() << std::endl;
	std::cerr << "t                      use count = " << t.use_count() << std::endl;
	std::cerr << "model                  use count = " << model.use_count() << std::endl;
	std::cerr << "mcmc                   use count = " << mcmc.use_count() << std::endl;
	std::cerr << "lsmove                 use count = " << lsmove.use_count() << std::endl;
	mcmc->debugUpdaterReport("before exiting run");
	
#if 0
	char answer = AskUser("\nPress y to quit...");
	if (answer == '\0')
		std::cerr << "Oops, something bad happened in the AskUser function." << std::endl;
	else if (answer == 'y')
		std::cerr << "\nYou pressed y so I'm quitting" << std::endl;
	else
		std::cerr << "\nYou pressed " << answer << " but I'm quitting anyway!" << std::endl;
#endif
	}
	
typedef boost::shared_ptr<std::string> StrShPtr;
	
class Test
	{
	public:
		Test(StrShPtr p) : myp(p)
			{
			}
		~Test()
			{
			//std::cerr << "Test dying...myp use count = " << myp.use_count() << std::endl;
			}
	private:
		StrShPtr myp;
	};
	
typedef boost::shared_ptr<Test> TestShPtr;
	
void test()
	{
	StrShPtr p(new std::string("hello"));
	std::cerr << "p use count = " << p.use_count() << std::endl;
	{
		TestShPtr t(new Test(p));
		std::cerr << "p use count after Test object created = " << p.use_count() << std::endl;
	}
	std::cerr << "p use count after Test object deleted = " << p.use_count() << std::endl;
	}

int main()
	{
	CREATE_MEMCHK
	//test();
	run();
 	MEMCHK_REPORT(std::cerr)
	return 0;
	}
	
//   12 bytes remaining from allocation at profiletest.cpp (57)		17	lot
//   48 bytes remaining from allocation at profiletest.cpp (61)		3	kappa_prior
//		1 is original new
//		1 is from transfer to hky model
//		1 is from transfer from model to updater
//   48 bytes remaining from allocation at profiletest.cpp (62)		6	state_freq_param_prior
//   48 bytes remaining from allocation at profiletest.cpp (63)		3	shape_prior
//   48 bytes remaining from allocation at profiletest.cpp (64)		4	edgelen_prior
//   48 bytes remaining from allocation at profiletest.cpp (65)		3	edgelen_hyperprior
//   308 bytes remaining from allocation at profiletest.cpp (76)	13	hky
//		1 is original new
//		1 is from copy to model
//		1 is from 
//   352 bytes remaining from allocation at profiletest.cpp (94)	10	likelihood
//   116 bytes remaining from allocation at profiletest.cpp (102)	20	t
