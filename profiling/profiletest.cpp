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
		model = hky;
#endif
	}
	unsigned nrates = ncat;
	unsigned nstates = 4;
	
	// Create a TreeLikelihood object for the model and copy data from data_matrix
	TreeLikelihoodShPtr likelihood(new TreeLikelihood(model));
	likelihood->copyDataFromDiscreteMatrix(*data_matrix);
	unsigned npatterns = likelihood->getNPatterns();
	std::cout << "Number of data patterns = " << npatterns << std::endl;
	std::cout << "Patterns found:" << std::endl;
	std::cout << likelihood->listPatterns(false) << std::endl;

	// Build a random tree
	t.reset(new Tree());
	TreeManip tree_manip(t);
	tree_manip.equiprobTree(ntax, lot, edgelen_prior, edgelen_prior);
	std::cout << "Starting tree = " << t->MakeNewick() << std::endl;

	// Prepare the tree
	likelihood->prepareForLikelihood(t);
	
	double lnL = likelihood->calcLnL(t);
	std::cout << "lnL of starting tree = " << lnL << std::endl;

	// Create chain, add parameters and moves
	mcmc.reset(new MCMCChainManager());
	mcmc->addMCMCUpdaters(model, t, likelihood, lot, 0, 1);	

	lsmove.reset(new LargetSimonMove());
	lsmove->setName("Larget-Simon move");
    lsmove->setWeight(100);
    lsmove->setTree(t);
    lsmove->setModel(model);
    lsmove->setTreeLikelihood(likelihood);
    lsmove->setLot(lot);
    lsmove->setLambda(0.2);
	mcmc->addMove(lsmove);
	mcmc->finalize();
	
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
	// 	for (unsigned k = 0; k < 10; ++k)
	// 		{
	// 		tree_manip.equiprobTree(ntax, lot, edgelen_prior, edgelen_prior);
	// 		likelihood->prepareForLikelihood(t);
	// 		double lnL = likelihood->calcLnL(t);
	// 		std::cout << "lnL of starting tree = " << lnL << std::endl;
	// 		}
		
	model->releaseUpdaters();
	}
	
int main()
	{
	CREATE_MEMCHK
	run();
 	MEMCHK_REPORT(std::cerr)
	return 0;
	}
	
