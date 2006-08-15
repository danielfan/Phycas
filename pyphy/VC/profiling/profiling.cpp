#include <algorithm>
#include <iostream>
#include "pyphy/VC/test/test.hpp"
#include "pyphy/src/cipres/CipresDataMatrixHelper.h"
#include "pyphy/src/cipres/cipres_nexus_reader.hpp"
#include "pyphy/src/tree_likelihood.hpp"
#include "pyphy/src/likelihood_models.hpp"
#include "pyphy/src/larget_simon_move.hpp"
#include "pyphy/src/mcmc_chain_manager.hpp"
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/oldphycas/characters_manager.hpp"

using namespace phycas;

int	CipresNexusReader::GetNChar()
	{
	return phoCharactersMgr->GetNumChars();
	}

int main()
	{
	// Read the data matrix
	CipresNexusReader reader = CipresNexusReader(-1);
	//std::string file_path = "../../pyphy/green.nex";
	std::string file_path = "../../pyphy/nyldna4.nex";
	reader.ReadFilePath(file_path);
	std::vector<std::string> taxon_names = reader.GetTaxLabels();
	unsigned ntax = (unsigned)taxon_names.size();
	unsigned nchar = reader.GetNChar();
	std::cout << "Data matrix " << file_path;
	std::cout << "\n  ntax  = " << ntax;
	std::cout << "\n  nchar = " << nchar;
	std::cout << "\n  taxon labels:\n\t";
	std::copy(taxon_names.begin(), taxon_names.end(), std::ostream_iterator<std::string>(std::cout, "\n\t"));
	std::cout << std::endl;

	// Get access to the embedded data
	CipresNative::DiscreteMatrix * data_matrix = createNativeDiscreteMatrix(reader, 0);
	assert(data_matrix != NULL);
	assert(ntax == data_matrix->getNTax());
	assert(nchar == data_matrix->getNChar());

	// Create a random number generator
	phycas::LotShPtr lot(new phycas::Lot());
	lot->SetSeed(13579);

	// Create priors
	ProbDistShPtr kappa_prior(new ExponentialDistribution(1.0));
	ProbDistShPtr basefreq_param_prior(new ExponentialDistribution(1.0));
	ProbDistShPtr shape_prior(new ExponentialDistribution(1.0));
	ProbDistShPtr edgelen_prior(new ExponentialDistribution(2.0));
	ProbDistShPtr edgelen_hyperprior(new InverseGammaDistribution(2.1, 0.909));

	// Replace random number generator in each prior with lot
	kappa_prior->SetLot(lot.get());
	basefreq_param_prior->SetLot(lot.get());
	shape_prior->SetLot(lot.get());
	edgelen_prior->SetLot(lot.get());
	edgelen_hyperprior->SetLot(lot.get());

	// Create a model
	unsigned ncat = 3;
	phycas::HKYShPtr hky(new phycas::HKY());
	hky->setKappa(4.0);
	hky->setKappaPrior(kappa_prior);
	hky->setNucleotideFreqs(0.1, 0.2, 0.3, 0.4);
	hky->setBaseFreqParamPrior(basefreq_param_prior);
	hky->setNGammaRates(ncat);
	hky->setShape(0.5);
	hky->setDiscreteGammaShapePrior(shape_prior);
	hky->setPriorOnShapeInverse(false);
	hky->setNotPinvarModel();
	hky->setEdgeLenPrior(edgelen_prior);
	hky->setEdgeLenHyperPrior(edgelen_hyperprior);
	phycas::ModelShPtr model = hky;
	unsigned nrates = ncat;
	unsigned nstates = 4;

	// Create a TreeLikelihood object for the model and copy data from data_matrix
	typedef boost::shared_ptr<phycas::TreeLikelihood> TreeLikelihoodShPtr; 
	TreeLikelihoodShPtr likelihood(new phycas::TreeLikelihood(model));
	likelihood->copyDataFromDiscreteMatrix(*data_matrix);
	unsigned npatterns = likelihood->getNPatterns();
	std::cout << "Number of data patterns = " << npatterns << std::endl;
	std::cout << "Patterns found:" << std::endl;
	std::cout << likelihood->listPatterns(false) << std::endl;

	// Build a random tree
	phycas::TreeShPtr t(new phycas::Tree());
	phycas::TreeManip tree_manip(t);
	tree_manip.randomTree(ntax, lot, edgelen_prior, false);
	std::cout << "Starting tree = " << t->MakeNewick() << std::endl;

	// Prepare the tree
	likelihood->prepareForLikelihood(t);

	// Create chain, add parameters and moves
	typedef boost::shared_ptr<phycas::MCMCChainManager> MCMCChainManagerShPtr; 
	MCMCChainManagerShPtr mcmc(new phycas::MCMCChainManager());
	mcmc->addMCMCUpdaters(model, t, likelihood, lot, false, 0, 1);
	typedef boost::shared_ptr<phycas::LargetSimonMove> LargetSimonMoveShPtr; 
	LargetSimonMoveShPtr lsmove(new phycas::LargetSimonMove());
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
		const phycas::MCMCUpdaterVect & all_updaters = mcmc->getAllUpdaters();
		for (phycas::MCMCUpdaterVect::const_iterator p = all_updaters.begin(); p != all_updaters.end(); ++p)
			{
			unsigned w = (*p)->getWeight();
			for (unsigned rep = 0; rep < w; ++rep)
				{
				(*p)->update();
				}
			}
		}
	std::cout << "Log-likelihood after 1 cycle = " << mcmc->getLastLnLike() << std::endl;

#if 0
	char answer = AskUser("\nPress y to quit...");
	if (answer == '\0')
		std::cerr << "Oops, something bad happened in the AskUser function." << std::endl;
	else if (answer == 'y')
		std::cerr << "\nYou pressed y so I'm quitting" << std::endl;
	else
		std::cerr << "\nYou pressed " << answer << " but I'm quitting anyway!" << std::endl;
#endif
 
	return 0;
	}
