/*
 *  beaglelib.h
 *  phycas
 *
 *  Created by Daniel on 2/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#if ! defined(BEAGLELIB_HPP)
#define BEAGLELIB_HPP

#include <vector>
#include "boost/shared_ptr.hpp"

namespace phycas {

class BeagleLib
{
	public:
	BeagleLib();
	~BeagleLib();
	
	void Init(unsigned nTaxa, unsigned nCat, unsigned nStates, unsigned nPatterns);
	void ListResources();
	void SetStateFrequencies(const std::vector<double> &freqs);
	void SetTipStates();
	void SetCategoryRatesAndWeights(const std::vector<double> &rates, const std::vector<double> &weights);
	void SetPatternWeights(const std::vector<double> &patternWeights);
	void SetEigenDecomposition(const std::vector<double> &eigenValues, const std::vector<double> &eigenVectors, const std::vector<double> &inverseEigenVectors);
	
	int	_instance;
	unsigned _nTaxa;
	unsigned _nCat;
	unsigned _nStates;
	unsigned _nPatterns;

	std::vector<int> _operations;
	std::vector<int> _pMatrixIndex;
	std::vector<double> _brLens;
};

typedef boost::shared_ptr<BeagleLib> BeagleLibShPtr;

}

#endif 