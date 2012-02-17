/*
 *  beaglelib.cpp
 *  phycas
 *
 *  Created by Daniel on 2/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */	
#include "beaglelib.hpp"
//BEAGLELIB
#include "libhmsbeagle/beagle.h"

using namespace phycas;

BeagleLib::BeagleLib():_instance(0), _nTaxa(0), _nCat(0), _nStates(0), _nPatterns(0) {
}

BeagleLib::~BeagleLib() {
	std::cerr << "Inside BeagleLib destructor.\n"; 
	if (_instance) {
		int code = beagleFinalizeInstance(_instance);				
		if(code != 0) {
			std::cout << "Fail to finalize instance.\n";
		}
	}
}

void BeagleLib::ListResources() {
	// list all the resources
	//
    BeagleResourceList* rsrcList;
    rsrcList = beagleGetResourceList();
	std::cout << "Available resources:\n";
    for(unsigned i = 0; i < rsrcList->length; ++i) {
		std::cout << "\tResource " << i << ":\n\t\tName : " << rsrcList->list[i].name << '\n';
		std::cout << "\t\tDesc : " << rsrcList->list[i].description << '\n';
    }    
}

void BeagleLib::Init(unsigned nTaxa, unsigned nCat, unsigned nStates, unsigned nPatterns) {
	_nTaxa     = nTaxa;
	_nCat      = nCat;
	_nStates   = nStates;
	_nPatterns = nPatterns;
	
	// initialize the instance details
	//
	BeagleInstanceDetails instDetails;
	
    // create instances of the BEAGLE library, GPU is preferred; only one instance is created, used for the site loglikelihood
	// calculation site-by-site since each site has its own omega
	//
	// This function creates a single instance of the BEAGLE library and can be called
	// multiple times to create multiple data partition instances each returning a unique
	// identifier.
	//
	_instance = beagleCreateInstance(
									(int)nTaxa,					// Number of tip data elements
									(int)(2*nTaxa-2 - nTaxa),	// Number of partials buffers to create
									(int)nTaxa,					// Number of compact state representation buffers to create
									nStates,					// Number of states in the continuous-time Markov chain
									(int)nPatterns,				// Number of site patterns to be handled by the instance;
									(int)nPatterns,				// Number of rate matrix eigen-decomposition, category weight, and state frequency buffers to allocate
									(int)((2*nTaxa-3)*1),		// Number of transition probability matrix buffers
									(int)nCat,					// Number of rate categories; nCat=1 for the codon model
									0,							// Number of scale buffers to create, ignored for auto scale or always scale                
									NULL,						// List of potential resources on which this instance is allowed; NULL implies no restriction
									0,							// Length of resourceList list
									0,							// Bit-flags indicating preferred implementation charactertistics
									BEAGLE_FLAG_PROCESSOR_GPU | BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_EIGEN_REAL | BEAGLE_FLAG_SCALING_ALWAYS,
									// Bit-flags indicating required implementation characteristics
									&instDetails);				// Pointer to return implementation and resource details
	if(_instance < 0) {
		std::cout << "Failed to obtain beagle instance.\n";
		exit(1);
	}
	
	// list the resource being used
	//
	std::cout << "Instance " << _instance << " using resource " << instDetails.resourceNumber << ":\n";
	std::cout << "\tRsrc Name : " << instDetails.resourceName << '\n';
	std::cout << "\tImpl Name : " << instDetails.implName << '\n';
	std::cout << "\tImpl Desc : " << instDetails.implDescription << '\n';
}

void BeagleLib::SetStateFrequencies(const std::vector<double> &freqs) {
	// This function copies a state frequency array into an instance buffer.
	// maybe need to set seqLens times, same as eigenBufferCount
	//
	int code = beagleSetStateFrequencies(
										 _instance,		// Instance number
										 0,				// Index of state frequencies buffer. eigenIndex
										 &freqs[0]);	// State frequencies array (stateCount)
	if(code != 0) {
		std::cout << "Failed to set state frequencies.\n";
		exit(1);
	}
}

void BeagleLib::SetTipStates() {
#if 0 
	for(int i = 0; i < nTaxa; ++i) {
		int code;
		code = beagleSetTipStates(
								  _instance,							// Instance number
								  i,								// Index of destination compactBuffer
								  &((digitSeqVec[i].second)[0]));	// Pointer to compact states
		if(code != 0) {
			std::cout << "Failed to set tip states.\n";
			exit(1);
		}
	}
#endif
}

void BeagleLib::SetCategoryRatesAndWeights(const std::vector<double> &rates, const std::vector<double> &weights) {
	// set up Gamma rates and probabilities for rate heterogeneity
	//
	int code;
	
	// This function sets the vector of category rates for an instance.
	//
	code = beagleSetCategoryRates(
								  _instance,	// Instance number
								  &rates[0]);	// Array containing categoryCount rate scalers
	if(code != 0) {
		std::cout << "Failed to set category rates.\n";
		exit(1);
	}
	
	// This function copies a category weights array into an instance buffer.
	// maybe need to set seqLens times, same as eigenBufferCount
	//
	code = beagleSetCategoryWeights(
									_instance,	// Instance number 
									0,			// Index of category weights buffer. eigenIndex
									&weights[0]);	// Category weights array (categoryCount)
	if(code != 0) {
		std::cout << "Failed to set category weights.\n";
		exit(1);
	}
}

void BeagleLib::SetPatternWeights(const std::vector<double> &patternWeights) {
	// This function sets the vector of pattern weights for an instance.
	//
	int code = beagleSetPatternWeights(
									   _instance,			// Instance number
									   &patternWeights[0]);	// Array containing patternCount weights
	if(code != 0) {
		std::cout << "Failed to set pattern weights.\n";
		exit(1);
	}
}

void BeagleLib::SetEigenDecomposition(const std::vector<double> &eigenValues, const std::vector<double> &eigenVectors, const std::vector<double> &inverseEigenVectors) {
	// set up an eigen-decomposition buffer
	// This function copies an eigen-decomposition into an instance buffer. beaglelib transoposes the input 
	// eigenvectors and inverse-eigenvectors, which is kind of waste in my case.
	//
	for (unsigned i = 0; i < _nPatterns; ++i) {
		int code = beagleSetEigenDecomposition(
											   _instance,					// Instance number
											   (int)i,				// Index of eigen-decomposition buffer
											   &eigenVectors[0],	// Flattened matrix (stateCount x stateCount) of eigen-vectors
											   &inverseEigenVectors[0],	// Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors
											   &eigenValues[0]);			// Vector of eigenvalues
		if(code != 0) {
			std::cout << "Failed to set eigendecomposition.\n";
			exit(1);
		}
	}
}

