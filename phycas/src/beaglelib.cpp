/*
 *  beaglelib.cpp
 *  phycas
 *
 *  Created by Daniel on 2/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */	
#include "beaglelib.hpp"
#include "libhmsbeagle/beagle.h"
#include "states_patterns.hpp"
#include "tip_data.hpp"
#include <boost/format.hpp>

using namespace phycas;

BeagleLib::BeagleLib():_instance(0), _nTaxa(0), _nCat(0), _nStates(0), _nPatterns(0), _logLikelihood(0.0), _reCalLogLikelihoodNeeded(true) {
}

BeagleLib::~BeagleLib() {
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
    for(int i = 0; i < rsrcList->length; ++i) {
		std::cout << "\tResource " << i << ":\n\t\tName : " << rsrcList->list[i].name << '\n';
		std::cout << "\t\tDesc : " << rsrcList->list[i].description << '\n';
    }    
}

//void BeagleLib::Init(unsigned nTaxa, unsigned nCat, unsigned nStates, unsigned nPatterns, unsigned debug) {
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
//	if (debug==0) {
		_instance = beagleCreateInstance(
										 _nTaxa,						// Number of tip data elements
										 (2*_nTaxa-2 - _nTaxa),		// Number of partials buffers to create
										 _nTaxa,						// Number of compact state representation buffers to create
										 _nStates,					// Number of states in the continuous-time Markov chain
										 _nPatterns,					// Number of site patterns to be handled by the instance;
										 1,							// Number of rate matrix eigen-decomposition, category weight, and state frequency buffers to allocate
										 (2*_nTaxa-3),				// Number of transition probability matrix buffers
										 _nCat,						// Number of rate categories; nCat=1 for the codon model
										 0,							// Number of scale buffers to create, ignored for auto scale or always scale                
										 NULL,						// List of potential resources on which this instance is allowed; NULL implies no restriction
										 0,							// Length of resourceList list
										 0,							// Bit-flags indicating preferred implementation charactertistics
#if 1
										 BEAGLE_FLAG_PROCESSOR_GPU | 
#else
										 BEAGLE_FLAG_PROCESSOR_CPU | 
#endif
										 BEAGLE_FLAG_PRECISION_SINGLE | 
										 BEAGLE_FLAG_EIGEN_REAL | 
										 BEAGLE_FLAG_SCALING_ALWAYS,	// Bit-flags indicating required implementation characteristics
										 &instDetails);				// Pointer to return implementation and resource details
//	}

//	if (debug==1) {
//		_instance = beagleCreateInstance(
//										 _nTaxa,						// Number of tip data elements
//										 (2*_nTaxa-2 - _nTaxa),		// Number of partials buffers to create
//										 _nTaxa,						// Number of compact state representation buffers to create
//										 _nStates,					// Number of states in the continuous-time Markov chain
//										 _nPatterns,					// Number of site patterns to be handled by the instance;
//										 1,							// Number of rate matrix eigen-decomposition, category weight, and state frequency buffers to allocate
//										 (2*_nTaxa-3),				// Number of transition probability matrix buffers
//										 _nCat,						// Number of rate categories; nCat=1 for the codon model
//										 0,							// Number of scale buffers to create, ignored for auto scale or always scale                
//										 NULL,						// List of potential resources on which this instance is allowed; NULL implies no restriction
//										 0,							// Length of resourceList list
//										 0,							// Bit-flags indicating preferred implementation charactertistics
//										 BEAGLE_FLAG_PROCESSOR_CPU | 
//										 BEAGLE_FLAG_PRECISION_SINGLE | 
//										 BEAGLE_FLAG_EIGEN_REAL | 
//										 BEAGLE_FLAG_SCALING_ALWAYS,	// Bit-flags indicating required implementation characteristics
//										 &instDetails);				// Pointer to return implementation and resource details
//	}
	
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

void BeagleLib::IndexNodes(TreeShPtr t) {
	TreeNode* nd      = t->GetLastPreorder();
	int internalIndex = _nTaxa;
	int tipIndex      = 0;
	
	for (; !nd->IsTipRoot(); nd = nd->GetNextPostorder()) {
		if (nd->IsTip()) {
			nd->SetTmp((double)tipIndex++);
		}
		else {
			PHYCAS_ASSERT(nd->IsInternal());			
			nd->SetTmp((double)internalIndex++);
		}
	}	
	t->GetRoot()->SetTmp((double)tipIndex);
}

void BeagleLib::SetTipStates(TreeShPtr t, unsigned whichSubset) {	
	IndexNodes(t);
	for (preorder_iterator node = t->begin(); node != t->end(); ++node) {
		if (node->IsTip()) {
			TipData*            td = node->GetTipData();
			const state_code_t* sc = td->getConstStateCodes(whichSubset);
			std::vector<int> v(_nPatterns, _nStates);
			for (int i = 0; i < _nPatterns; ++i) {
				if (sc[i] > _nStates-1) {
					v[i] = _nStates;
				}
				else {
					v[i] = sc[i];
 				}
				int code = beagleSetTipStates(
											  _instance,			// Instance number
											  (int)node->GetTmp(),	// Index of destination compactBuffer
											  &v[0]);				// Pointer to compact states
				if (code != 0) {
					std::cout << "Failed to set tip states.\n";
					exit(1);
				}				
			}
		}
	}
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
									_instance,		// Instance number 
									0,				// Index of category weights buffer. eigenIndex
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
	int code = beagleSetEigenDecomposition(
										   _instance,				// Instance number
										   0,						// Index of eigen-decomposition buffer
										   &eigenVectors[0],		// Flattened matrix (stateCount x stateCount) of eigen-vectors
										   &inverseEigenVectors[0],	// Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors
										   &eigenValues[0]);		// Vector of eigenvalues
	if(code != 0) {
		std::cout << "Failed to set eigendecomposition.\n";
		exit(1);
	}
}

void BeagleLib::DefineOperations(TreeShPtr t, double edgelenScaler) {
	// define operation, scale index
	//  * Operations list is a list of 7-tuple integer indices, with one 7-tuple per operation.
	//  * Format of 7-tuple operation: {destinationPartials,
	//	*                               destinationScaleWrite,
	//	*                               destinationScaleRead,
	//	*                               child1Partials,
	//	*                               child1TransitionMatrix,
	//	*                               child2Partials,
	//	*                               child2TransitionMatrix}
	//
	
	_operations.clear();
	_pMatrixIndex.clear();
	_brLens.clear();
	int internalIndex = _nTaxa;
	TreeNode* nd      = t->GetLastPreorder();
	
	for (; !nd->IsTipRoot(); nd = nd->GetNextPostorder()) {
		if (nd->IsTip()) {
			_pMatrixIndex.push_back((int)nd->GetTmp());
			_brLens.push_back((nd->GetEdgeLen())*edgelenScaler);
		}
		else {
			PHYCAS_ASSERT(nd->IsInternal());
			
			// destination partial to be calculated
			//			
			_operations.push_back(internalIndex);
			_pMatrixIndex.push_back(internalIndex);
			_brLens.push_back((nd->GetEdgeLen())*edgelenScaler);
			nd->SetTmp(internalIndex++);
			
			// destination scaling buffer index to write to
			//
			_operations.push_back(BEAGLE_OP_NONE);
			
			// destination scaling buffer index to read from
			//
			_operations.push_back(BEAGLE_OP_NONE);
			
			// left child partial index
			//
			int leftChildIndex = (int)nd->GetLeftChild()->GetTmp();
			_operations.push_back(leftChildIndex);
			
			// left child transition matrix index
			//
			_operations.push_back(leftChildIndex);
			
			// right child partial index
			//
			int rightChildIndex = (int)nd->GetLeftChild()->GetRightSib()->GetTmp();
			_operations.push_back(rightChildIndex); // TODO assumes binary tree
			
			// right child transition matrix index
			//
			_operations.push_back(rightChildIndex);
		}
	}
	
	_pMatrixIndex[_pMatrixIndex.size()-1] = (int)t->GetRoot()->GetTmp();
}

double BeagleLib::CalcLogLikelihood(TreeShPtr t) {
	if (!_reCalLogLikelihoodNeeded) {
		return _logLikelihood;
	}
	else {
		int list_len = (int)(2*_nTaxa-3);
		int code;
		
		// update transition probability matrix
		// Calculate a list of transition probability matrices
		// 
		code = beagleUpdateTransitionMatrices(
											  _instance,			// Instance number
											  0,					// Index of eigen-decomposition buffer
											  &_pMatrixIndex[0],	// List of indices of transition probability matrices to update
											  NULL,					// List of indices of first derivative matrices to update (NULL implies no calculation)
											  NULL,					// List of indices of second derivative matrices to update (NULL implies no calculation)
											  &_brLens[0],			// List of edge lengths with which to perform calculations
											  list_len);			// Length of lists
		if(code != 0) {
			std::cout << "Failed to update transition matrices.\n";
			exit(1);
		}
		

//		//debug
//		for(int i = 0; i < (2*_nTaxa-3); ++i) {
//			std::vector<double> outMatrix(4*4, 0.0);
//			beagleGetTransitionMatrix(_instance, i, &outMatrix[0]);
//			std::cerr << "transition matrix " << i << ":\n";
//			for(unsigned i = 0; i < 4; ++i) {
//				for(unsigned j = 0; j < 4; ++j)
//					std::cerr << outMatrix[i*4+j] << '\t';
//				std::cerr << '\n';
//			}	
//		}
//		char ch;
//		std::cin >> ch;
		
		
		// Calculate or queue for calculation partials using a list of operations
		//
		int totalOperations = (int)(_operations.size()/7);
		code = beagleUpdatePartials(
									_instance,							// Instance number
									(BeagleOperation*)&_operations[0],	// List of 7-tuples specifying operations
									totalOperations,					// Number of operations
									BEAGLE_OP_NONE);					// Index number of scaleBuffer to store accumulated factors
		if(code != 0) {
			std::cout << "Failed to update partials.\n";
			exit(1);
		}
		
		// This function integrates a list of partials at a parent and child node with respect
		// to a set of partials-weights and state frequencies to return the log likelihood
		// and first and second derivative sums
		//
		int stateFrequencyIndex  = 0;
		int categoryWeightsIndex = 0;	
		int cumulativeScalingIndex = BEAGLE_OP_NONE;
		int childIndex = (int)t->GetRoot()->GetTmp();
		int parentIndex = (int)t->GetRoot()->GetLeftChild()->GetTmp();
		_logLikelihood = 0.0;
		code = beagleCalculateEdgeLogLikelihoods(
												 _instance,					// Instance number
												 &parentIndex,				// List of indices of parent partialsBuffers
												 &childIndex,				// List of indices of child partialsBuffers
												 &childIndex,				// List indices of transition probability matrices for this edge
												 NULL,						// List indices of first derivative matrices
												 NULL,						// List indices of second derivative matrices
												 &categoryWeightsIndex,		// List of weights to apply to each partialsBuffer
												 &stateFrequencyIndex,		// List of state frequencies for each partialsBuffer. There should be one set for each of parentBufferIndices.
												 &cumulativeScalingIndex,	// List of scaleBuffers containing accumulated factors to apply to each partialsBuffer. There should be one index for each of parentBufferIndices.
												 1,							// Number of partialsBuffers
												 &_logLikelihood,			// Pointer to destination for resulting log likelihood
												 NULL,						// Pointer to destination for resulting first derivative
												 NULL);						// Pointer to destination for resulting second derivative
		if(code != 0) {
			std::cout << "Failed to calculate edge logLikelihoods in CalcLogLikelihood. The error code is " << code << ".\n";
			exit(1);
		}
		
		return _logLikelihood;
	}
}
