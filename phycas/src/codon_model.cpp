/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <cmath>
#include <iostream>
#include "ncl/nxsallocatematrix.h"
#include "phycas/src/likelihood_models.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

static char * codon_code_to_triplet[] =
	{
	"AAA",
	"AAC",
	"AAG",
	"AAT",
	"ACA",
	"ACC",
	"ACG",
	"ACT",
	"AGA",
	"AGC",
	"AGG",
	"AGT",
	"ATA",
	"ATC",
	"ATG",
	"ATT",
	"CAA",
	"CAC",
	"CAG",
	"CAT",
	"CCA",
	"CCC",
	"CCG",
	"CCT",
	"CGA",
	"CGC",
	"CGG",
	"CGT",
	"CTA",
	"CTC",
	"CTG",
	"CTT",
	"GAA",
    "GAC",
	"GAG",
	"GAT",
	"GCA",
	"GCC",
	"GCG",
	"GCT",
	"GGA",
	"GGC",
	"GGG",
	"GGT",
	"GTA",
	"GTC",
	"GTG",
	"GTT",
	"TAC",
	"TAT",
	"TCA",
	"TCC",
	"TCG",
	"TCT",
	"TGC",
	"TGG",
	"TGT",
	"TTA",
	"TTC",
	"TTG",
	"TTT"
	};

#if 0
static char * codon_code_to_aa_name[] =
	{
	"Lys", //  0 AAA
	"Asn", //  1 AAC
	"Lys", //  2 AAG
	"Asn", //  3 AAT
	"Thr", //  4 ACA
	"Thr", //  5 ACC
	"Thr", //  6 ACG
	"Thr", //  7 ACT
	"Arg", //  8 AGA
	"Ser", //  9 AGC
	"Arg", // 10 AGG
	"Ser", // 11 AGT
	"Ile", // 12 ATA
	"Ile", // 13 ATC
	"Met", // 14 ATG
	"Ile", // 15 ATT
	"Gln", // 16 CAA
	"His", // 17 CAC
	"Gln", // 18 CAG
	"His", // 19 CAT
	"Pro", // 20 CCA
	"Pro", // 21 CCC
	"Pro", // 22 CCG
	"Pro", // 23 CCT
	"Arg", // 24 CGA
	"Arg", // 25 CGC
	"Arg", // 26 CGG
	"Arg", // 27 CGT
	"Leu", // 28 CTA
	"Leu", // 29 CTC
	"Leu", // 30 CTG
	"Leu", // 31 CTT
	"Glu", // 32 GAA
	"Asp", // 33 GAC
	"Glu", // 34 GAG
	"Asp", // 35 GAT
	"Ala", // 36 GCA
	"Ala", // 37 GCC
	"Ala", // 38 GCG
	"Ala", // 39 GCT
	"Gly", // 40 GGA
	"Gly", // 41 GGC
	"Gly", // 42 GGG
	"Gly", // 43 GGT
	"Val", // 44 GTA
	"Val", // 45 GTC
	"Val", // 46 GTG
	"Val", // 47 GTT
	"Tyr", // 48 TAC
	"Tyr", // 49 TAT
	"Ser", // 50 TCA
	"Ser", // 51 TCC
	"Ser", // 52 TCG
	"Ser", // 53 TCT
	"Cys", // 54 TGC
	"Trp", // 55 TGG
	"Cys", // 56 TGT
	"Leu", // 57 TTA
	"Phe", // 58 TTC
	"Leu", // 59 TTG
	"Phe"  // 60 TTT
	};
#endif

//    0    Ala     A       Alanine
//    1    Arg     R       Arginine
//    2    Asn     N       Asparagine
//    3    Asp     D       Aspartic acid (Aspartate)
//    4    Cys     C       Cysteine
//    5    Gln     Q       Glutamine
//    6    Glu     E       Glutamic acid (Glutamate)
//    7    Gly     G       Glycine
//    8    His     H       Histidine
//    9    Ile     I       Isoleucine
//   10    Leu     L       Leucine
//   11    Lys     K       Lysine
//   12    Met     M       Methionine
//   13    Phe     F       Phenylalanine
//   14    Pro     P       Proline
//   15    Ser     S       Serine
//   16    Thr     T       Threonine
//   17    Trp     W       Tryptophan
//   18    Tyr     Y       Tyrosine
//   19    Val     V       Valine

static unsigned codon_code_to_aa_code[] =
	{
	11, //  0 AAA "Lys"
	2, //  1 AAC "Asn"
	11, //  2 AAG "Lys"
	2, //  3 AAT "Asn"
	16, //  4 ACA "Thr"
	16, //  5 ACC "Thr"
	16, //  6 ACG "Thr"
	16, //  7 ACT "Thr"
	1, //  8 AGA "Arg"
	15, //  9 AGC "Ser"
	1, // 10 AGG "Arg"
	15, // 11 AGT "Ser"
	9, // 12 ATA "Ile"
	9, // 13 ATC "Ile"
	12, // 14 ATG "Met"
	9, // 15 ATT "Ile"
	5, // 16 CAA "Gln"
	8, // 17 CAC "His"
	5, // 18 CAG "Gln"
	8, // 19 CAT "His"
	14, // 20 CCA "Pro"
	14, // 21 CCC "Pro"
	14, // 22 CCG "Pro"
	14, // 23 CCT "Pro"
	1, // 24 CGA "Arg"
	1, // 25 CGC "Arg"
	1, // 26 CGG "Arg"
	1, // 27 CGT "Arg"
	10, // 28 CTA "Leu"
	10, // 29 CTC "Leu"
	10, // 30 CTG "Leu"
	10, // 31 CTT "Leu"
	6, // 32 GAA "Glu"
    3, // 33 GAC "Asp"
	6, // 34 GAG "Glu"
	3, // 35 GAT "Asp"
	0, // 36 GCA "Ala"
	0, // 37 GCC "Ala"
	0, // 38 GCG "Ala"
	0, // 39 GCT "Ala"
	7, // 40 GGA "Gly"
	7, // 41 GGC "Gly"
	7, // 42 GGG "Gly"
	7, // 43 GGT "Gly"
	19, // 44 GTA "Val"
	19, // 45 GTC "Val"
	19, // 46 GTG "Val"
	19, // 47 GTT "Val"
	18, // 48 TAC "Tyr"
	18, // 49 TAT "Tyr"
	15, // 50 TCA "Ser"
	15, // 51 TCC "Ser"
	15, // 52 TCG "Ser"
	15, // 53 TCT "Ser"
	4, // 54 TGC "Cys"
	17, // 55 TGG "Trp"
	4, // 56 TGT "Cys"
	10, // 57 TTA "Leu"
	13, // 58 TTC "Phe"
	10, // 59 TTG "Leu"
	13  // 60 TTT "Phe"
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 61, the codon frequencies to 1/61, and `kappa' and `omega' both to 1.0.
*/
Codon::Codon()
  : Model(61), kappa(1.0), omega(1.0), kappa_fixed(false), omega_fixed(false)
	{
	is_codon_model = true;
	state_repr.reserve(64);
	char * bases[] = {"A", "C", "G", "T"};
	for (unsigned i = 0; i < 4; ++i)
		{
		for (unsigned j = 0; j < 4; ++j)
			{
			for (unsigned k = 0; k < 4; ++k)
				{
				std::string s = str(boost::format("%c%c%c") % bases[i] % bases[j] % bases[k]);
				state_repr.push_back(s);
				}
			}
		}
	// ignore stop codons
	state_repr.erase(state_repr.begin()+48);
	state_repr.erase(state_repr.begin()+49);
	state_repr.erase(state_repr.begin()+54);

	updateQMatrix();

	//std::ofstream outf("debug_codon_identities.txt");
	//outf << "Codon identities:\n\n";
	//unsigned z = 0;
	//for (std::vector<std::string>::const_iterator it = state_repr.begin(); it != state_repr.end(); ++it)
	//	{
	//	outf << (z++) << '\t' << (*it) << '\n'; 
	//	}
	//outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "HKY85", "HKY85+G", "HKY85+I" or "HKY85+G+I".
*/
std::string Codon::getModelName() const
	{
	std::string s = "Codon";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any 
|	parameters related to rate heterogeneity. This function then adds additional codon-model-specific parameters to the 
|	supplied `parameters' vector. This incudes the codon frequencies as well as the transition/transversion rate 
|	ratio kappa and the nonsynonymous/synonymous rate ratio omega.
*/
void Codon::createParameters(
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the vector of edge length hyperparameters to fill */
  MCMCUpdaterVect & parameters) const		/**< is the vector of model-specific parameters to fill */
	{
	Model::createParameters(t, edgelens, edgelen_hyperparams, parameters);

	PHYCAS_ASSERT(!kappa_param);
	kappa_param = MCMCUpdaterShPtr(new KappaParam());
	kappa_param->setName("trs/trv rate ratio");
	kappa_param->setStartingValue(4.0);
	kappa_param->setTree(t);
	kappa_param->setPrior(kappa_prior);
	if (kappa_fixed)
		kappa_param->fixParameter();
	parameters.push_back(kappa_param);

	PHYCAS_ASSERT(!omega_param);
	omega_param = MCMCUpdaterShPtr(new OmegaParam());
	omega_param->setName("nonsynon./synon. rate ratio");
	omega_param->setStartingValue(1.0);
	omega_param->setTree(t);
	omega_param->setPrior(omega_prior);
	if (omega_fixed)
		omega_param->fixParameter();
	parameters.push_back(omega_param);

	PHYCAS_ASSERT(freq_params.empty());

	for (unsigned i = 0; i < 61; ++i)
		{
		MCMCUpdaterShPtr state_freq_param = MCMCUpdaterShPtr(new StateFreqParam(i));
		std::string s = str(boost::format("freq. for codon %s") % state_repr[i]);
		state_freq_param->setName(s);
		state_freq_param->setTree(t);
		state_freq_param->setStartingValue(1.0);
		state_freq_param->setPrior(freq_param_prior);
		if (state_freq_fixed)
			state_freq_param->fixParameter();
		parameters.push_back(state_freq_param);
		freq_params.push_back(state_freq_param);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by 
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for 
|	compatability with MrBayes). The Codon model provide additional columns for kappa, omega, the codon frequencies, 
|	the gamma shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if an invariable 
|	sites model is being used)
*/
std::string Codon::paramHeader() const
	{
	std::string s = std::string("Gen\tLnL\tTL\tkappa\tomega\t");
	for (std::vector<std::string>::const_iterator it = state_repr.begin(); it != state_repr.end(); ++it)
		{
		s += "freq";
		s += (*it); 
		s += '\t';
		}
	if (is_flex_model)
		{
		s += "\tncat";
		}
	else if (num_gamma_rates > 1)
		{
		s += "\tshape";
		}
	if (is_pinvar_model)
		s += "\tpinvar";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific 
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes).
*/
std::string Codon::paramReport(
  unsigned ndecimals) const /**< floating point precision to use */
	{
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = boost::str(boost::format(fmt) % kappa);
	s += boost::str(boost::format(fmt) % omega);
	//std::string s = boost::str(boost::format("\t%.5f\t%.5f\t") % kappa % omega);
	for (unsigned i = 0; i < 61; ++i)
		{
		s += str(boost::format(fmt) % state_freqs[i]);
		//s += str(boost::format("%.5f\t") % state_freqs[i]);
		}
	if (is_flex_model)
		{
		s += str(boost::format("%d\t") % num_gamma_rates);
		//s += str(boost::format("\t%d") % num_gamma_rates);
		}
	else if (num_gamma_rates > 1)
		{
		s += str(boost::format(fmt) % gamma_shape);
		//s += str(boost::format("\t%.5f") % gamma_shape);
		}
	if (is_pinvar_model)
        {
		s += str(boost::format(fmt) % pinvar);
		//s += str(boost::format("\t%.5f") % pinvar);
        }
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to true. The fixParameter member function of the KappaParam object is either 
|	called immediately (if `kappa_param' is a valid pointer) or is called in createParameters (when `kappa_param' is 
|	first assigned).
*/
void Codon::fixKappa()
	{
	kappa_fixed = true;
	if (kappa_param)
		kappa_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to false. The freeParameter member function of the KappaParam object is called 
|	immediately if `kappa_param' is a valid pointer.
*/
void Codon::freeKappa()
	{
	kappa_fixed = false;
	if (kappa_param)
		kappa_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `omega_fixed' to true. The fixParameter member function of the OmegaParam object is either 
|	called immediately (if `omega_param' is a valid pointer) or is called in createParameters (when `omega_param' is 
|	first assigned).
*/
void Codon::fixOmega()
	{
	omega_fixed = true;
	if (omega_param)
		omega_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `omega_fixed' to false. The freeParameter member function of the OmegaParam object is called 
|	immediately if `omega_param' is a valid pointer.
*/
void Codon::freeOmega()
	{
	omega_fixed = false;
	if (omega_param)
		omega_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa'.
*/
double Codon::getKappa()
 	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' to supplied value `k'. Throws XLikelihood exception if `k' is less than or equal to 0.0.
*/
void Codon::setKappa(double k)
 	{
	if (k <= 0.0)
		throw XLikelihood();
	kappa = k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `omega'.
*/
double Codon::getOmega()
 	{
	return omega;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `omega' to supplied value `w'. Throws XLikelihood exception if `w' is less than or equal to 0.0.
*/
void Codon::setOmega(double w)
 	{
	if (w <= 0.0)
		throw XLikelihood();
	omega = w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa_prior'.
*/
ProbDistShPtr Codon::getKappaPrior()
 	{
	return kappa_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setKappaPrior(ProbDistShPtr d)
 	{
	kappa_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `omega_prior'.
*/
ProbDistShPtr Codon::getOmegaPrior()
 	{
	return omega_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `omega_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setOmegaPrior(ProbDistShPtr d)
 	{
	omega_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr Codon::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from 
|	the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
void Codon::calcPMat(double * * pMat, double edgeLength) const
	{
	updateQMatrix();
	q_matrix.recalcPMat(pMat, edgeLength);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Codon::calcUniformizationLambda() const
	{
    assert(0);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double Codon::calcLMat(double * * lMat) const
	{
    std::cerr << "Error in Codon::calcLMat: q_matrix does not yet have the required recalcLMat function" << std::endl;
    assert(0);
	//updateQMatrix();
	//q_matrix.recalcLMat(lMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function 
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double Codon::calcUMat(double * * uMat) const
	{
    std::cerr << "Error in Codon::calcUMat: q_matrix does not yet have the required recalcUMat function" << std::endl;
    assert(0);
	//updateQMatrix();
	//q_matrix.recalcUMat(uMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all state frequencies to 1/num_states. The base class version is 
|	called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member 
|	knows about the change in state frequencies.
*/
void Codon::setAllFreqsEqual()
	{
	Model::setAllFreqsEqual();
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the 61 codon state frequencies to the values expected based on the four base frequencies provided (all four 
|	values provided should greater than or equal to 0.0, but do not need to sum to 1.0). The frequency of a codon is, 
|	almost, the product of the three component base frequencies: the "almost" qualification being needed because the
|	three stop codons are not included, so each codon frequency must be corrected by dividing by the sum of the 61 
|	non-stop three-nucleotide products. For example, if the specified base freqencies were 0.1, 0.2, 0.3 and 0.4, then 
|	the frequency of the ACT codon would be the product (0.1)*(0.2)*(0.4) = 0.008, divided by the sum of all 61 such 
|	products, which in this case is 0.972, yielding 0.008/0.971 = 0.00823. Note state frequencies set using this 
|	function will be obliterated unless state frequencies are fixed using the fixStateFreqs method.
*/
inline void Codon::setNucleotideFreqs(
  double freqA,				/**< the new value of `state_freq_unnorm'[0] (i.e. frequency of base A) */
  double freqC,				/**< the new value of `state_freq_unnorm'[1] (i.e. frequency of base C) */
  double freqG,				/**< the new value of `state_freq_unnorm'[2] (i.e. frequency of base G) */
  double freqT)				/**< the new value of `state_freq_unnorm'[3] (i.e. frequency of base T/U) */
	{
	PHYCAS_ASSERT(num_states == 4);
	PHYCAS_ASSERT(freqA >= 0.0);
	PHYCAS_ASSERT(freqC >= 0.0);
	PHYCAS_ASSERT(freqG >= 0.0);
	PHYCAS_ASSERT(freqT >= 0.0);
	PHYCAS_ASSERT(state_freq_unnorm.size() == 61);
	state_freq_unnorm[0]  = freqA*freqA*freqA;	// 0 AAA
	state_freq_unnorm[1]  = freqA*freqA*freqC;	// 1 AAC
	state_freq_unnorm[2]  = freqA*freqA*freqG;	// 2 AAG
	state_freq_unnorm[3]  = freqA*freqA*freqT;	// 3 AAT
	state_freq_unnorm[4]  = freqA*freqC*freqA;	// 4 ACA
	state_freq_unnorm[5]  = freqA*freqC*freqC;	// 5 ACC
	state_freq_unnorm[6]  = freqA*freqC*freqG;	// 6 ACG
	state_freq_unnorm[7]  = freqA*freqC*freqT;	// 7 ACT
	state_freq_unnorm[8]  = freqA*freqG*freqA;	// 8 AGA
	state_freq_unnorm[9]  = freqA*freqG*freqC;	// 9 AGC
	state_freq_unnorm[10] = freqA*freqG*freqG;	// 10 AGG
	state_freq_unnorm[11] = freqA*freqG*freqT;	// 11 AGT
	state_freq_unnorm[12] = freqA*freqT*freqA;	// 12 ATA
	state_freq_unnorm[13] = freqA*freqT*freqC;	// 13 ATC
	state_freq_unnorm[14] = freqA*freqT*freqG;	// 14 ATG
	state_freq_unnorm[15] = freqA*freqT*freqT;	// 15 ATT
	state_freq_unnorm[16] = freqC*freqA*freqA;	// 16 CAA
	state_freq_unnorm[17] = freqC*freqA*freqC;	// 17 CAC
	state_freq_unnorm[18] = freqC*freqA*freqG;	// 18 CAG
	state_freq_unnorm[19] = freqC*freqA*freqT;	// 19 CAT
	state_freq_unnorm[20] = freqC*freqC*freqA;	// 20 CCA
	state_freq_unnorm[21] = freqC*freqC*freqC;	// 21 CCC
	state_freq_unnorm[22] = freqC*freqC*freqG;	// 22 CCG
	state_freq_unnorm[23] = freqC*freqC*freqT;	// 23 CCT
	state_freq_unnorm[24] = freqC*freqG*freqA;	// 24 CGA
	state_freq_unnorm[25] = freqC*freqG*freqC;	// 25 CGC
	state_freq_unnorm[26] = freqC*freqG*freqG;	// 26 CGG
	state_freq_unnorm[27] = freqC*freqG*freqT;	// 27 CGT
	state_freq_unnorm[28] = freqC*freqT*freqA;	// 28 CTA
	state_freq_unnorm[29] = freqC*freqT*freqC;	// 29 CTC
	state_freq_unnorm[30] = freqC*freqT*freqG;	// 30 CTG
	state_freq_unnorm[31] = freqC*freqT*freqT;	// 31 CTT
	state_freq_unnorm[32] = freqG*freqA*freqA;	// 32 GAA
	state_freq_unnorm[33] = freqG*freqA*freqC;	// 33 GAC
	state_freq_unnorm[34] = freqG*freqA*freqG;	// 34 GAG
	state_freq_unnorm[35] = freqG*freqA*freqT;	// 35 GAT
	state_freq_unnorm[36] = freqG*freqC*freqA;	// 36 GCA
	state_freq_unnorm[37] = freqG*freqC*freqC;	// 37 GCC
	state_freq_unnorm[38] = freqG*freqC*freqG;	// 38 GCG
	state_freq_unnorm[39] = freqG*freqC*freqT;	// 39 GCT
	state_freq_unnorm[40] = freqG*freqG*freqA;	// 40 GGA
	state_freq_unnorm[41] = freqG*freqG*freqC;	// 41 GGC
	state_freq_unnorm[42] = freqG*freqG*freqG;	// 42 GGG
	state_freq_unnorm[43] = freqG*freqG*freqT;	// 43 GGT
	state_freq_unnorm[44] = freqG*freqT*freqA;	// 44 GTA
	state_freq_unnorm[45] = freqG*freqT*freqC;	// 45 GTC
	state_freq_unnorm[46] = freqG*freqT*freqG;	// 46 GTG
	state_freq_unnorm[47] = freqG*freqT*freqT;	// 47 GTT
	state_freq_unnorm[48] = freqT*freqA*freqC;	// 48 TAC
	state_freq_unnorm[49] = freqT*freqA*freqT;	// 49 TAT
	state_freq_unnorm[50] = freqT*freqC*freqA;	// 50 TCA
	state_freq_unnorm[51] = freqT*freqC*freqC;	// 51 TCC
	state_freq_unnorm[52] = freqT*freqC*freqG;	// 52 TCG
	state_freq_unnorm[53] = freqT*freqC*freqT;	// 53 TCT
	state_freq_unnorm[54] = freqT*freqG*freqC;	// 54 TGC
	state_freq_unnorm[55] = freqT*freqG*freqG;	// 55 TGG
	state_freq_unnorm[56] = freqT*freqG*freqT;	// 56 TGT
	state_freq_unnorm[57] = freqT*freqT*freqA;	// 57 TTA
	state_freq_unnorm[58] = freqT*freqT*freqC;	// 58 TTC
	state_freq_unnorm[59] = freqT*freqT*freqG;	// 59 TTG
	state_freq_unnorm[60] = freqT*freqT*freqT;	// 60 TTT
	normalizeFreqs();
	q_matrix.setStateFreqs(state_freqs);

	//std::ofstream outf("debug_codon_frequencies.txt");
	//outf << "Codon state frequencies:\n\n";
	//unsigned z = 0;
	//for (std::vector<double>::const_iterator it = state_freqs.begin(); it != state_freqs.end(); ++it)
	//	{
	//	outf << (z++) << '\t' << (*it) << '\n'; 
	//	}
	//outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses current values of `omega' and `kappa' to update the relative rate (rr) vector inside `q_matrix'.
*/
void Codon::updateQMatrix() const
	{
	std::vector<double> rel_rate_vect;
	rel_rate_vect.reserve(1830);

	//std::ofstream outf("debug_qmatrix.txt");
	//outf << "QMatrix calculation (kappa = " << kappa << ", omega = " << omega << ")" << std::endl;
	//outf << "from_codon" << '\t' << "to_codon" << '\t' << "from_aa" << '\t' << "to_aa" << '\t' << "nchanges" << '\t';
	//outf << "type" << '\t' << "rate" << std::endl;

	for (unsigned row = 0; row < 61; ++row)
		{
		unsigned from_aa_state = codon_code_to_aa_code[row];
		for (unsigned col = row + 1; col < 61; ++col)
			{
			unsigned to_aa_state = codon_code_to_aa_code[col];
			unsigned nchanges = 0;
			char from_base = '\0';
			char to_base = '\0';
			for (unsigned k = 0; k < 3; ++k)
				{
				if (codon_code_to_triplet[row][k] != codon_code_to_triplet[col][k])
					{
					from_base = codon_code_to_triplet[row][k];
					to_base = codon_code_to_triplet[col][k];
					++nchanges;
					}
				}

			//outf << codon_code_to_triplet[row] << '\t' << codon_code_to_triplet[col] << '\t';
			//outf << codon_code_to_aa_name[row] << '\t' << codon_code_to_aa_name[col] << '\t';
			//outf << nchanges << '\t';

			double rel_rate = 0.0;
			//assert(nchanges > 0);
            PHYCAS_ASSERT(nchanges > 0);
			if (nchanges == 1)
				{
				// Determine transition vs. transversion
				bool transversion = true;
				if ((from_base == 'A' && to_base == 'G') || (from_base == 'C' && to_base == 'T') || (from_base == 'G' && to_base == 'A') || (from_base == 'T' && to_base == 'C'))
					transversion = false;

				if (from_aa_state == to_aa_state)
					{
					// synonymous substitution
					rel_rate = (transversion ? 1.0 : kappa);
					//outf << (transversion ? "SynTrv" : "SynTrs") << '\t' << rel_rate;
					}
				else
					{
					// nonsynonymous substitution
					rel_rate = (transversion ? omega : omega*kappa);
					//outf << (transversion ? "NonTrv" : "NonTrs") << '\t' << rel_rate;
					}
				}
			//else
			//	outf << "" << '\t' << rel_rate;

			//outf << std::endl;
			
			rel_rate_vect.push_back(rel_rate);
			}
		}
	q_matrix.setRelativeRates(rel_rate_vect);
	q_matrix.setStateFreqs(state_freqs);

	//outf.close();
	}
