/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/oldphycas/distribution_command_param.hpp"
using std::vector;
using std::pair;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::isalpha;
#endif

bool DistributionCmdOption::ReadBernoulli(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kBernoulli;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts BERNOULLI({PROBABILITY =} 0.3) or BERNOULLI(mean = .3)
	//
	vector<const VecDistParamDesc *> parseOptions;
	const DistParamDesc probP = DistParamDesc("Probability", kProbParam, 0, true, true);
	const DistParamDesc meanB = DistParamDesc("MEAN", kProbParam, 0, false, false);
	VecDistParamDesc preferred, meanForm;
	preferred.push_back(probP);
	meanForm.push_back(meanB);
	parseOptions.push_back(&preferred);
	parseOptions.push_back(&meanForm);
	
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting Bernoulli(Probability = <number between 0 and 1>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	return true;
	}
	
bool DistributionCmdOption::ReadBinomial(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kBinomial;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token, "(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts BINOMIAL({PROBABILITY =} 0.3, {NTRIALS =} 10) or BINOMIAL(NTRIALS = 10, PROBABILITY = 0.3) 
	//
	const DistParamDesc probFirst = DistParamDesc("Probability", kProbParam, 0, true, true);
	const DistParamDesc nTrialsSecond = DistParamDesc("Ntrials", kPosIntegerParam, 0, true, true);
	const DistParamDesc probSecond = DistParamDesc("Probability", kProbParam, 0, true, false);
	const DistParamDesc nTrialsFirst = DistParamDesc("Ntrials", kPosIntegerParam, 0, true, false);
	VecDistParamDesc preferred, backward;
	preferred.push_back(probFirst);
	preferred.push_back(nTrialsSecond);
	backward.push_back(nTrialsFirst);
	backward.push_back(probSecond);
	vector<const VecDistParamDesc *> parseOptions;
	parseOptions.push_back(&preferred);
	parseOptions.push_back(&backward);
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting Binomial(Probability = <number between 0 and 1>, NTrials = <integer>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	if (classOfDistribution & kZeroToOneBit) 
		{
		if (currentValue->discreteVariables[0] > 1)
			{
			errState |= NxsCmdOption::query_for_error;
			errSnippet = "Expecting a ";
			errSnippet << GetExpectedDistributionString() << " (the NTrials parameter for a Binomial distribution must be 1)";
			return false;
			}
		}
	return true;
	}
	
bool DistributionCmdOption::ReadExponential(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kExponential;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts Exponential({HAZARD =} 0.333333) or Exponential(Mean = 3) 
	//
	const DistParamDesc hazardParam = DistParamDesc("Hazard", kPositiveParam, 0, true, true);
	const DistParamDesc meanParam = DistParamDesc("MEAN", kPositiveParam, 0, false, false);
	VecDistParamDesc preferred, meanForm;
	preferred.push_back(hazardParam);
	meanForm.push_back(meanParam);
	vector<const VecDistParamDesc *> parseOptions;
	parseOptions.push_back(&preferred);
	parseOptions.push_back(&meanForm);
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting Exponential(Hazard = <number>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	return true;
	}
	
bool DistributionCmdOption::ReadGamma(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kGamma;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts GAMma({SHape = } 1, {SCale = } 2), GAMma(SCale = 2, SHape = 1), GAMma(MEAN = 1, STDDEV = 2),  GAMma(MEAN = 1, VAR = 2)
	//
	const DistParamDesc shapeParam = DistParamDesc("SHape", kPositiveParam, 0, true, true);
	const DistParamDesc scaleParam = DistParamDesc("SCale", kPositiveParam, 1, true, true);
	const DistParamDesc shapeBackwardParam = DistParamDesc("SHape", kPositiveParam, 0, true, false);
	const DistParamDesc scaleBackwardParam = DistParamDesc("SCale", kPositiveParam, 1, true, false);
	const DistParamDesc meanParam = DistParamDesc("MEAN", kPositiveParam, 0, false, false);
	const DistParamDesc stdDevParam = DistParamDesc("STDDEV", kNonNegativeParam, 1, false, false);
	const DistParamDesc varParam = DistParamDesc("VARiance", kNonNegativeParam, 1, false, false);
	VecDistParamDesc paramForm,  paramBackwardForm, meanSDForm, meanVarForm;
	paramForm.push_back(shapeParam);
	paramForm.push_back(scaleParam);
	paramBackwardForm.push_back(scaleBackwardParam);
	paramBackwardForm.push_back(shapeBackwardParam);
	meanSDForm.push_back(meanParam);
	meanSDForm.push_back(stdDevParam);
	meanVarForm.push_back(meanParam);
	meanVarForm.push_back(varParam);
	vector<const VecDistParamDesc *> parseOptions;
	parseOptions.push_back(&paramForm);
	parseOptions.push_back(&paramBackwardForm);
	parseOptions.push_back(&meanSDForm);
	parseOptions.push_back(&meanVarForm);
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting Gamma(Shape = <number>, Scale = <number>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	return true;	
	}

bool DistributionCmdOption::ReadInvGamma(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kInverseGamma;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts INVGAMma({SHape = } 1, {SCale = } 2), INVGAMma(SCale = 2, SHape = 1), INVGAMma(MEAN = 1, STDDEV = 2),  INVGAMma(MEAN = 1, VAR = 2)
	//
	const DistParamDesc shapeParam = DistParamDesc("SHape", kPositiveParam, 0, true, true);
	const DistParamDesc scaleParam = DistParamDesc("SCale", kPositiveParam, 1, true, true);
	const DistParamDesc shapeBackwardParam = DistParamDesc("SHape", kPositiveParam, 0, true, false);
	const DistParamDesc scaleBackwardParam = DistParamDesc("SCale", kPositiveParam, 1, true, false);
	const DistParamDesc meanParam = DistParamDesc("MEAN", kPositiveParam, 0, false, false);
	const DistParamDesc stdDevParam = DistParamDesc("STDDEV", kNonNegativeParam, 1, false, false);
	const DistParamDesc varParam = DistParamDesc("VARiance", kNonNegativeParam, 1, false, false);
	VecDistParamDesc paramForm,  paramBackwardForm, meanSDForm, meanVarForm;
	paramForm.push_back(shapeParam);
	paramForm.push_back(scaleParam);
	paramBackwardForm.push_back(scaleBackwardParam);
	paramBackwardForm.push_back(shapeBackwardParam);
	meanSDForm.push_back(meanParam);
	meanSDForm.push_back(stdDevParam);
	meanVarForm.push_back(meanParam);
	meanVarForm.push_back(varParam);
	vector<const VecDistParamDesc *> parseOptions;
	parseOptions.push_back(&paramForm);
	parseOptions.push_back(&paramBackwardForm);
	parseOptions.push_back(&meanSDForm);
	parseOptions.push_back(&meanVarForm);
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting InverseGamma(Shape = <number>, Scale = <number>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	if (currentValue->meaningOfVariable == DistributionDescription::kParam)
		{
		PHYCAS_ASSERT(currentValue->contVariables.size() == 2);
		if (currentValue->contVariables[0] <= 2)
			{
			errState = query_for_error;
			errSnippet = "the shape parameter of the InverseGamma distribution must be greater than two";
			return false;
			}
		}
	return true;
	}
	
bool DistributionCmdOption::ReadDirichlet(
  NxsToken &token,	/* pointer to token that should hold the next value of the setting in string form */
  bool			asBeta)
	{
	currentValue->family = DistributionDescription::kDirichlet;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	if (nVariates != 0)
		{
		//	accepts Dirichlet(0.333333 , # , #, ...)
		//
		
		//	we'll accept dirichlet's for single parameter distributions (e.g. you might want to put a beta on pinv)
		//	if we're reading a single variate distribution, expect 2 dirichlet parameters
		//
		unsigned nDirParams = (nVariates == 1 ? 2 : nVariates);
		VecDistParamDesc dirForm;
		for (unsigned i = 0; i < nDirParams; ++i)
			dirForm.push_back(DistParamDesc(string(), kPositiveParam, i, true, true));
		vector<const VecDistParamDesc *> parseOptions;
		parseOptions.push_back(&dirForm);
		if (ReadParams(parseOptions, token) == NULL)
			{
			errState = query_for_error;
			errSnippet.clear();
			errSnippet << "Expecting ";
			if (asBeta)
				errSnippet << "Beta";
			else
				errSnippet << "Dirichlet";
			errSnippet << "(<number>";
			for (unsigned j = 0; j < nVariates; ++j)
				 errSnippet << ", <number>";
			errSnippet << ") but encountered " << distCmdOptcmdRead << "...";
			return false;
			}
		}
	else
		{
		PHYCAS_ASSERT(!asBeta);
		// nVariates unknown read any number or params (> 1)
		//
		unsigned nRead = 0;
		while (! (token.GetTokenReference() == ')'))
			{
			distCmdOptcmdRead << token.GetTokenReference();
			double d;
			//POL 8-July-2005 added ! before IsADouble in if statement below 
			if (!IsADouble(token.GetTokenReference(), &d)  || d < 0.0)
			 	{
				errState = query_for_error;
				errSnippet.clear();
				errSnippet << "Expecting Dirichlet(<number>, <number> ...) but encountered " << distCmdOptcmdRead << "...";
				return false;
				}
			currentValue->contVariables.push_back(d);
			++token;
			if (token.GetTokenReference() == ',')
				{
				distCmdOptcmdRead << " , ";	
				++token;
				}
			++nRead;	
			}
		if (nRead < 2)
			{
			errState = query_for_error;
			errSnippet.clear();
			errSnippet << "Expecting at least two parameters for the Dirichlet distribution, but encountered " << nRead;
			return false;
			}
		}
	return true;
	}
	
string DistributionCmdOption::GetExpectedDistributionString() const
	{
	string retStr;
	retStr << 'a' << (nVariates == 1 ?" single-variate" : " multi-variate") << " distribution";
	
	if (classOfDistribution == kZeroOrOne)
		retStr << " for a 0/1 quantity";
	else if (classOfDistribution & kZeroToOneBit)
		{
		if (classOfDistribution & kSumToOneConstraintBit)
			retStr << " for positive numbers that sum to 1";
		else
			retStr << " for numbers between 0 and 1";
		}
	else if (classOfDistribution & kNonNegativeBit)
		{
		retStr << " for positive ";
		if (classOfDistribution == kDiscreteNonNegative)
			retStr << "integers";
		else
			retStr << "numbers";
		}
	else if (!(classOfDistribution & kContinuous))
		retStr << " for integers";
	
	return retStr;
	}  

VecString DistributionCmdOption::GetValidArgument()
	{
	return VecString(1, "Some Distribution"); //@ need valid arg format for DistributionCmdOption
	}
bool DistributionCmdOption::CheckClassOfDistribution(DistributionDescription::DistrFamily f)
  	{
	switch (f)
		{
		case DistributionDescription::kBernoulli:
		case DistributionDescription::kBinomial:
			if ( !(classOfDistribution & kDiscrete) || (nVariates > 1))
				{
				errState |= NxsCmdOption::query_for_error;
				errSnippet = "Expecting a ";
				errSnippet << GetExpectedDistributionString() << " (the " << DistributionDescription::GetNameOfDistribution(f) << " is a discrete " << DistributionDescription::GetFundamentalRange(f) << " distribution ).";
				return false;
				}
			return true;
		case DistributionDescription::kUniform:
		case DistributionDescription::kGamma:
		case DistributionDescription::kExponential:
		case DistributionDescription::kInverseGamma:
			if ( !(classOfDistribution & kContinuous) || (nVariates > 1) || ((classOfDistribution & kZeroToOneBit) && f != DistributionDescription::kUniform))
				{
				errState |= NxsCmdOption::query_for_error;
				errSnippet = "Expecting a ";
				errSnippet << GetExpectedDistributionString() << " (the " << DistributionDescription::GetNameOfDistribution(f) << " is a continuous " << DistributionDescription::GetFundamentalRange(f) << " distribution ).";
				return false;
				}
			return true;
		case DistributionDescription::kBeta:
			if ( !(classOfDistribution & kContinuous) || (nVariates > 2))
				{
				errState |= NxsCmdOption::query_for_error;
				errSnippet = "Expecting a ";
				errSnippet << GetExpectedDistributionString() << " (the " << DistributionDescription::GetNameOfDistribution(f) << " is a continuous " << DistributionDescription::GetFundamentalRange(f) << " distribution ).";
				return false;
				}
			return true;
		case DistributionDescription::kDirichlet:
			if ( !(classOfDistribution & kContinuous))
				{
				errState |= NxsCmdOption::query_for_error;
				errSnippet = "Expecting a ";
				errSnippet << GetExpectedDistributionString() << " (the " << DistributionDescription::GetNameOfDistribution(f) << " is a continuous, multivariate " << DistributionDescription::GetFundamentalRange(f) << " distribution ).";
				return false;
				}
			return true;
		case DistributionDescription::kPartOfJoint:
			break;
		}
	errState |= unrecognized;
	return false;			
  	}

DistributionCmdOption::DistributionCmdOption(
  const string &n,
  DistributionDescription *ddToModify, 
  const string &def, 
  ClassOfDistribution c, 
  unsigned numVariates,
  DistributionProvider dp, 
  bool persis, 
  CmdPermissionLevel perm)
  	:SimpleCmdOptionInterface<DistributionDescription>( n, ddToModify, *ddToModify, true,  persis, perm),
  	knownDistNameProvider(dp),
  	classOfDistribution(c),
  	nVariates(numVariates)
  	{
  	PHYCAS_ASSERT(ddToModify);
  	if (!def.empty() && ReadStringAsValue(def))
		defaultValue = *currentValue;
	}
  
  
/*----------------------------------------------------------------------------------------------------------------------
|
|	When adding a new distribution, remember to protect its name in the distributionmanager.cpp's reservedNames so that
|	the user can't apply a generic name to a specific distribution.
*/
bool DistributionCmdOption::ReadValue(
  NxsToken &token, 
  bool equalsAlreadyRead) /* true if an equals is NOT expected */
	{
	if (!equalsAlreadyRead && !EatEqualsThenAdvance(token))
		return false;

	currentValue->Clear();
	if (knownDistNameProvider)
		{
		const DistributionDescription *knownDist  = knownDistNameProvider(token.GetTokenReference().c_str());
		if (knownDist != NULL)
			{
			*currentValue = *knownDist;
			CheckClassOfDistribution(currentValue->family);
			if (WasValidRead())
				{
				++token;
				return true;
				}
			return false;
			}
		}
	distCmdOptcmdRead = token.GetTokenReference();
	if (token.IsAbbreviation("Uniform"))
		{
		if (!ReadUniform(token))
			return false;
		}
	else if (token.IsAbbreviation("BERnoulli"))
		{
		if (!ReadBernoulli(token))
			return false;
		}
	else if (token.IsAbbreviation("BInomial"))
		{
		if (!ReadBinomial(token))
			return false;
		}
	else if (token.IsAbbreviation("Exponential"))
		{
		if (!ReadExponential(token))
			return false;
		}
	else if (token.IsAbbreviation("Gamma"))
		{
		if (!ReadGamma(token))
			return false;
		}
	else if (token.IsAbbreviation("InvGamma") || token.IsAbbreviation("InverseGamma"))
		{
		if (!ReadInvGamma(token))
			return false;
		}
	else if (token.IsAbbreviation("Dirichlet"))
		{
		if (!ReadDirichlet(token))
			return false;
		}
	else if (token.IsAbbreviation("BETa"))
		{
		if (!ReadBeta(token))
			return false;
		}
	else
		{
		errState = unrecognized;
		return false;
		}
	if (EatWordThenAdvance(token, ")", ncl::kStringRespectCase))
		return WasValidRead();
	return false;
	}
  	
const DistributionCmdOption::VecDistParamDesc *DistributionCmdOption::ReadParams(
  const vector<const VecDistParamDesc *> &inVec, 
  NxsToken &token, 
  unsigned index)
	{
	vector<const VecDistParamDesc *> stillViable;
	vector<const VecDistParamDesc *>::const_iterator vecIt = inVec.begin();
	DistParamType pt = kUndefinedParam;
	PHYCAS_ASSERT(token.GetTokenReference().length() > 0);
	//	Narrow down which of the ways of specifying the distribution is being used.
	// we should be able to determine the type of the next token here.  If none are consistent with the token stream, return NULL
	//
	if (isalpha(token.GetTokenReference()[0]))
		{
		distCmdOptcmdRead << token.GetTokenReference();	
		for (; vecIt != inVec.end(); ++vecIt)
			{
			if (index < (*vecIt)->size())
				{
				const DistParamDesc &dpd = (**vecIt)[index];
				if (token.IsAbbreviation(dpd.name))
					{
					PHYCAS_ASSERT(pt == kUndefinedParam || pt == dpd.type);	// we shouldn't have multiple matches with different types of params
					pt = dpd.type;
					stillViable.push_back(*vecIt);
					}
				}
			}
		if (pt == kUndefinedParam)
			return NULL;
		if (!AdvanceThenEatEqualsThenAdvance(token))
			{
			distCmdOptcmdRead << token.GetTokenReference();	
			return NULL;
			}
		distCmdOptcmdRead << " = ";	
		}
	else
		{
		for (; vecIt != inVec.end(); ++vecIt)
			{
			if (index < (*vecIt)->size() )
				{
				const DistParamDesc &dpd = (**vecIt)[index];
				if (dpd.nameIsOptional)
					{
					PHYCAS_ASSERT(pt == kUndefinedParam || pt == dpd.type);	// we shouldn't have multiple matches with different types of params
					pt = dpd.type;
					stillViable.push_back(*vecIt);
					}
				}
			}
		if (pt == kUndefinedParam)
			return NULL;
		}
	//	Read the token into a local temp variable
	//
	distCmdOptcmdRead << token.GetTokenReference();	
	long l = 0;
	double d = 0.0;
	switch(pt)
		{
		case (kIntegerParam) :
		case (kPosIntegerParam):
			if (!IsALong(token.GetTokenReference(), &l))
				 return NULL;
			if (pt == kPosIntegerParam && l < 1)
				return NULL;
			break;
		case (kUnboundedParam):
		case (kPositiveParam): 
		case (kNonNegativeParam): 
		case (kProbParam):
			if (!IsADouble(token.GetTokenReference(), &d))
				 return NULL;
			if (pt == kPositiveParam && d <= 0.0)
				 return NULL;
			if (pt == kNonNegativeParam && d < 0.0)
				 return NULL;
			if (pt == kProbParam && (d < 0.0 || d > 1.0))
				 return NULL;
			break;
		case (kUndefinedParam):
			PHYCAS_ASSERT(false);
			return NULL;
		}
	++token;
	PHYCAS_ASSERT(!stillViable.empty());
	const DistributionCmdOption::VecDistParamDesc *returnVec = NULL;
	//  recurse if there are more parameters to read)
	if (stillViable[0]->size() > index + 1)
		{
		if (token.GetTokenReference() == ',')
			{
			distCmdOptcmdRead << " , ";	
			++token;
			}
		returnVec = ReadParams(stillViable, token, index + 1);
		}
	else
		{
		PHYCAS_ASSERT (stillViable.size() == 1); //it shouldn't be possible reach the end of the param list and still not be certain which description is being used.
		returnVec = stillViable[0];
		currentValue->meaningOfVariable = ((*returnVec)[index].isParam ? DistributionDescription::kParam: DistributionDescription::kMoment);
		}
	if (returnVec != NULL)
		{
		const DistParamDesc &dpd = (*returnVec)[index];
		PHYCAS_ASSERT( (currentValue->meaningOfVariable == DistributionDescription::kMoment) && (!(dpd.isParam)) || (currentValue->meaningOfVariable == DistributionDescription::kParam) && (dpd.isParam));	// we shouldn't allow mixing of descriptions via params and via moments
		//	square all standard deviations to get the variance
		//
		if (currentValue->meaningOfVariable == DistributionDescription::kMoment && EqualsCaseInsensitive(dpd.name, "STDDEV"))
			{
			d *= d;
			}
		if (pt == kIntegerParam || pt == kPosIntegerParam)
			{
			while (currentValue->discreteVariables.size() <= dpd.paramIndex)
				currentValue->discreteVariables.push_back(0);
			PHYCAS_ASSERT(currentValue->discreteVariables[dpd.paramIndex] == 0);
			currentValue->discreteVariables[dpd.paramIndex] = l;
			}
		else
			{
			while (currentValue->contVariables.size() <= dpd.paramIndex)
				currentValue->contVariables.push_back(0.0);
			PHYCAS_ASSERT(currentValue->contVariables[dpd.paramIndex] == 0.0);
			currentValue->contVariables[dpd.paramIndex] = d;
			}
		}
	return returnVec;
	}


bool DistributionCmdOption::ReadBeta(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	if (nVariates > 2)
		{
		errState = query_for_error;
		errSnippet = "Expecting a ";
		errSnippet << GetExpectedDistributionString() << " (the Beta distribution is a single variate distribution) ";
		}
	
	currentValue->family = DistributionDescription::kBeta;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts BETa({Alpha = }0.333333 , {Beta = }1.2), BETa(Beta = 1.2, Alpha = 0.3), BETa(MEAN = 1.0, STDDEV = 1.0),  BETa(MEAN = 1.0, VAR = 1.0)
	//
		
	const DistParamDesc alphaParam("Alpha", kPositiveParam, 0, true, true);
	const DistParamDesc betaParam("Beta", kPositiveParam, 1, true, true);
	const DistParamDesc alphaBackwardParam("Alpha", kPositiveParam, 0, true, false);
	const DistParamDesc betaBackwardParam("Beta", kPositiveParam, 1, true, false);
	const DistParamDesc meanParam = DistParamDesc("MEAN", kProbParam, 0, false, false);
	const DistParamDesc stdDevParam = DistParamDesc("STDDEV", kProbParam, 1, false, false);
	const DistParamDesc varParam = DistParamDesc("VARiance", kProbParam, 1, false, false);
	VecDistParamDesc paramForm,  paramBackwardForm, meanSDForm, meanVarForm;
	paramForm.push_back(alphaParam);
	paramForm.push_back(betaParam);
	paramBackwardForm.push_back(alphaBackwardParam);
	paramBackwardForm.push_back(alphaBackwardParam);
	meanSDForm.push_back(meanParam);
	meanSDForm.push_back(stdDevParam);
	meanVarForm.push_back(meanParam);
	meanVarForm.push_back(varParam);
	vector<const VecDistParamDesc *> parseOptions;
	parseOptions.push_back(&paramForm);
	parseOptions.push_back(&paramBackwardForm);
	parseOptions.push_back(&meanSDForm);
	parseOptions.push_back(&meanVarForm);
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting  Beta (Alpha = <number>, Beta = <number>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	if (currentValue->meaningOfVariable == DistributionDescription::kMoment)
		{
		PHYCAS_ASSERT(currentValue->contVariables.size() == 2);
		if (currentValue->contVariables[1] > (currentValue->contVariables[0] - currentValue->contVariables[0]*currentValue->contVariables[0]))
			{
			errState = query_for_error;
			errSnippet = "the variance of a Beta distribution cannot be greater than (mean - mean*mean)";
			return false;
			}
		}
	
	return true;
	}

bool DistributionCmdOption::ReadUniform(
  NxsToken &token)	/* pointer to token that should hold the next value of the setting in string form */
	{
	currentValue->family = DistributionDescription::kUniform;
	if (!CheckClassOfDistribution(currentValue->family)  || !AdvanceThenEatWordThenAdvance(token,"(", ncl::kStringRespectCase))
		return false;
	distCmdOptcmdRead << '(';	
	//	accepts U({lower =} 0, {upper =} 1)
	//	(commas are optional)
	vector<const VecDistParamDesc *> parseOptions;
	const DistParamDesc lowerB = DistParamDesc("LowerBound", kUnboundedParam, 0, true, true);
	const DistParamDesc upperB = DistParamDesc("UpperBound", kUnboundedParam, 1, true, true);
	/*
		DistParamDesc unNamedBound = DistParamDesc(string(), kUnboundedParam, true);
		VecDistParamDesc unNamed, named, firstNamed, secondNamed;
		unNamed.push_back(unNamedBound); unNamed.push_back(unNamedBound);
		firstNamed.push_back(lowerB); firstNamed.push_back(unNamedBound);
		secondNamed.push_back(unNamedBound); secondNamed.push_back(upperB);
		vector<const VecDistParamDesc *> parseOptions;
		parseOptions.push_back(&unNamed);
		parseOptions.push_back(&named);
		parseOptions.push_back(&firstNamed);
		parseOptions.push_back(&secondNamed);
	*/
	VecDistParamDesc named;
	named.push_back(lowerB);
	named.push_back(upperB);
	parseOptions.push_back(&named);
	
	if (ReadParams(parseOptions, token) == NULL)
		{
		errState = query_for_error;
		errSnippet.clear();
		errSnippet << "Expecting Uniform(LowerBound = <number>, UpperBound = <number>) but encountered " << distCmdOptcmdRead << "...";
		return false;
		}
	PHYCAS_ASSERT(currentValue->contVariables.size() == 2);
	PHYCAS_ASSERT(currentValue->discreteVariables.empty());
	if (currentValue->contVariables[0] >= currentValue->contVariables[1])
		{
		errState = illegal_range;
		errSnippet.clear();
		errSnippet << currentValue->contVariables[0];
		errSnippet << '-';
		errSnippet << currentValue->contVariables[1];
		return false;
		}
	if (classOfDistribution & kZeroToOneBit) 
		{
		if (currentValue->contVariables[0] < 0.0 || currentValue->contVariables[1] > 1.0)
			{
			errState |= NxsCmdOption::query_for_error;
			errSnippet = "Expecting a ";
			errSnippet << GetExpectedDistributionString() << " (the ";
			if (currentValue->contVariables[0] < 0.0)
				errSnippet << "LowerBound of " << currentValue->contVariables[0];
			else
				errSnippet << "UpperBound of " << currentValue->contVariables[1];
			errSnippet << " is illegal)";
			return false;
			}
		}
	return true;
	}

