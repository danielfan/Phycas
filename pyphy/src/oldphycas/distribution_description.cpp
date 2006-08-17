//#include "phycas/force_include.h"
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/oldphycas/distribution_description.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
using std::string;

DistributionDescription::DistributionDescription()
	:family(kUniform),
	meaningOfVariable(kParam), 
	jointDist(NULL),
	indexInJointDist(0)
	{
	}
	
void DistributionDescription::Clear()
	{
	family = kUniform;
	meaningOfVariable = kParam; 
	delete jointDist;
	jointDist = NULL;
	indexInJointDist = 0;
	contVariables.clear();
	discreteVariables.clear();
	}
	

string DistributionDescription::GetValueInStringForm() const
	{
	if (!IsValid())
		return string("not valid");
	string retStr;
	retStr << GetDistributionName() << '(';
	unsigned i; 
	switch (family)
		{
		case kPartOfJoint:
			PHYCAS_ASSERT(jointDist != NULL);
			retStr << "The " << GetOrderString(indexInJointDist) << " element in";
			retStr += jointDist->GetValueInStringForm();
			break;
		case kBernoulli:
			retStr << contVariables[0];
			break;
		case kBinomial:
			retStr << contVariables[0] << ',' << discreteVariables[0];
			break;
		case kDirichlet:
			retStr << contVariables[0];
			for (i = 1; i < contVariables.size(); ++i)
				retStr << ',' << contVariables[i];	
			break;
		case kUniform:
			retStr << contVariables[0] << ',' << contVariables[1];
			break;
		case kGamma:
		case kBeta:
		case kInverseGamma:
			if (meaningOfVariable == kParam)
				retStr << contVariables[0] << ',' << contVariables[1];
			else
				retStr << "Mean=" << contVariables[0] << ", Variance =" << contVariables[1];
			break;
		case kExponential:
			if (meaningOfVariable == kParam)
				retStr << contVariables[0];
			else
				retStr << "Mean=" << 1.0/contVariables[0];
			break;
		default:
			return string();
		}
	retStr << ')';
	return retStr;
	}



MVProbDistPtr DistributionDescription::CreateMultiVariateProbabilityDistribution() const
	{
#	if defined NCL_NXS_THROW_UNDEFINED
		if (!IsValid())
			throw NxsX_UndefinedException("Invalid Distribution ", __FILE__, __LINE__); 
#	endif			
	MVProbDistPtr	mvProbDist;
	switch (family)
		{
		case kPartOfJoint:
 			return jointDist->CreateMultiVariateProbabilityDistribution();
		case kDirichlet:
			if (mvProbDist == MVProbDistPtr())
				mvProbDist  = MVProbDistPtr(new phycas::DirichletDistribution(contVariables));
			break;
		default:
			PHYCAS_ASSERT(0);
			return MVProbDistPtr();
		}
	return mvProbDist;
	}
	
ProbDistPtr DistributionDescription::CreateProbabilityDistribution() const
	{
	ProbDistPtr probDist;
	if (IsMultivariate())	
		return ProbDistPtr();
#	if defined NCL_NXS_THROW_UNDEFINED
		if (!IsValid())
			throw NxsX_UndefinedException("Invalid Distribution ", __FILE__, __LINE__); 
#	endif			
	double firstParam, secondParam;
	switch (family)
		{
		case kBernoulli:
			probDist = ProbDistPtr(new phycas::BernoulliDistribution(contVariables[0]));
			break;

		case kBinomial:
			probDist = ProbDistPtr(new phycas::BinomialDistribution(discreteVariables[0], contVariables[0]));	//new_distr
			break;

		case kUniform:
			probDist = ProbDistPtr(new phycas::UniformDistribution(contVariables[0], contVariables[1]));
			break;

		case kBeta:
			if (meaningOfVariable == kParam)
				{
				firstParam = contVariables[0];
				secondParam = contVariables[1];
				probDist = ProbDistPtr(new phycas::BetaDistribution(firstParam, secondParam));
				}
			else
				{
				probDist = ProbDistPtr(new phycas::BetaDistribution(1.0, 1.0));
				probDist->SetMeanAndVariance(contVariables[0], contVariables[1]);
				}
			break;

		case kGamma:
			if (meaningOfVariable == kParam)
				{
				firstParam = contVariables[0];
				secondParam = contVariables[1];
				}
			else
				{
				secondParam = contVariables[1]/contVariables[0];
				firstParam = contVariables[0]/secondParam;
				}
			probDist = ProbDistPtr(new phycas::GammaDistribution(firstParam, secondParam));
			break;

		case kExponential:
			if (meaningOfVariable == kParam)
				firstParam = contVariables[0];
			else
				firstParam = 1.0/contVariables[0];
			probDist = ProbDistPtr(new phycas::ExponentialDistribution(firstParam));
			break;

		case kInverseGamma:
			if (meaningOfVariable == kParam)
				{
				firstParam = contVariables[0];
				secondParam = contVariables[1];
				}
			else
				{
				firstParam = 2.0 + (contVariables[0] * contVariables[0] /contVariables[1]);
				secondParam = 1.0 / ((firstParam- 1.0) * contVariables[0]);
				}
			probDist = ProbDistPtr(new phycas::InverseGammaDistribution(firstParam, secondParam));
			break;

		default:
			PHYCAS_ASSERT(0);
			return ProbDistPtr();
		}
	return probDist;
	}
	
bool DistributionDescription::IsValid() const
	{
	unsigned i; 
	switch (family)
		{
		case kPartOfJoint:
			if (jointDist == NULL || !jointDist->IsValid())
 				return false;
 			return (indexInJointDist < jointDist->GetNVariates());
		case kBernoulli:
			if (contVariables.empty())
				return false;
			return (contVariables[0] >= 0.0 && contVariables[0] <=1.0);
		case kBeta:
			if (contVariables.size() != 2)
				return false;
			if (contVariables[0] < 0.0 || contVariables[1] < 0.0)
				return false;
			if (meaningOfVariable == kMoment)
				{
				if (contVariables[0] > 1.0 || contVariables[1] > 1.0 || (contVariables[0]*contVariables[0] + contVariables[1]) > contVariables[0])
					return false;
				}
			return true;
		case kDirichlet:
			if (contVariables.size() < 2)
				return false;
			for (i = 0; i < contVariables.size(); ++i)
				{
				if (contVariables[i] <0.0)
					return false;
				}
			return true;
		case kBinomial:
			if (contVariables.empty() || discreteVariables.empty())
				return false;
			if (discreteVariables[0] < 1)
				return false;
			return (contVariables[0] >= 0.0 && contVariables[0] <=1.0);

		case kUniform:
			if (contVariables.size() < 2)
				return false;
			return (contVariables[0] < contVariables[1]);

		case kGamma:
			if (contVariables.size() < 2)
				return false;
			return (contVariables[0] > 0.0 && contVariables[1] > 0.0);

		case kExponential:
			if (contVariables.empty())
				return false;
			return (contVariables[0] > 0.0);

		case kInverseGamma:
			if (contVariables.size() < 2)
				return false;
			if (meaningOfVariable == kParam && contVariables[0] <= 2.0 )
				return false;
			return (contVariables[0] > 0.0 && contVariables[1] > 0.0);
			
		}
	return false;
	}


DistributionDescription &DistributionDescription::Copy(const DistributionDescription &r)
	{
	family = r.family;
	meaningOfVariable = r.meaningOfVariable; 
	if (r.jointDist != NULL)
		{
		if (jointDist == NULL)
			jointDist = new DistributionDescription(*r.jointDist);
		else
			*jointDist = *r.jointDist;
		}
	else
		{
		if (jointDist != NULL)
			{
			delete jointDist;
			jointDist = NULL;
			}
		}
	indexInJointDist = r.indexInJointDist;
	contVariables = r.contVariables;
	discreteVariables = r.discreteVariables;
	return *this;
	}

