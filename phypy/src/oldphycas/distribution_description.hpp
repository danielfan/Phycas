#if !defined NCL_DISTRIBUTION_DESCRIPRIPTION_ITEM_H
#define NCL_DISTRIBUTION_DESCRIPRIPTION_ITEM_H

#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/misc/string_extensions.hpp"
#include "phypy/src/ncl/output/temp_basic_output_operators.hpp"
#include "phypy/src/probability_distribution.hpp"

typedef boost::shared_ptr<phycas::ProbabilityDistribution> ProbDistPtr; 
typedef boost::shared_ptr<phycas::MultivariateProbabilityDistribution> MVProbDistPtr; 

class DistributionDescription // shouldn't we have a separate types for discrete and continuous
	{
	public:
		enum DistrFamily {
			kPartOfJoint,
			kBernoulli,
			kBeta,
			kBinomial,
			kDirichlet,
			kExponential,
			kGamma,
			kInverseGamma,
			kUniform
			};
		enum VarMeaning {	/* the stored variables either represent the moment of the distribution or the paramters*/
			kParam  		= 0x01,	/* variables store the parameters of the distribution */
			kMoment			= 0x02	/* variables store the moments of the distribution */
			};

		DistributionDescription();
		DistributionDescription(const DistributionDescription &);	
		~DistributionDescription();
		void				Clear(); 
		ProbDistPtr			CreateProbabilityDistribution() const; 
		MVProbDistPtr		CreateMultiVariateProbabilityDistribution() const; 
		std::string 		GetValueInStringForm() const;
		std::string			GetDistributionName() const;
		unsigned			GetNVariates() const;
		unsigned			GetVariateIndex() const;
		STATELESS_FUNC std::string	GetNameOfDistribution(DistrFamily);
		STATELESS_FUNC std::string 	GetFundamentalRange(DistrFamily);
		bool				IsDiscrete() const;
		bool				IsContinuous() const {return !IsDiscrete();}
		bool				IsMultivariate() const;	
		bool				IsPartOfJointMultivariate() const;
		bool				IsValid() const;
		void				SetJointMultivariate(const DistributionDescription &d, unsigned ind);
		DistributionDescription &operator=(const DistributionDescription &);	
	protected:
		
		DistributionDescription &Copy(const DistributionDescription &r);	//used in copy const and = op

		DistrFamily 				family;
		VarMeaning  				meaningOfVariable;
		std::vector<double> 		contVariables;
		std::vector<long>	 		discreteVariables;
		DistributionDescription	  *	jointDist;
		unsigned 					indexInJointDist;
		//@ to do store a vector of DistributionDescription in case the user specifies a distribution for a parameter
	
	private :
	
		friend class DistributionCmdOption;
	};

template<class OUT_STREAM>
class GenericPrinterClass<kConciseOutStyle, DistributionDescription, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const DistributionDescription & dd)
			{
			outStream << dd.GetValueInStringForm();
			}
	};

inline DistributionDescription::~DistributionDescription()
	{
	delete jointDist;
	}
	
inline bool DistributionDescription::IsPartOfJointMultivariate() const
	{
	PHYCAS_ASSERT(family != kPartOfJoint || jointDist != NULL);
	return (family == kPartOfJoint);
	}
	
inline unsigned	DistributionDescription::GetVariateIndex() const
	{
	return indexInJointDist;
	}
	
inline void	DistributionDescription::SetJointMultivariate(
  const DistributionDescription &d, 
  unsigned ind)
	{
	Clear();
	jointDist = new DistributionDescription(d);
	indexInJointDist = ind;
	family = kPartOfJoint;
	}
	
inline bool DistributionDescription::IsMultivariate() const
	{
	return (family == kDirichlet);
	}

inline unsigned	DistributionDescription::GetNVariates() const
	{
	if (IsMultivariate())
		{
		PHYCAS_ASSERT(family == kDirichlet);
		return (unsigned)contVariables.size();
		}
	return 1;
	}	
	
inline DistributionDescription::DistributionDescription(const DistributionDescription &r)
	:jointDist(NULL)	//need to initialize joint Dist to NULL so Copy doesn't delete unallocated memory!
	{
	Copy(r);
	}
	
inline DistributionDescription &DistributionDescription::operator=(const DistributionDescription &r)
	{
	return Copy(r);
	}
	
inline std::string DistributionDescription::GetNameOfDistribution(DistrFamily f)
	{
	std::string s;
	switch (f)
		{
		case DistributionDescription::kBernoulli: return "Bernoulli";
		case DistributionDescription::kBeta: return "Beta";
		case DistributionDescription::kBinomial: return "Binomial";
		case DistributionDescription::kDirichlet: return "Dirichlet";
		case DistributionDescription::kUniform: return "Uniform";
		case DistributionDescription::kGamma: return "Gamma";
		case DistributionDescription::kExponential: return "Exponential";
		case DistributionDescription::kInverseGamma: return "InverseGamma";
		case DistributionDescription::kPartOfJoint: break;
		}
	PHYCAS_ASSERT(0);
	return std::string();
	}

inline std::string DistributionDescription::GetFundamentalRange(DistrFamily f)
	{
	switch (f)
		{
		case DistributionDescription::kBernoulli: 
		case DistributionDescription::kBeta: 
		case DistributionDescription::kDirichlet: 
													return "0-1";
		case DistributionDescription::kUniform: 	return "unbounded";
		case DistributionDescription::kBinomial: 
		case DistributionDescription::kGamma: 
		case DistributionDescription::kExponential:
		case DistributionDescription::kInverseGamma: 
													return "positive";
		case DistributionDescription::kPartOfJoint: break;
		}
	PHYCAS_ASSERT(0);
	return std::string();
	}

inline bool DistributionDescription::IsDiscrete() const
	{
	return (family == kBernoulli || family == kBinomial);
	}

inline std::string DistributionDescription::GetDistributionName() const
	{
	if (family == kPartOfJoint)
		{
		PHYCAS_ASSERT(jointDist != NULL);
		std::string s;
		s << '{' << indexInJointDist << '}' << jointDist->GetDistributionName();
		}
	return GetNameOfDistribution(family);
	}

#endif

