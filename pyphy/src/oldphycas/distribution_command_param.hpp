#if !defined NCL_NXS_DISTRIBUTION_CMD_OPTION_H
#define NCL_NXS_DISTRIBUTION_CMD_OPTION_H

#include "pyphy/src/ncl/command/nxs_cmd_param.hpp"
#include "pyphy/src/oldphycas/distribution_description.hpp"
#include<boost/function.hpp>

/**
 *	In the definitions below:
 *		- capitalized letters indicate required portions of the string
 *		- {} indicates a part of the definition that can be omitted
 *		- <range-type> indicates a variable of the type described by range-type, so <probaility> means a [0.0, 1.0] (a number between 0 and 1, including the endpoints)
 *	Category	range			nVar 								Definition format
 *	Discrete	0 or 1			1									BERNoulli({Probability =} <probability>) or BERNoulli(Mean = <probability>)
 *	Discrete	0 to Ntrials	1									BINOMial({Probability =} <probability>, {Ntrials =} <pos. int.>) or BINOMial(Ntrials = <pos. int.>, Probability = <probability>)
 *	Continuous	Non-neg			1									EXPonential({Hazard =} <pos. real>) or EXPonential(Mean = <pos. real>) 
 *	Continuous	Non-neg			1									GAMma({SHape = } <pos. real>, {SCale = } <pos. real>), GAMma(SCale = <pos. real>, SHape = <pos. real>), GAMma(MEAN = <pos. real>, STDDEV = <pos. real>),  or GAMma(MEAN = <pos. real>, VAR = <pos. real>)
 *	Continuous	Non-neg			1									INVGAMma({SHape = } <real on (2.0, DBL_MAX)>, {SCale = } <pos. real>), INVGAMma(SCale = <real on (2.0, DBL_MAX)>, SHape = <pos. real>), INVGAMma(MEAN = <pos. real>, STDDEV = <pos. real>), or INVGAMma(MEAN = <pos. real>, VAR = <pos. real>)
 *	Continuous  Sum-to-1	len(paramlist) - 1 or len(paramlist)	DIRichlet(<pos. real>, <pos. real>, ... <pos. real>)
 *	Continuous  Sum-to-1 		1 or 2								BETa({Alpha = } <pos. real> , {Beta = } <pos. real>), BETa(Beta = <pos. real>, Alpha = <pos. real>), or BETa(MEAN = <probability>, STDDEV = <pos. real>)
 *	Continuous	Lower-Upper		1									Uniform({Lower =} <real>, {Upper =} <real on (Lower, DBL_MAX)>)
 *
 */
class DistributionCmdOption : public SimpleCmdOptionInterface<DistributionDescription>
	{
	public:
		typedef boost::function1<const DistributionDescription *,const std::string &> DistributionProvider;
		enum ClassOfDistribution 
				{
				kDiscrete				= 0x01, /*unbounded discrete */
				kContinuous				= 0x02, /*unbounded continuous */
				kDiscreteOrContinuous 	= 0x03, /*any distribution */
				kNonNegativeBit			= 0x04,
				kDiscreteNonNegative	= kNonNegativeBit + kDiscrete,
				kContinuousNonNegative	= kNonNegativeBit + kContinuous,
				kNonNegative			= kNonNegativeBit + kDiscreteOrContinuous, /*discrete or continuous positive distribution */
				kZeroToOneBit			= 0x08, 
				kZeroOrOne				= kZeroToOneBit + kDiscrete, /*discrete and 0-1 */
				kContinuousZeroToOne	= kZeroToOneBit + kContinuous, 
				k_ZeroToOne  			= kZeroToOneBit + kDiscreteOrContinuous, /*discrete or continuous 0-1 distribution */
				kSumToOneConstraintBit	= 0x10,	/*DON'T uses this bit alone!!! use kSumToOne (which has the kContinuous,kZeroToOne and  kSumToOneConstraintBits all on) */
				kSumToOne				= 0x1A
				};
		DistributionCmdOption(const std::string &n, DistributionDescription *ddToModify, const std::string &def, ClassOfDistribution c, unsigned numVariates, DistributionProvider dp, bool persis, CmdPermissionLevel perm);
		std::string		GetCurrentValueAsString() const;
		std::string		GetDisplayType(bool includeIndefiniteArticle = false,  bool plural = false) const;
		std::string		GetExpectedDistributionString() const;
		VecString		GetValidArgument();
		bool			IsCurrentlyValid();
		bool			ReadValue(NxsToken &token, bool equalsAlreadyRead);
		void			RevertValueBecauseCommandFailed();
		void 			StorePreviousValue();
		virtual void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
	private:
		bool			CheckClassOfDistribution(DistributionDescription::DistrFamily f);
		bool			ReadBernoulli(NxsToken &token);
		bool			ReadBeta(NxsToken &token);
		bool			ReadBinomial(NxsToken &token);
		bool			ReadDirichlet(NxsToken &token, bool asBeta = false);
		bool			ReadExponential(NxsToken &token);
		bool			ReadGamma(NxsToken &token);
		bool			ReadInvGamma(NxsToken &token);
		bool			ReadUniform(NxsToken &token);
		
		enum DistParamType
			{
			kUndefinedParam, 
			kIntegerParam,
			kPosIntegerParam, 
			kUnboundedParam, 
			kPositiveParam, 
			kNonNegativeParam, 
			kProbParam
			};
		struct DistParamDesc
			{	
				std::string		name;
				DistParamType 	type;
				bool			isParam;
				bool 			nameIsOptional;
				unsigned 		paramIndex;
				
			
				DistParamDesc(const std::string &n, DistParamType t, unsigned i, bool p, bool o = false)
					:name(n),
					type(t),
					isParam(p),
					nameIsOptional(o),
					paramIndex(i)
					{}
			};
		typedef std::vector<DistParamDesc> VecDistParamDesc;
		
		const VecDistParamDesc *ReadParams(const std::vector<const VecDistParamDesc*> &, NxsToken &token, unsigned index = 0);
		DistributionProvider 	knownDistNameProvider;	//can be NULL
		ClassOfDistribution 	classOfDistribution;
		unsigned				nVariates;			/* 0 means accept any single or multi-variate, otherwise indicates the the number of variates to be described */
		unsigned				nVariatesBefCmd;	/* 0 means accept any single or multi-variate, otherwise indicates the the number of variates to be described */

		std::string distCmdOptcmdRead;//used in composing error messages
				
	};
	

inline bool DistributionCmdOption::IsCurrentlyValid()
	{
	return currentValue->IsValid();
	}
	
inline std::string 	DistributionCmdOption::GetDisplayType(bool includeIndefiniteArticle,  bool plural) const
	{
	return StandardNoun("distribution", includeIndefiniteArticle, plural);
	}

	
inline void DistributionCmdOption::StorePreviousValue()
	{
	nVariatesBefCmd = nVariates;
	SimpleCmdOptionInterface<DistributionDescription>::StorePreviousValue();
	}
	
inline void DistributionCmdOption::RevertValueBecauseCommandFailed()
	{
	nVariates = nVariatesBefCmd;
	SimpleCmdOptionInterface<DistributionDescription>::RevertValueBecauseCommandFailed();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string DistributionCmdOption::GetCurrentValueAsString() const
	{
	return currentValue->GetValueInStringForm();
	}

	
#endif
