#ifndef NCL_NXS_ADVANCED_CMD_OPTION_H
#define NCL_NXS_ADVANCED_CMD_OPTION_H
#include <boost/bind.hpp>
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_primitive_cmd_param.hpp"
#include "ncl/misc/nxs_index_set.hpp"
	
/*----------------------------------------------------------------------------------------------------------------------
|	A CmdOption that will read in NEXUS set descriptions.
|	To do this it must be supplied with functors that can tell it the maximum number (max index + 1) to accept, how to translate a 
|	label to an index and how to translate a label to a set.
|	depending on the bool shouldExpectName in the constructor, the NxsIndexSetCmdOption can be told to read  
|	"name = set-description" as opposed just "set-descripion"
*/
class NxsIndexSetCmdOption: public SimpleCmdOptionInterface<NxsIndexSet>
	{
	public:
		class	NxsX_MaxExceeded
			{
			public :
				unsigned maxV;
				unsigned suppliedVal;
			
				NxsX_MaxExceeded(unsigned m, unsigned s) : maxV(m), suppliedVal(s){}
			};
		typedef boost::function0< unsigned >	MaxNumberSource;
		typedef boost::function1<unsigned, const std::string &>	LabelDecoder;
		typedef boost::function1<const NxsIndexSet *, const std::string &>	SetDecoder;
		
		
		void		AllowIndices(bool v) {allowNumbers = v;}
		std::string GetCurrentValueAsString() const;
		std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
		unsigned	GetMaxIndex() const {return maxIndex;}
		std::string GetPrefix() const {return pref;}
		VecString   GetValidArgument();
		bool		IsCurrentlyValid();
		bool		ReadValue(NxsToken &token, bool equalsAlreadyRead);
		void		ReturnValueToDefault();
		void 		RefreshMaxIndex();
		void		WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
					NxsIndexSetCmdOption(const std::string &n, NxsIndexSet *manipVal, const std::string &def, MaxNumberSource  maxProv, LabelDecoder lab , SetDecoder setP, bool shouldExpectName, bool  persist, std::string prefixWord, CmdPermissionLevel pLevel );
		
	protected:
	
		
		LabelDecoder	   	labelProvider;
		SetDecoder			setProvider;
		MaxNumberSource 	maxNumberProv;
		unsigned			maxIndex;
		const bool			expectingName;	/* */
		bool				allowNumbers;
		const std::string 	pref;// name prefix, ( e.g. tax, char, or tree)
		const std::string		defStr;
	private:
		unsigned 	GetIndexFromLabel(const std::string & nextTok);
		unsigned 	GetIndexFromToken(const std::string & nextTok);
		bool		ReadVectorFormOfSet(NxsToken &token);
		unsigned	ReadStride(NxsToken &token, unsigned maxStride);
		bool		ReadStandardFormOfSet(NxsToken &token);
	};
	
typedef boost::shared_ptr< NxsIndexSetCmdOption > SetCmdOptionShPtr;

inline std::string NxsIndexSetCmdOption::GetDisplayType(
  bool includeIndefiniteArticle,
  bool plural) const
	{
	std::string p = pref;
	p << "Set";	
	return StandardNoun(p,includeIndefiniteArticle, plural);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses the maxNumberProv to get the largest index accepted
*/
inline void NxsIndexSetCmdOption::RefreshMaxIndex()
	{
	if (maxNumberProv)
		{
		maxIndex = maxNumberProv();
  		--maxIndex;
  		}
  	}
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string NxsIndexSetCmdOption::GetCurrentValueAsString() const
	{
	std::string s;
	s << *currentValue;
	return s;
	}

template <class Mgr>
NxsIndexSetCmdOption *InstantiateSetOptionGeneric(const std::string &n, NxsIndexSet *manip, const std::string &d, Mgr *mgrPtr, bool  persist, const std::string &pref, CmdPermissionLevel pLevel);

template <class Mgr>
inline NxsIndexSetCmdOption *InstantiateSetOptionGeneric(const std::string &n, NxsIndexSet *manip, const std::string &d, Mgr *mgrPtr, bool  persist, const std::string &pref, CmdPermissionLevel pLevel)
	{
	NxsIndexSetCmdOption::MaxNumberSource mn = NxsIndexSetCmdOption::MaxNumberSource(boost::bind(&Mgr::GetSize, mgrPtr));
	NxsIndexSetCmdOption::LabelDecoder indS = NxsIndexSetCmdOption::LabelDecoder(boost::bind(&Mgr::FindIndex, mgrPtr, _1));
	NxsIndexSetCmdOption::SetDecoder setS = NxsIndexSetCmdOption::SetDecoder(boost::bind(&Mgr::GetSet, mgrPtr, _1));
	const bool expectName = false; // should be an argument
	return new NxsIndexSetCmdOption(n, manip, d, mn, indS, setS, expectName, persist, pref, pLevel);
	}

#endif
