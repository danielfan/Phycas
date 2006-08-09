#ifndef NCL_NXS_RESTRICTED_STRING_CMD_OPTION_H
#define NCL_NXS_RESTRICTED_STRING_CMD_OPTION_H

#include "ncl/nxs_defs.hpp"
#include "ncl/misc/nxs_index_set.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
#include <boost/function.hpp>
/*----------------------------------------------------------------------------------------------------------------------
|	a NxsStringCmdOption that does not allow the value to be a single punctuation character (the punctuation list is 
|	(){}"'-[]/\,;:=*`+<>_	
*/
class NxsWordCmdOption : public NxsStringCmdOption
	{
	public :
		virtual bool	IsCurrentlyValid();
		NxsWordCmdOption( const std::string &n, std::string *manipVal, std::string def, bool validityIsLabile, bool persist, CmdPermissionLevel pLevel );
	};	

/*----------------------------------------------------------------------------------------------------------------------
|	A NxsWordCmdOption that will not accept some reserved words.  
|	The list of illegal words can be supplied as a vector of strings, a |-separated string, or (if the list changes)
|	as a shared_ptr to a ValueHolder<VecString>
*/
class RestrictNameCmdOption : public NxsWordCmdOption
	{
	public:
		typedef boost::function0<VecString> VecStringSource;
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, VecStringSource tabooProvider, bool persist, CmdPermissionLevel pLevel );
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, VecString tabooVec, bool persist, CmdPermissionLevel pLevel );
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, const char *tabooStr, bool persist, CmdPermissionLevel pLevel );
		bool	IsCurrentlyValid();
		void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;		
	protected:
		VecStringSource tabooListProvider;
		mutable VecString 		tabooList;
	private:
		inline void		RefreshTabooList() const;
	};

typedef boost::shared_ptr< RestrictNameCmdOption > RestrictNameCmdOptionShPtr;

inline void RestrictNameCmdOption::RefreshTabooList() const
	{
	if (tabooListProvider)
		tabooList = tabooListProvider();
	}


#endif

