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

#ifndef NCL_NXS_META_CMD_OPTIONS_H
#define NCL_NXS_META_CMD_OPTIONS_H
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/command/nxs_command_decl.hpp"
#include "phycas/src/ncl/command/nxs_cmd_param.hpp"
#include <boost/mem_fn.hpp>

#if defined TEMPLATED_MIXED_CMD_OPTION
#	error TEMPLATED version is old (does not handle a vector of sub-options)
	/*----------------------------------------------------------------------------------------------------------------------
	|	used to hold user options that can be one of multiple types (e.g. kappa = 3.0 or kappa = estimate)
	|	
	*/
	template <class FirstCmdOpt, class SecondCmdOpt> class NxsMixedCmdOption : public NxsCmdOption
		{
		
		public :
			virtual std::string 	GetCurrentValueAsString() const = 0;
			virtual VecString		GetLegalChoices() const;
			virtual std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
			virtual bool		AbbreviateChoices(); //base class does nothing - only choice commands do something substantive
					bool		IsCurrentlyValid();
			virtual void		StorePreviousValue(readStatus);
			virtual	bool		ReadValue(NxsToken &token, bool equalsAlreadyRead = false);
			void				ReturnValueToDefault();
			virtual void		RevertValueBecauseCommandFailed();
			void				TryToRead(NxsToken &token, bool equalsAlreadyRead);
			
			NxsMixedCmdOption(const std::string &n, boost::shared_ptr<FirstCmdOpt> firstSub,  boost::shared_ptr<SecondCmdOpt> sectSub,  bool persist, CmdPermissionLevel pLevel ); 
			
		protected:	
			

			boost::shared_ptr<FirstCmdOpt> firstSubOption;
			boost::shared_ptr<SecondCmdOpt> secondSubOption;
			bool 				 useFirst;
			bool				 useFirstBefCmd;
			bool				 defaultUseFirst;
		};

	template <class FirstCmdOpt, class SecondCmdOpt> inline NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::NxsMixedCmdOption(
	  const std::string &n, 
	  boost::shared_ptr<FirstCmdOpt>  firstSub, 
	  boost::shared_ptr<SecondCmdOpt>  sectSub,  
	  bool persist, 
	  CmdPermissionLevel pLevel)
	  	:NxsCmdOption(n, persist, pLevel),
	  	firstSubOption(firstSub),
	  	secondSubOption(sectSub),
	  	useFirst(true),
	  	useFirstBefCmd(true),
	  	defaultUseFirst(true)
	  	{
	  	}
	  	
	template <class FirstCmdOpt, class SecondCmdOpt> std::string NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::GetCurrentValueAsString() const
		{
		if (useFirst)
			return firstSubOption->GetCurrentValueAsString();
		else
			return secondSubOption->GetCurrentValueAsString();
		}
			
	template <class FirstCmdOpt, class SecondCmdOpt> std::string NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::GetDisplayType(bool indef,
	  bool plural) const
		{
		std::string retVal = firstSubOption->GetDisplayType(indef, plural);
		retVal << " or " << secondSubOption->GetDisplayType(indef, plural); 
		return retVal;
		}
			
	template <class FirstCmdOpt, class SecondCmdOpt> bool NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::AbbreviateChoices() 
		{
		return (firstSubOption->AbbreviateChoices() && secondSubOption->AbbreviateChoices());
		}	
		
	template <class FirstCmdOpt, class SecondCmdOpt> bool NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt>::IsCurrentlyValid()
		{
		NxsCmdOption *subToUse = ( useFirst ? firstSubOption.get() : secondSubOption.get());
		if (!subToUse->IsCurrentlyValid())
			{
			const std::string &eSnip = subToUse->GetErrorSnippet();
			FlagError(subToUse->GetErrorState(), eSnip);
			return false;
			}
		return true;
		}
		
	template <class FirstCmdOpt, class SecondCmdOpt> void NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::StorePreviousValue(readStatus nextStatus)
		{
		useFirstBefCmd = useFirst;
		firstSubOption->PrepareToRead();
		secondSubOption->PrepareToRead();
		}
		
	template <class FirstCmdOpt, class SecondCmdOpt> bool NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt>::ReadValue(
	 NxsToken &token,	/* the stream of tokens that are being read */
	 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	  	{
	  	useFirst = true;
	  	NxsTokenizerState b = token.GetTokenizerState();
	  	if (!firstSubOption->ReadValue(token, equalsAlreadyRead))
	  		{
	  		token.SeekTokenizerState(b);
	  		if (!secondSubOption->ReadValue(token, equalsAlreadyRead))
	  			return false;
	  		else
	  			useFirst = false;
	  		}
	  	return WasValidRead();
	  	}
	  	
	template <class FirstCmdOpt, class SecondCmdOpt> void NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::RevertValueBecauseCommandFailed()
		{
		useFirst = useFirstBefCmd;
		firstSubOption->RevertBecauseCommandFailed();
		secondSubOption->RevertBecauseCommandFailed();
		}
		
	template <class FirstCmdOpt, class SecondCmdOpt> void NxsMixedCmdOption<FirstCmdOpt, SecondCmdOpt> ::ReturnValueToDefault()
		{
		useFirst = defaultUseFirst;
		if (useFirst)
			firstSubOption->ReturnToDefault();
		else
			secondSubOption->ReturnToDefault();
		}
#else

	/*----------------------------------------------------------------------------------------------------------------------
	|	used to hold user options that can be one of multiple types (e.g. kappa = 3.0 or kappa = estimate)
	|	
	*/
	class NxsMixedCmdOption : public NxsCmdOption
		{
		
		public :
							NxsMixedCmdOption(const std::string &n, unsigned *optIndex, unsigned defIndex, VecNxsCmdOptionShPtr suboptionVec, bool persist, CmdPermissionLevel pLevel ); 
			std::string 	GetCurrentValueAsString() const;
			VecString		GetLegalChoices() const;
			std::string		GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
			VecString		GetValidArgument();
			bool			AbbreviateChoices(); //base class does nothing - only choice commands do something substantive
			bool			IsCurrentlyValid();
			void			StorePreviousValue();
			bool			ReadValue(NxsToken &token, bool equalsAlreadyRead = false);
			void			ReturnValueToDefault();
			void			RevertValueBecauseCommandFailed();
			void			TryToRead(NxsToken &token, bool equalsAlreadyRead);
			virtual void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;

		protected:	
			unsigned * const 	subOptionIndex;
			VecNxsCmdOptionShPtr 	subOptVector;
			unsigned 			indexBefCmd;
			unsigned			defaultIndex;
		};

	

#endif	//TEMPLATED_MIXED_CMD_OPTION

#endif

