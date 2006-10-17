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

#ifndef NCL_NXS_ARRAY_CMD_OPTIONS_H
#define NCL_NXS_ARRAY_CMD_OPTIONS_H
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_command_decl.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
template <class T> class NxsArrayCmdOption: public NxsCmdOption
	{
	
	public :
		typedef boost::shared_ptr< SimpleCmdOptionInterface<T> > ElementReaderPtr;
		NxsArrayCmdOption(const std::string &n, std::vector<T> * manipVar, ElementReaderPtr elReader,  int elementsToRead, bool requireParens, bool persist, CmdPermissionLevel pLevel ); 
		std::string	GetCurrentValueAsString() const;
		std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
		bool		IsCurrentlyValid();
		void		StorePreviousValue();
		bool		ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		void		ReturnValueToDefault();
		void		RevertValueBecauseCommandFailed();
		
		
		
	protected:	
		bool 		ReadNextValue(NxsToken &token);
		
		std::vector<T>			valBefCmd;
		std::vector<T>			defVal;
		std::vector<T>		   *currentValue;
		int					numExpectedElements;	//0 means read an undetermined number
		bool				expectingParens;
		ElementReaderPtr 	elementReader;
	};

template <class T> NxsArrayCmdOption<T>::NxsArrayCmdOption(
  const std::string &n, 
  std::vector<T> * manipVar, 
  ElementReaderPtr elReader,  
  int elementsToRead, 
  bool requireParens, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:NxsCmdOption(n, NxsCmdOption::kDefCouldBeValid, elReader->CheckIfOld(), persist, pLevel),
  	elementReader(elReader), 
  	currentValue(manipVar),
  	numExpectedElements(elementsToRead),
  	expectingParens(requireParens)
  	{
  	valBefCmd = defVal = *currentValue;
  	}

	
template <class T> std::string NxsArrayCmdOption<T>::GetCurrentValueAsString() const
	{
	if (!currentValue->empty())
		{
		std::string retStr = "(";
		for (unsigned i = 0; i < currentValue->size(); ++i)
			{
			elementReader->SetValue((*currentValue)[i]);
			elementReader->GetCurrentValueAsString();
			}
		retStr << ')';
		return retStr;
		}
	else
		return "()";
	}
	
template <class T> std::string NxsArrayCmdOption<T>::GetDisplayType(
  bool includeIndefiniteArticle,
  bool plural) const
	{
	std::string retStr;
	if (includeIndefiniteArticle && !plural)
		retStr << "a ";
	retStr << "list";
	if (plural)
		retStr << 's';
	retStr << " of " << elementReader->GetDisplayType(false, true);
	return retStr;
	}
	
template <class T> bool NxsArrayCmdOption<T>::IsCurrentlyValid()
	{
	T singleElementStorage = elementReader->GetValue();
	for (unsigned i = 0; i < currentValue->size(); ++i)
		{
		elementReader->SetValue((*currentValue)[i]);
		if (!elementReader->IsCurrentlyValid())
			{
			elementReader->SetValue(singleElementStorage);
			return false;
			}
		}
	elementReader->SetValue(singleElementStorage);
	return true;
	}
	
template <class T> void NxsArrayCmdOption<T>::StorePreviousValue()
	{
	valBefCmd = *currentValue;
	elementReader->PrepareToRead();
	}
	
template <class T> void NxsArrayCmdOption<T>::ReturnValueToDefault()
	{
	*currentValue = defVal;
	elementReader->ReturnToDefault();
	}
	
template <class T> void NxsArrayCmdOption<T>::RevertValueBecauseCommandFailed()
	{
	*currentValue = valBefCmd;
	elementReader->RevertBecauseCommandFailed();
	}
	
template <class T> bool NxsArrayCmdOption<T>::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	currentValue->clear();
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		if (expectingParens && EatWordThenAdvance(token, "(", ncl::kStringRespectCase))
			{
			if (numExpectedElements > 0)
				{
				for (unsigned nel = 0; nel < numExpectedElements; ++nel)
					{
					if (!ReadNextValue(token))
						return false;
					}
				}
			else
				{
				if (expectingParens)
					{
					while (token.GetTokenReference() != ')')
						{
						if (!ReadNextValue(token))
							return false;
						}
					}
				else
					{
					NXS_ASSERT(!HadError());
					for (;;)
						{
						if (token.GetTokenReference() == ';')
							break;
						NxsTokenizerState b = token.GetTokenizerState();
						if (!ReadNextValue(token))
							{
							if (GetErrorState() == unrecognized || GetErrorState() == unrecog_labile)
								{
								errState = no_err;
								token.SeekTokenizerState(b);
								}
							else
								return false;
							break;
							}
						}
					}
				}
			if (!HadError() && expectingParens)
				 {
				 if (!EatWordThenAdvance(token, ")", ncl::kStringRespectCase))
				 	return false;
				 }
			}
		return WasValidRead();
		}
	return false;
	}
	
template <class T> bool NxsArrayCmdOption<T>::ReadNextValue(
  NxsToken &token)
	{
	if (!elementReader->ReadValue(token, true))
		{
		FlagError(elementReader->GetErrorState(), elementReader->GetErrorSnippet());
		return false;
		}
	currentValue->push_back(elementReader->GetValue());
	return true;
	}
#endif