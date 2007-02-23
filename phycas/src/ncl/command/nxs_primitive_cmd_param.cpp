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

//#include "phycas/force_include.h"
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/command/nxs_primitive_cmd_param.hpp"
using std::string;

BoolCmdOption::BoolCmdOption(
  const string &n, 
  bool *manipVal,
  bool def, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<bool>(n, manipVal, def, false, persist, pLevel),
  	recognizedAsConverse(false)
  	{
  	converseAbbrev << "NO" << abbreviation;
  	}


bool BoolCmdOption::CanReadKeyword(
  const string &s, 
  int permLevel) const
	{
	if (HasPermission(permLevel))
		{
		if (IsCapAbbreviation(s, GetAbbreviation()))
			{
			recognizedAsConverse = false;
			return true;
			}
		if (IsCapAbbreviation(s, converseAbbrev))
			{
			recognizedAsConverse = true;
			return true;
			}
		}
	return false;
	}
	
string BoolCmdOption::GetConverseName() const
	{
	string r;
	r << "NO" <<  GetName();
	return r;
	}
		
	
	
void BoolCmdOption::SetConverseAbbreviation(
  const string &ca)
	{
	converseAbbrev = ca;
	}
				

bool BoolCmdOption::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	if (recognizedAsConverse)
		{
		*currentValue = false;
		recognizedAsConverse = false;	// set recognizedAsConverse back to false so that ReadValue has the normal behavious by default
		return WasValidRead();
		}
	else
		{
		if (equalsAlreadyRead || token.GetTokenReference() == '=')
			{
			if (!equalsAlreadyRead)
				++token;
			if (IsCapAbbreviation(token.GetTokenReference(), "Yes") || IsCapAbbreviation(token.GetTokenReference(), "True"))
				*currentValue  = true;
			else if (IsCapAbbreviation(token.GetTokenReference(), "No") || IsCapAbbreviation(token.GetTokenReference(), "False"))
				*currentValue  = false;
			else
				return FlagError(NxsCmdOption::unrecognized);
			if (!HadError())
				++token;
			}
		else
			*currentValue = true;
		}
	return WasValidRead();
	}

QuotedStringCmdOption::QuotedStringCmdOption(const string &n, string *manipVal, string def, bool persist, CmdPermissionLevel pLevel)
	:NxsStringCmdOption(n, manipVal, def, false, persist, pLevel)
	{
	}

bool QuotedStringCmdOption::ReadValue( 
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		if (token.GetTokenReference() == '\"')
			{
			token.ReadAllTokensUntil('\"');
			*currentValue = token.GetTokenReference();
			if (WasValidRead())
				{
				++token;
				return true;
				}
			}
		else
			FlagError(unexpected_char, "\"");
		}	
	return false;
	}


