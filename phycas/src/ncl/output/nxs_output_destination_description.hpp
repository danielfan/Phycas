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

#if ! defined(NCL_NXS_OUTPUT_DESTINATION_DESCRIPTION_HPP)
#define NCL_NXS_OUTPUT_DESTINATION_DESCRIPTION_HPP
#include "phycas/src/ncl/misc/nxs_file_path.hpp"

class NxsOutputDestinationDescription
	{
	public:
		enum OutputDestinationType
			{
			kSuppressOutput,
			kDirectOutputToNCLStream,
			kDirectOutputToFile,
			kDirectOutputToBoth
			};
		enum ValidNCLOutputRedirect
			{
			kRedirectToOutputStream		= 0,
			kRedirectToCommentStream	= 1,
			kRedirectToErrorStream		= 2,
			kRedirectToPlotStream		= 3
			};
		STATIC_CONST_ACCESSOR std::string TranslateRedirectionStreamName(ValidNCLOutputRedirect);
		STATIC_CONST_ACCESSOR const VecString 	& getLegalRedirectionStreamsNames();
		
		NxsOutputDestinationDescription();
		NxsOutputDestinationDescription(ValidNCLOutputRedirect o);
		NxsOutputDestinationDescription(const std::string &fn, bool appendIfPresent = false, bool replaceIfPresent = false);
		NxsOutputDestinationDescription(ValidNCLOutputRedirect o, const std::string &fn, bool appendIfPresent = false, bool replaceIfPresent = false);
		NxsOutputDestinationDescription &operator=(const NxsOutputDestinationDescription &other);
		bool IsSuppressOutput() const
			{
			return outDestinationType == kSuppressOutput;
			}
		bool IsRedirectToNCLStream() const
			{
			return outDestinationType == kDirectOutputToNCLStream ||  outDestinationType == kDirectOutputToBoth;
			}
		bool IsRedirectToFile() const
			{
			return outDestinationType == kDirectOutputToFile ||  outDestinationType == kDirectOutputToBoth;
			}
		std::string GetRedirectionStreamName(unsigned i) const
			{
			return TranslateRedirectionStreamName(redirectStreams[i]);
			}
		std::string GetNexusDescription() const;
		
		typedef std::vector<ValidNCLOutputRedirect> VecRedirection;
		
		OutputDestinationType   outDestinationType;
		VecRedirection			redirectStreams;
		mutable NxsOutFilePath  filePath; 
	};

#endif
