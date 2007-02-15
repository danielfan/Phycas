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
#include "phycas/src/ncl/output/nxs_output.hpp"
#include "phycas/src/ncl/output/nxs_output_stream_wrapper.hpp"
#include "phycas/src/ncl/output/nxs_output_destination_description.hpp"
#include "phycas/src/ncl/misc/string_extensions.hpp"
using std::ofstream;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::size_t; // @VKJ 11/23 needed to compile
#endif

const VecString & NxsOutputDestinationDescription::getLegalRedirectionStreamsNames()
	{
	STATIC_CONST const VecString allowedRedirectionStreamsNamesVec = SplitString("Output|Comment|Error|Plot", '|');
	return allowedRedirectionStreamsNamesVec;
	}

string NxsOutputDestinationDescription::TranslateRedirectionStreamName(ValidNCLOutputRedirect vnor)
	{
	const VecString & allowedRedirectionStreamsNamesVec = NxsOutputDestinationDescription::getLegalRedirectionStreamsNames();
	return allowedRedirectionStreamsNamesVec[(size_t) vnor];
	}

NxsOutputStreamWrapper::~NxsOutputStreamWrapper()
	{
	ncl::flush(*this);
	delete filePtr;
	}
		


NxsOutputStreamWrapperShPtr NxsOutputManager::GetGenericOutputStreamShPtr(const NxsOutputDestinationDescription & nodd)
	{
	return NxsOutputStreamWrapperShPtr(new NxsOutputStreamWrapper(nodd, *this));
	}

NxsOutputStreamWrapper::NxsOutputStreamWrapper(const NxsOutputDestinationDescription & nodd, NxsOutputManager &outMgr)
	:plotPtr(NULL),
	filePtr(NULL)
	{
	if (nodd.IsSuppressOutput())
		return;
	if (nodd.IsRedirectToFile())
		ResetFile(nodd.filePath);
	if (nodd.IsRedirectToNCLStream())
		{
		std::set<NxsWritableStream *> streamSet;
		typedef NxsOutputDestinationDescription::VecRedirection::const_iterator RedirectIterator;
		for (RedirectIterator rIt = nodd.redirectStreams.begin(); rIt != nodd.redirectStreams.end(); ++rIt)
			{
			switch (*rIt)
				{
				case (NxsOutputDestinationDescription::kRedirectToOutputStream):
					streamSet.insert(outMgr.GetOutputStreamPtr());
					break;
				case (NxsOutputDestinationDescription::kRedirectToCommentStream):
					streamSet.insert(outMgr.GetOutputCommentStreamPtr());
					break;
				case (NxsOutputDestinationDescription::kRedirectToErrorStream):
					streamSet.insert(outMgr.GetErrorStreamPtr());
					break;
				case (NxsOutputDestinationDescription::kRedirectToPlotStream):
					SetPlotStreamPtr(outMgr.GetPlotStreamPtr());
					break;
				}
			}
		std::copy(streamSet.begin(), streamSet.end(), std::back_inserter(streamPtrs));
		}
	}

void NxsOutputStreamWrapper::ResetFile(NxsOutFilePath & nofp)
	{
	delete filePtr;
	ofstream * fp = new ofstream();
	if (nofp.Open(*fp))
		filePtr = fp;
	else
		delete fp;
	}

