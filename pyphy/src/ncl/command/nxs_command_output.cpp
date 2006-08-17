// This file is included by nxs_characters_block.cpp to avoid a bizarre anonymous namespace multiple definition link error that TL is getting on Mac 10.3.9 (gcc 3.3)
//#if defined (INCLUDE_TO_AVOID_LINK_ERROR)

//#include "phycas/force_include.h"
#include "pyphy/src/ncl/command/nxs_command_output.hpp"
#include "pyphy/src/ncl/command/nxs_choice_cmd_param.hpp"
#include "pyphy/src/ncl/command/nxs_file_cmd_param.hpp"
#include "pyphy/src/ncl/command/nxs_mixed_cmd_param.hpp"
#include "pyphy/src/ncl/command/nxs_restricted_string_cmd_param.hpp"
#include "pyphy/src/oldphycas/distribution_command_param.hpp"
#include "pyphy/src/ncl/command/nxs_set_cmd_param.hpp"

#if defined(NCL_SOCKET_IO)

	void NxsChoiceCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const
		{
		RefreshChoices();
		NxsSaxOutputWrapper tiWrapper(outStream, "choice_type_info");
		for (VecString::const_iterator cIt = choices.begin(); cIt != choices.end(); ++cIt)
			{
			NxsSaxOutputWrapper cWrapper(outStream, "choice");
			cWrapper << *cIt;
			}
		}
		
	void DistributionCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		const ClassOfDistribution distClass = this->classOfDistribution;
		const std::string s = (distClass & kDiscrete ? (distClass & kContinuous ? "Any" : "Discrete") : "Continuous");
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("distrib_class", s));
		if (this->nVariates != 0)
			{
			std::string nv;
			nv << this->nVariates;
			attr.push_back(NxsSAXAttribute("num_variates", nv));
			}
		NxsSaxOutputWrapper tiWrapper(outStream, "distribution_type_info", attr);
			// currently ncl does not support min and max other than 0 and 1, so this code doesn't handle 
			// all of the elements in the command_state schema (e.g. we don't worry about min_val, max_val) - we never get "Bounded" or "Unbounded" range constraints
		std::string range_constraint_enum;
		if (distClass == kDiscreteOrContinuous || distClass == kDiscrete || distClass == kContinuous)
			range_constraint_enum = "Any"; //@ ncl doesn't discriminate between accepting distributions with any range and _requiring_ unbounded distribution
		else 
			{   // order is important in this if/else (e.g. SumToOne overrides ZeroToOne)
			if (distClass & kSumToOneConstraintBit)
				range_constraint_enum = "SumToOne";
			else if (distClass & kZeroToOneBit)
				range_constraint_enum = "ZeroToOne";
			else if (distClass & kNonNegativeBit)
				range_constraint_enum = "NonNegative";
			else
				{
				NXS_ASSERT(0);
				range_constraint_enum = "Any";
				}
			}
		VecNxsSAXAttribute rcAttr(1, NxsSAXAttribute("constraint", range_constraint_enum));
		NxsSaxOutputWrapper rcWrapper(outStream, "range_constraint", rcAttr);
		}
		
	void NxsInFileCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "infile_type_info", attr);
		}
	void NxsOutFileCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "outfile_type_info", attr);
		}
	
	void NxsOutputCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		NxsSaxOutputWrapper tiWrapper(outStream, "output_type_info");
		if (true)
			{
			VecNxsSAXAttribute attr(1, NxsSAXAttribute("suppress", currentValue->IsSuppressOutput() ? "true" : "false"));
			NxsSaxOutputWrapper tiWrapper(outStream, "default", attr);
			if (currentValue->IsRedirectToNCLStream())
				{
				for (UInt i = 0; i < currentValue->redirectStreams.size();++i)
					{
					NxsSaxOutputWrapper tiWrapper(outStream, "redirect");
					outStream << currentValue->GetRedirectionStreamName(i);
					}
				}
			if (currentValue->IsRedirectToFile())
				{
				// will be a for loop if we allow multiple files
				attr[0] = NxsSAXAttribute("append", currentValue->filePath.GetAppend() ? "true" : "false");
				attr.push_back(NxsSAXAttribute("replace", currentValue->filePath.GetReplace() ? "true" : "false"));
				attr.push_back(NxsSAXAttribute("path", currentValue->filePath.GetFullName()));
				NxsSaxOutputWrapper tiWrapper(outStream, "file", attr);
				}
			}
		}
		
	void NxsMixedCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream&) const
		{
		}
	void BoolCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream ) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "bool_type_info", attr);
		}
	void RestrictNameCmdOption::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		RefreshTabooList();
		NxsSaxOutputWrapper tiWrapper(outStream, "restricted_string_type_info");
		for (VecString::const_iterator cIt = tabooList.begin(); cIt != tabooList.end(); ++cIt)
			{
			NxsSaxOutputWrapper cWrapper(outStream, "disallowed_value");
			cWrapper << *cIt;
			}
		}
		
	template <typename T>
	void NumberCmdOption<T>::WriteTypeInfoStateElementPrivate(NxsXMLSocketOutputStream & outStream, const std::string & s) const
		{
		VecNxsSAXAttribute attr(3);
		std::string minVStr, maxVStr;
		minVStr << DescribeValue(kMinValueDescription);
		maxVStr << DescribeValue(kMaxValueDescription);
		attr[0] = NxsSAXAttribute("default", GetCurrentValueAsString());
		attr[1] = NxsSAXAttribute("min_val", minVStr);
		attr[2] = NxsSAXAttribute("max_val", maxVStr);
		NxsSaxOutputWrapper tiWrapper(outStream, s, attr);
		}
	
	template <>
	void NumberCmdOption<UInt>::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		WriteTypeInfoStateElementPrivate(outStream, "integer_type_info");
		}
		
	template <>
	void NumberCmdOption<int>::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		WriteTypeInfoStateElementPrivate(outStream, "integer_type_info");
		}
	template <>
	void NumberCmdOption<double>::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		WriteTypeInfoStateElementPrivate(outStream, "double_type_info");
		}
	template <>
	void SimpleCmdOption<char>::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "char_type_info", attr);
		}

	template <>
	void SimpleCmdOption<std::string>::WriteTypeInfoStateElement(NxsXMLSocketOutputStream & outStream) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "string_type_info", attr);
		}

	void NxsIndexSetCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const
		{
		VecNxsSAXAttribute attr(1, NxsSAXAttribute("default", GetCurrentValueAsString()));
		NxsSaxOutputWrapper tiWrapper(outStream, "set_type_info", attr);
		}

#else //if defined(NCL_SOCKET_IO)

	void NxsChoiceCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream & ) const
		{
		}
		
	void DistributionCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	void NxsInFileCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	void NxsOutFileCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
	void NxsMixedCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	void BoolCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
	
	void RestrictNameCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	template <typename T>
	void NumberCmdOption<T>::WriteTypeInfoStateElementPrivate(NxsHiddenQueryStream &, const std::string &) const
		{
		}
	
	template <>
	void NumberCmdOption<UInt>::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	template <>
	void NumberCmdOption<int>::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	template <>
	void NumberCmdOption<double>::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
		
	template <>
	void SimpleCmdOption<char>::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}

	template <>
	void SimpleCmdOption<std::string>::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}

	void NxsIndexSetCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream &) const
		{
		}
	void NxsOutputCmdOption::WriteTypeInfoStateElement(NxsHiddenQueryStream & ) const
		{
		}

#endif// defined(NCL_SOCKET_IO)

//// end INCLUDE_TO_AVOID_LINK_ERROR hack
//#endif
