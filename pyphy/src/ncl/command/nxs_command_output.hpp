#if !defined(NXS_COMMAND_OUTPUT_HPP)
#define NXS_COMMAND_OUTPUT_HPP

#include "phypy/src/ncl/output/nxs_output.hpp"
#include "phypy/src/ncl/output/temp_basic_output_operators.hpp"
#include "phypy/src/ncl/command/nxs_command_manager.hpp"
#include "phypy/src/ncl/command/nxs_command.hpp"
#include "phypy/src/ncl/command/nxs_cmd_param.hpp"
#include "phypy/src/ncl/command/nxs_primitive_cmd_param.hpp"
#include "phypy/src/ncl/nxs_basic_manager.hpp"

#if defined(NCL_SOCKET_IO)
typedef std::pair<UInt, const NxsCmdOption *> PairPlacementParam;
template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, PairPlacementParam, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PairPlacementParam & ppp)
			{
			const NxsCmdOption & cmd_param = *ppp.second;
			VecNxsSAXAttribute atts(2);
			const std::string & n = cmd_param.GetName();
			atts[0] = NxsSAXAttribute("label", n);
			const bool isAvailable = cmd_param.IsAvailable();
			atts[1] = NxsSAXAttribute("available", (isAvailable ? "true" : "false"));
			if (n.empty())
				{
				std::string s;
				s << ppp.first;
				atts.push_back(NxsSAXAttribute("placement", s));
				}
			NxsSaxOutputWrapper nsow(outStream, "cmd_param", atts);
			if (isAvailable)
				cmd_param.WriteTypeInfoStateElement(outStream);
			}
	};
	
template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, NxsCommand, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const NxsCommand & cmd)
			{
			VecNxsSAXAttribute atts(2);
			atts[0] = NxsSAXAttribute("label", cmd.GetName());
			const bool isAvailable = cmd.IsAvailable();
			atts[1] = NxsSAXAttribute("available", (isAvailable ? "true" : "false"));
			NxsSaxOutputWrapper nsow(outStream, "command", atts);
			if (isAvailable)
				{
				unsigned i = 1;
				const VecNxsCmdOptionShPtr & p = cmd.GetUnnamedSettings();
				for (VecNxsCmdOptionShPtr::const_iterator pIt = p.begin(); pIt != p.end(); ++pIt)
					{
					if (!(*pIt)->IsParserSpecificOption())
						{
						Emit<FORMAT_HINT>(outStream, PairPlacementParam(i, pIt->get())) << '\n';
						++i;
						}
					}
				const VecNxsCmdOptionShPtr & k = cmd.GetKeywords();
				for (VecNxsCmdOptionShPtr::const_iterator kIt = k.begin(); kIt != k.end(); ++kIt)
					Emit<FORMAT_HINT>(outStream, PairPlacementParam(i, kIt->get())) << '\n';
				}
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kCommandStateKnownSetOutStyle, NxsIndexSet, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const NxsIndexSet & s)
			{
			VecNxsSAXAttribute attrs(1, NxsSAXAttribute("label", s.GetName()));
			NxsSaxOutputWrapper setWrapper(outStream, "known_set", attrs);
			if (true)   //scope for NxsSaxOutputWrapper construct, desctruct
				{
				NxsSaxOutputWrapper setWrapper(outStream, "members");
				std::string stringToPrint;
				outStream << s.GetNexusDescription();
				/*
				NxsIndexSet::const_iterator sIt = s.begin();
				if (sIt != s.end())
					{
					outStream << (*sIt + 1);
					for (sIt++; sIt != s.end(); ++sIt)
						outStream << ',' << (*sIt + 1);
					} */
				}
			}
	};
						
template<class OUT_STREAM>
class GenericPrinterClass<kCommandStateOutSyle, NxsCommandManager::IndexManagerInfo, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const NxsCommandManager::IndexManagerInfo & stringBasicListMgr)
			{
			NxsBasicListManager * blMgr = stringBasicListMgr.second; 
			if (blMgr != NULL)
				{
				UInt maxNumber = blMgr->GetSize();
				std::string n;
				n << maxNumber;
				VecNxsSAXAttribute attrs(1, NxsSAXAttribute("max_index", n));
				NxsSaxOutputWrapper tiWrapper(outStream, stringBasicListMgr.first , attrs);
				attrs.push_back(NxsSAXAttribute("index", "0"));
				for (unsigned i = 0; i < maxNumber; ++i)
					{
					const std::string & usl = blMgr->GetUserSuppliedLabel(i);
					if (!usl.empty())
						{
						n.clear();
						n << (i+1);
						attrs[0] = NxsSAXAttribute("label", usl);
						attrs[1] = NxsSAXAttribute("index", n);
						NxsSaxOutputWrapper labelWrapper(outStream, "index_label",  attrs);
						}
					}
				VecString sn = blMgr->GetSetNames();
				for (VecString::const_iterator snIt = sn.begin(); snIt != sn.end(); ++snIt)
					{
					const NxsIndexSet * s = blMgr->GetSet(*snIt);
					if (s != NULL)
						Emit<kCommandStateKnownSetOutStyle>(outStream, *s);
					}
				}
			}
	};

template<>
class GenericPrinterClass<kCommandStateOutSyle, NxsCommandManager, NxsXMLSocketOutputStream>
	{
	public:
		GenericPrinterClass(NxsXMLSocketOutputStream & outStream, const NxsCommandManager & cmdMgr)
			{
			NxsSaxOutputWrapper nsow(outStream, "command_state");
			const VecNxsCommandShPtr & cmds = cmdMgr.GetAllCommands();
			for (VecNxsCommandShPtr::const_iterator cIt = cmds.begin(); cIt != cmds.end(); ++cIt)
				Emit<kCommandStateOutSyle>(outStream, **cIt) << '\n';
			const NxsCommandManager::VecIndexManagers & im = cmdMgr.GetCommandEnvironments();
			for (NxsCommandManager::VecIndexManagers::const_iterator imIt = im.begin(); imIt != im.end(); ++imIt)
				Emit<kCommandStateOutSyle>(outStream, *imIt);
			outStream << '\n';
			}
	};
#endif

template<class OUT_STREAM>
class GenericPrinterClass<kCommandStateOutSyle, NxsCommandManager, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & , const NxsCommandManager & )
			{
			// NxsCommandManager prints its state to the NxsXMLSocketOutputStream only
			
			}
	};

#endif

