#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_taxa_block.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "ncl/nxs_exception.hpp"
#include "ncl/output/nxs_output.hpp"
using std::string;


/*----------------------------------------------------------------------------------------------------------------------
|	Checks that Dimensions and TaxLabels commands are in agreement.
|	Alerts the Taxa Manager that another taxa block has been read via NewBlockRead.  
*/
CmdResult NxsTaxaBlock::EndEncountered()
	{
	return taxaMgr.NewBlockRead(this);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Creates command objects for reading the DIMENSIONS and TAXLABELS commands.
*/
void NxsTaxaBlock::AddRecognizedCommands()	
	{
	typedef NxsExternalCommandParser<NxsTaxaBlock>							TaxaBlockParser;
	typedef NxsManualCommand<TaxaBlockParser, NxsDummyExecuteCommandCallback> TaxaBlockNoExecuteParsedCommand;
	typedef boost::shared_ptr<TaxaBlockNoExecuteParsedCommand> 				TaxaBlockNoExecuteParsedCommandPtr;
	
	TaxaBlockParser dimenParser = TaxaBlockParser(this, &NxsTaxaBlock::ParseDimensions);
	TaxaBlockNoExecuteParsedCommandPtr dimenCommand =  TaxaBlockNoExecuteParsedCommandPtr(new TaxaBlockNoExecuteParsedCommand("Dimensions", dimenParser, NxsDummyExecuteCommandCallback()));
	dimenCommand->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceAlways);
	AddCommand(dimenCommand);
	NxsCommandShPtr taxLabelsCommand = CreateTaxLabelsCommand();
	taxLabelsCommand->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceAlways);
	AddCommand(taxLabelsCommand);
	}

NxsCommandShPtr NxsTaxaLabelReader::CreateTaxLabelsCommand()
	{
	typedef NxsExternalCommandParser<NxsTaxaLabelReader> TaxaLabelParser;
	typedef NxsManualCommand<TaxaLabelParser, NxsDummyExecuteCommandCallback> TaxaLabelNoExecuteParsedCommand;
	typedef boost::shared_ptr<TaxaLabelNoExecuteParsedCommand> 				TaxaLabelNoExecuteParsedCommandPtr;
	TaxaLabelParser taxLabelsParser = TaxaLabelParser(this, &NxsTaxaLabelReader::ParseTaxLabels);
	return TaxaLabelNoExecuteParsedCommandPtr(new TaxaLabelNoExecuteParsedCommand("TaxLabels", taxLabelsParser, NxsDummyExecuteCommandCallback()));
	}

bool NxsTaxaBlock::ParseDimensions(NxsToken &token)	
	{
	++token;
	token.ThrowIfNot("NTAX");
	++token;
	token.ThrowIfNot("=");
	++token;
	unsigned u;
	if (!IsAnUnsigned(token.GetTokenReference(), &u)|| u == UINT_MAX || u == 0)
		token.ThrowUnexpectedTokenNxsException("a positive number");
	assert(taxLabels.empty());
	taxLabels.clear();
	VecString tmp;
	FillVectorWithNumbers(tmp, 1, u + 1);
	taxLabels.AppendLabels(tmp);
	++token;
	token.ThrowIfNot(";");
	return true;
	}

bool NxsTaxaLabelReader::ParseTaxLabels(NxsToken &token)	
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::isdigit;
#	endif
	const unsigned ntax = (unsigned)taxLabels.size();
	if (ntax == 0)
		throw NxsException("The NTAX must be read from the DIMENSIONS command before the TAXLABELS command.", token);
	++token;
	for (unsigned i = 0; i < ntax; ++i)
		{
		const std::string & newLabel = token.GetTokenReference();
		if (!newLabel.empty())
			{
			if (IsLegalNexusWord(newLabel))
				{
				if (isdigit(newLabel.at(0)))
					{
					unsigned u;
					if (IsAnUnsigned(newLabel, &u) && u != (i+1))
						throw NxsException("Numbers are not a legal NEXUS taxon labels.", token);
					}
				const unsigned ind = taxLabels.FindIndex(newLabel, false);
				if (ind != UINT_MAX)
					{
					string s;
					s << "The label " << newLabel << " cannot be used for taxon " << i + i<< " it was already used for taxon "<< ind + 1 << " (error in the TAXLABELS command).";
					throw NxsException(s, token);
					}	
				taxLabels.SetLabel(i, newLabel);
				}
			else
				{
				string s;
				s << newLabel << " is not a legal NEXUS taxon label.";
				throw NxsException(s, token);
				}
			}
		++token;
		}
	token.ThrowIfNot(";");
	labelsRead = true;
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	calls base-class NxsCommandManagerBlock(OutputManager *) constructor.  Initializes the empty taxa block.
*/
NxsTaxaBlock::NxsTaxaBlock(NxsTaxaManager & inTaxaMgr) 
  	: NxsCommandManagerBlock("TAXA"),
	taxaMgr(inTaxaMgr)
  	{
	InitializeRecognizedCommands();
   	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Changes the Label of taxon i to the string s.
*/
void NxsTaxaBlock::ChangeLabel(
  unsigned i, 
  const string &s)
	{
	if (i >= (unsigned) taxLabels.size())
		throw NxsX_IndexNotFound();
	taxLabels.SetLabel(i, s);
	}

	
