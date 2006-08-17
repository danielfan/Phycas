//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/characters/nxs_characters_block.hpp"
#include "pyphy/src/ncl/characters/nxs_characters_manager.hpp"
#include "pyphy/src/ncl/command/nxs_auto_command.hpp"
#include "pyphy/src/ncl/command/nxs_primitive_cmd_param.hpp"
#include <boost/bind.hpp>
#include "pyphy/src/ncl/command/nxs_choice_cmd_param.hpp"
#include "pyphy/src/ncl/misc/nxs_test_impl.hpp"
#include "pyphy/src/ncl/command/nxs_set_cmd_param.hpp"
#include "pyphy/src/ncl/characters/nxs_characters_block_cmd_settings.hpp"
using std::string;

NxsCommandShPtr NxsAlternativeTaxaBlock::GetTaxLabelsCommand()
	{
	NxsCommandShPtr taxLabs = CreateTaxLabelsCommand();
	taxLabs->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	//	taxlabels only allowed if newtaxa is specified
	//
	string failMsg;
	failMsg  << "The DIMENSIONS command (specifying " << GetNumTaxaName() << ") must come before TAXLABELS";
	StoredMsgSource newTaxCheckMsg(failMsg, string());
	boost::function0<bool> newTaxCallBack = boost::bind(&NxsAlternativeTaxaBlock::CreatedNewTaxa, this);
	NxsTestShPtr newTaxCheck = NxsTestShPtr(new BoolFuncTest(newTaxCallBack, newTaxCheckMsg));
	taxLabs->AddCommandAvailabilityReq(newTaxCheck);
	taxLabs->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	
#	if defined(NDEBUG)
		taxLabs->FinishedAddingOptions();
#	else
		bool check = taxLabs->FinishedAddingOptions();
        NXS_ASSERT(check);
#	endif

	return taxLabs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a Dimension command to be added to the CommandManger interface of a NxsBlock derived from this 
|	instantion of NxsAlternativeTaxaBlock.
|	the command will only have a newtax and numtax options (the names of these options come from GetNewTaxaName() and 
|	GetNumTaxaName(), so they don't have to be "newtax")
|	IMPORTANT: AutoCommand::FinishedAddingOptions() must still be called (this allows other DIMENSIONS options to be added)
*/
NxsCommandShPtr NxsAlternativeTaxaBlock::CreateUnfinishedDimensionsCommand(NxsTestShPtr noAdvancedCmdsTest)
	{
	
	//	new taxa default should be false for CHARACTERS block and true for DATA block
	//	Note that NewTaxa is persistent because of its strange default behaviour (it is returned to whatever the correct default value is in Reset and CanReadBlock)
	//
	bool newTaxaDef = GetNewTaxaDefault();
	typedef NxsExecuteCommandCallback<NxsAlternativeTaxaBlock, NxsDimensionsSettings> DimensionsSettingsCmdCallback;
	typedef boost::shared_ptr<NxsDimensionsSettings> DimensionsSettingsPtr;
	DimensionsSettingsPtr settingsStruct = DimensionsSettingsPtr(new NxsDimensionsSettings());
	dimensionSettings = settingsStruct;	// NxsAlternativeTaxaBlock keeps a copy of the SettingStruct pointer
	DimensionsSettingsCmdCallback vDimentionsCmdCallback = DimensionsSettingsCmdCallback(this, &NxsAlternativeTaxaBlock::HandleDimensions, settingsStruct);
	typedef NxsAutoCommand<DimensionsSettingsCmdCallback>  DimensionsCommand;
	typedef boost::shared_ptr<DimensionsCommand> DimensionsCommandPtr;
	NxsCommandShPtr dimen = NxsCommandShPtr(new DimensionsCommand("Dimensions", vDimentionsCmdCallback));
	//													Name			Var							default		Min		Max				 persistent	 permission
	UIntCmdOption * ntaxCI = new UIntCmdOption	(	GetNumTaxaName(), 	&settingsStruct->nTaxa, 	0, 			1, 	UINT_MAX, false,kCmdPermBasicUser);
	BoolCmdOption * newTaxCI = new BoolCmdOption	(	GetNewTaxaName(),	&settingsStruct->newTaxa,	newTaxaDef,				  true, kCmdPermBasicUser); // see note above on persistence of newtax
	string testMsg;
	testMsg << GetNewTaxaName() << " = true was not specified before " << GetNumTaxaName();
	StoredMsgSource newTaxTestMsg(testMsg, string());
	BoolVarTest *barenewTaxTestPtr = new BoolVarTest(NxsVariableWrapper<bool>(&settingsStruct->newTaxa), newTaxTestMsg);
	NxsTestShPtr newTaxTest = NxsTestShPtr(barenewTaxTestPtr);
	ntaxCI->AddRequirement(newTaxTest);
	
	dimen->AddCommandAvailabilityReq(noAdvancedCmdsTest);
	dimen->AddKeyword(NxsCmdOptionShPtr(ntaxCI));
	dimen->AddKeyword(NxsCmdOptionShPtr(newTaxCI));
	return dimen;
	}
	
void NxsCharactersBlock::AddRecognizedCommands()
	{
	commands.clear();
	string testMsg;
	testMsg << "a " << GetAdvancedCommandName() << " command has been read";
	StoredMsgSource advancedTestMsg(testMsg, string());
	boost::function0<bool> advancedCallBack = boost::bind(&NxsAlternativeTaxaBlock::NoAdvancedCommandsHaveBeenRead, this);
	BoolFuncTest *bareAdvancedTestPtr = new BoolFuncTest(NxsFuncWrapper<bool>(advancedCallBack), advancedTestMsg);
	NxsTestShPtr noAdvancedCmdsTest = NxsTestShPtr(bareAdvancedTestPtr);
	
	AddDimensionsCommand(noAdvancedCmdsTest);
	
	
	testMsg.clear();
	testMsg << "the DIMENSIONS command has not been read";
	StoredMsgSource noDimensionTestMsg(testMsg, string());
	boost::function0<bool> noDimensionCallBack = boost::bind(&NxsAlternativeTaxaBlock::DimensionsHasBeenRead, this);
	BoolFuncTest *bareNoDimTestPtr = new BoolFuncTest(NxsFuncWrapper<bool>(noDimensionCallBack), noDimensionTestMsg);
	NxsTestShPtr someCharsExpectedTest = NxsTestShPtr(bareNoDimTestPtr);
	
	AddFormatCommand(noAdvancedCmdsTest);
	AddEliminateCommand(noAdvancedCmdsTest, someCharsExpectedTest);
	//	several commands require that nchar != 0 (this will only occur if the dimensions command wasn't read)  Create a test for this condition
	//
	

	typedef NxsExternalCommandParser<NxsCharactersBlock>	CharsBlockParser;
	typedef NxsManualCommand<CharsBlockParser, NxsDummyExecuteCommandCallback> CharsBlockParsedCommand;
	typedef boost::shared_ptr<CharsBlockParsedCommand> 	 CharsBlockParsedCommandPtr;
	
	CharsBlockParser matrixParser = CharsBlockParser(this, &NxsCharactersBlock::ParseMatrix);
	CharsBlockParsedCommandPtr matrixCmd = CharsBlockParsedCommandPtr(new CharsBlockParsedCommand("Matrix", matrixParser, NxsDummyExecuteCommandCallback()));
	matrixCmd->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	matrixCmd->AddCommandAvailabilityReq(someCharsExpectedTest);
	AddCommand(matrixCmd);
	
	AddCommand(CreateTaxLabelsCommand());
	
	CharsBlockParser vCSLParser = CharsBlockParser(this, &NxsCharactersBlock::ParseCharStateLabels);
	CharsBlockParsedCommandPtr vCSLCmd = CharsBlockParsedCommandPtr(new CharsBlockParsedCommand(GetCharStateLabelsCmdName(), vCSLParser, NxsDummyExecuteCommandCallback()));
	vCSLCmd->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	vCSLCmd->AddCommandAvailabilityReq(someCharsExpectedTest);
	AddCommand(vCSLCmd);
	
	CharsBlockParser vCLParser = CharsBlockParser(this, &NxsCharactersBlock::ParseCharLabels);
	CharsBlockParsedCommandPtr vCLCmd = CharsBlockParsedCommandPtr(new CharsBlockParsedCommand(GetCharLabelsCmdName(), vCLParser, NxsDummyExecuteCommandCallback()));
	vCLCmd->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	vCLCmd->AddCommandAvailabilityReq(someCharsExpectedTest);
	AddCommand(vCLCmd);
	
	
	CharsBlockParser vSLParser = CharsBlockParser(this, &NxsCharactersBlock::ParseStateLabels);
	CharsBlockParsedCommandPtr vSLCmd = CharsBlockParsedCommandPtr(new CharsBlockParsedCommand(GetStateLabelsCmdName(), vSLParser, NxsDummyExecuteCommandCallback()));
	vSLCmd->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	vSLCmd->AddCommandAvailabilityReq(someCharsExpectedTest);
	AddCommand(vSLCmd);
	}

bool NxsCharactersBlock::AddEliminateCommand(NxsTestShPtr noAdvancedCmdsTest, NxsTestShPtr somCharsExpectedTest)
	{
	typedef NxsExecuteCommandCallback<NxsCharactersBlock, NxsEliminateCmdSettings> EliminateSettingsCmdCallback;
	typedef boost::shared_ptr<NxsEliminateCmdSettings> EliminateSettingsPtr;
	EliminateSettingsPtr settingsStruct = EliminateSettingsPtr(new NxsEliminateCmdSettings());
	EliminateSettingsCmdCallback vEliminateCmdCallback = EliminateSettingsCmdCallback(this, &NxsCharactersBlock::HandleEliminate, settingsStruct);
	
	typedef NxsAutoCommand<EliminateSettingsCmdCallback>  EliminateCommand;
	typedef boost::shared_ptr<EliminateCommand> EliminateCommandPtr;
	EliminateCommandPtr eliminate = EliminateCommandPtr(new EliminateCommand("Eliminate", vEliminateCmdCallback));
	NxsIndexSetCmdOption *elimSetCI = InstantiateCharSetOption(string(), &settingsStruct->toEliminate, string(), &charactersMgr, false, kCmdPermBasicUser); 
	eliminate->AddUnnamedSetting(NxsCmdOptionShPtr(elimSetCI));
	eliminate->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	eliminate->AddSettingReq(noAdvancedCmdsTest);
	eliminate->AddSettingReq(somCharsExpectedTest);
	bool check = eliminate->FinishedAddingOptions();
	NXS_ASSERT(check);
	AddCommand(eliminate);
	return check;
	}


bool NxsCharactersBlock::AddDimensionsCommand(NxsTestShPtr noAdvancedCmdsTest)
	{
	bool check = true;
	
	//	new taxa default should be false for CHARACTERS block and true for DATA block
	//	Note that NewTaxa is persistent because of its strange default behaviour (it is returned to whatever the correct
	//	default value is in Reset and CanReadBlock
	//
	NxsCommandShPtr dimen = CreateUnfinishedDimensionsCommand(noAdvancedCmdsTest);
	//											Name			Var										default		Min		Max		  persistent	 permission
	UIntCmdOption *nCharCI = new UIntCmdOption(GetNumCharsName(),  &dimensionSettings->secondDimension, 0, 			1, 		UINT_MAX, false, kCmdPermBasicUser );
	dimen->AddKeyword(NxsCmdOptionShPtr(nCharCI));
	check = dimen->FinishedAddingOptions();
	NXS_ASSERT(check);
	
	AddCommand(dimen);
	dimen->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceAlways);
	return check;
	}
	

bool NxsCharactersBlock::AddFormatCommand(NxsTestShPtr noAdvancedCmdsTest)
	{
	bool check = true;
	
	typedef NxsExecuteCommandCallback<NxsCharactersBlock, NxsFormatCmdSettings> CharFormatSettingsCmdCallback;
	typedef boost::shared_ptr<NxsFormatCmdSettings> CharFormatSettingsPtr;
	CharFormatSettingsPtr settingsStruct = CharFormatSettingsPtr(new NxsFormatCmdSettings());
	CharFormatSettingsCmdCallback vFormatCmdCallback = CharFormatSettingsCmdCallback(this, &NxsCharactersBlock::HandleFormat, settingsStruct);
	
	typedef NxsAutoCommand<CharFormatSettingsCmdCallback>  CharFormatCommand;
	typedef boost::shared_ptr<CharFormatCommand> CharFormatCommandPtr;
	
	CharFormatCommandPtr vFormatCommand = CharFormatCommandPtr(new CharFormatCommand("Format", vFormatCmdCallback));
	//																		Name			Var											default			Min		Max														persistent	 permission
	NxsChoiceCmdOption *dataTypeCI =  new NxsChoiceCmdOption(	"Datatype", 	&settingsStruct->dataTypeIndex, &settingsStruct->dataTypeName,  "STANDARD", 	"STANDARD|DNA|RNA|NUCLEOTIDE|PROTEIN|CONTINUOUS", false, false, kCmdPermBasicUser);
	BoolCmdOption  *resCaseCI =  new BoolCmdOption(		"RespectCase",	&settingsStruct->respectingCase, false, false, kCmdPermBasicUser);
	NxsCharCmdOption *missingCI = new NxsCharCmdOption(		"Missing",		&settingsStruct->missingSymbol, '?', false, false, kCmdPermBasicUser);
	NxsCharCmdOption *gapCI = new NxsCharCmdOption(			"Gap",			&settingsStruct->gapSymbol, '\0', false, false, kCmdPermBasicUser);
	QuotedStringCmdOption *symbolsCI = new QuotedStringCmdOption("Symbols",	 &settingsStruct->rawSymbols, string(), false, kCmdPermBasicUser); 
	QuotedStringCmdOption *equateCI = new QuotedStringCmdOption("Equate",	&settingsStruct->rawEquate,	string(), false, kCmdPermBasicUser); 
	NxsCharCmdOption *matchCharCI = new NxsCharCmdOption(		"MatchChar",	&settingsStruct->matchSymbol, '\0', false, false, kCmdPermBasicUser);
	BoolCmdOption  *labelsCI =  new BoolCmdOption(		"Labels",		&settingsStruct->labels, true, false, kCmdPermBasicUser);
	BoolCmdOption  *transposeCI =  new BoolCmdOption(	"Transpose",	&settingsStruct->transposing, false, false, kCmdPermBasicUser);
	BoolCmdOption  *interleaveCI =  new BoolCmdOption(	"Interleave",	&settingsStruct->interleaving, false, false, kCmdPermBasicUser);
	tokensCI = NxsCmdOptionShPtr(new BoolCmdOption(	"Tokens",		&settingsStruct->tokens, false, false, kCmdPermBasicUser));
	NxsChoiceCmdOption *itemsCI = new NxsChoiceCmdOption	(	"Items", 		&settingsStruct->itemsIndex, &settingsStruct->itemTypeName, "STATES", "STATES|MIN|MAX|MEDIAN|AVERAGE|VARIANCE|STDERROR|SAMPLESIZE",false, false, kCmdPermBasicUser);
	//	REMEMBER TO UPDATE HandleFormat if you change the order of the choicess !!!
	statesFormatCI = NxsCmdOptionShPtr(new NxsChoiceCmdOption(	"StatesFormat", &settingsStruct->stateFormIndex, &settingsStruct->stateFormName, "STATESPRESENT", "STATESPRESENT|INDIVIDUALS|COUNT|FREQUENCY",false, false, kCmdPermBasicUser));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(dataTypeCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(resCaseCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(missingCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(gapCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(symbolsCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(equateCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(matchCharCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(labelsCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(transposeCI));
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(interleaveCI));
	vFormatCommand->AddKeyword(tokensCI);
	vFormatCommand->AddKeyword(NxsCmdOptionShPtr(itemsCI));
	vFormatCommand->AddKeyword(statesFormatCI);
	// still need to check legality of missing, gap and match
	//
	check = vFormatCommand->FinishedAddingOptions();
	NXS_ASSERT(check);
	vFormatCommand->SetNumberOfTimesAllowedPerBlock(NxsCommand::kOnceMax);
	vFormatCommand->AddSettingReq(noAdvancedCmdsTest);
	AddCommand(vFormatCommand);
	return check;

	}

