//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"
#include "ncl/command/nxs_primitive_cmd_param.hpp"
#include "ncl/command/nxs_command.hpp"
#include "ncl/misc/string_extensions.hpp"

using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::size_t;
#endif

NxsOutputDestinationDescription &NxsOutputDestinationDescription::operator=(const NxsOutputDestinationDescription &other)
	{
	outDestinationType = other.outDestinationType;
	redirectStreams = other.redirectStreams;
	filePath = other.filePath;
	return *this;
	}
NxsOutputDestinationDescription::NxsOutputDestinationDescription()
	:outDestinationType(kSuppressOutput)
	{
	assert(true);
	}
NxsOutputDestinationDescription::NxsOutputDestinationDescription(ValidNCLOutputRedirect o)
	:outDestinationType(kDirectOutputToNCLStream),
	redirectStreams(1,o)
	{
	assert(true);
	}
NxsOutputDestinationDescription::NxsOutputDestinationDescription(const std::string &fn, bool appendIfPresent, bool replaceIfPresent)
	:outDestinationType(kDirectOutputToFile),
	filePath(fn)
	{
	if (appendIfPresent)
		filePath.SetAppend();
	else if (replaceIfPresent)
		filePath.SetReplace();
	}
NxsOutputDestinationDescription::NxsOutputDestinationDescription(ValidNCLOutputRedirect o, const std::string &fn, bool appendIfPresent, bool replaceIfPresent)
	:outDestinationType(kDirectOutputToBoth),
	redirectStreams(1, o),
	filePath(fn)
	{
	if (appendIfPresent)
		filePath.SetAppend();
	else if (replaceIfPresent)
		filePath.SetReplace();
	}
string NxsOutputCmdOption::GetCurrentValueAsString() const
	{
	return currentValue->GetNexusDescription();
	}
	
string NxsOutputDestinationDescription::GetNexusDescription() const
	{
	string s = "(";
	if (outDestinationType == kSuppressOutput)  
		s << "Suppress"; 
	else
		{
		if (outDestinationType == kDirectOutputToNCLStream || outDestinationType == kDirectOutputToBoth)
			{
			for (VecRedirection::const_iterator rsIt = redirectStreams.begin(); rsIt != redirectStreams.end(); ++rsIt)
				s << "RedirectTo = " << TranslateRedirectionStreamName(*rsIt) << ' ';
			}
		if (outDestinationType == kDirectOutputToFile || outDestinationType == kDirectOutputToBoth)
			{
			s << " File = " << filePath.GetFullName();
			if (filePath.GetAppend())
				s << " Append";
			else if (filePath.GetReplace())
				s << " Replace";
			}
		}
	s << ')';
	return s;
	}

bool NxsOutputCmdOption::InterpretRedirectString(const string & s)
	{
	const VecString & allowedRedirectionStreamsNamesVec = NxsOutputDestinationDescription::getLegalRedirectionStreamsNames();
	for (unsigned i = 0; i < allowedRedirectionStreamsNamesVec.size(); ++i)
		{
		if (IsCapAbbreviation(s, allowedRedirectionStreamsNamesVec[i]))
			{
			currentValue->redirectStreams.push_back(static_cast<NxsOutputDestinationDescription::ValidNCLOutputRedirect>(i));
			return true;
			}
		}
	return false;
	}
	
bool NxsOutputCmdOption::ReadValue(
  NxsToken &token, 
  bool equalsAlreadyRead)
	{
	currentValue->redirectStreams.clear();
	if (!equalsAlreadyRead && !EatEqualsThenAdvance(token))
		return false;
	currentValue->outDestinationType = NxsOutputDestinationDescription::kSuppressOutput;
	if (IsCapAbbreviation(token.GetTokenReference(), "Suppress"))
		{
		++token;
		return true;
		}
	if (!EatWordThenAdvance(token, "(", ncl::kStringRespectCase))
		return false;
	if (IsCapAbbreviation(token.GetTokenReference(), "Suppress"))
		return AdvanceThenEatWordThenAdvance(token, ")", ncl::kStringRespectCase);
	for (;;)
		{
		if (IsCapAbbreviation(token.GetTokenReference(), "Redirectto"))
			{
			if (currentValue->outDestinationType == NxsOutputDestinationDescription::kSuppressOutput)
				currentValue->outDestinationType = NxsOutputDestinationDescription::kDirectOutputToNCLStream;
			else if (currentValue->outDestinationType == NxsOutputDestinationDescription::kDirectOutputToFile)
				currentValue->outDestinationType = NxsOutputDestinationDescription::kDirectOutputToBoth;
			if (!AdvanceThenEatEqualsThenAdvance(token))
				return false;
			if (!InterpretRedirectString(token.GetTokenReference()))
				{
				const std::string allowedRedirectionStreamsNames = Join(NxsOutputDestinationDescription::getLegalRedirectionStreamsNames(), std::string("|"));
				return FlagError(NxsCmdOption::miss_wd, allowedRedirectionStreamsNames);
				}
			++token;
			}
		else if (IsCapAbbreviation(token.GetTokenReference(), "File"))
			{
			if (currentValue->outDestinationType == NxsOutputDestinationDescription::kSuppressOutput)
				currentValue->outDestinationType = NxsOutputDestinationDescription::kDirectOutputToFile;
			else if (currentValue->outDestinationType == NxsOutputDestinationDescription::kDirectOutputToNCLStream)
				currentValue->outDestinationType = NxsOutputDestinationDescription::kDirectOutputToBoth;
			else if (currentValue->outDestinationType != NxsOutputDestinationDescription::kDirectOutputToBoth)
				return FlagError(NxsCmdOption::miss_wd, ")");
			if (!AdvanceThenEatEqualsThenAdvance(token))
				return false;
			currentValue->filePath.SetReplace(false); // replace is not persistent
			NxsOutFilePath newlyRead(token.GetTokenReference());
			AssignAllFieldsExceptAppAndReplace(&currentValue->filePath, newlyRead); 
			if (!WasValidRead())
				return false;
			++token;
			currentValue->filePath.sourceOfName = NxsFilePath::kUserSupplied;
			if (IsCapAbbreviation(token.GetTokenReference(), "Append"))
				{
				currentValue->filePath.SetAppend();
				++token;
				}
			else if (IsCapAbbreviation(token.GetTokenReference(), "Replace"))
				{
				currentValue->filePath.SetReplace();
				currentValue->filePath.SetAppend(false);
				++token;
				}
			}
		if (EatWordThenAdvance(token, ")", ncl::kStringRespectCase))
			{
			if (currentValue->outDestinationType == NxsOutputDestinationDescription::kSuppressOutput)
				return FlagError(NxsCmdOption::miss_wd, "Suppress|RedirectTo|File");
			return true;
			}
		else if (!IsCapAbbreviation(token.GetTokenReference(), "File") && !IsCapAbbreviation(token.GetTokenReference(), "RedirectTo"))
			{
			if (currentValue->outDestinationType == NxsOutputDestinationDescription::kSuppressOutput)
				return FlagError(NxsCmdOption::miss_wd, "Suppress|RedirectTo|File");
			else
				return FlagError(NxsCmdOption::miss_wd, "RedirectTo|File");
			}
		else
			ClearErrorState();
		}
	}
	
void NxsOutputCmdOption::StorePreviousValue()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	previousSourceOfSetting = (unsigned) currentValue->filePath.sourceOfName;
	valBefCmd = *currentValue;
	}
	
void NxsOutputCmdOption::ReturnValueToDefault()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	*currentValue = valBefCmd;
	currentValue->filePath.sourceOfName = NxsFilePath::kDefaultOpt;
	}

void NxsOutputCmdOption::RevertValueBecauseCommandFailed() //need to use AssignAllFieldsExceptAppAndReplace 
	{
	*currentValue = valBefCmd;
	currentValue->filePath.sourceOfName = (previousSourceOfSetting == NxsFilePath::kDefaultOpt ? NxsFilePath::kDefaultOpt : NxsFilePath::kUserSupplied);
	}

bool NxsOutputCmdOption::IsCurrentlyValid()
	{
	if (currentValue->outDestinationType == NxsOutputDestinationDescription::kDirectOutputToFile && !currentValue->filePath.IsLegalFileName())
		return FlagError(query_for_error, MakeStrPrintF("The file %s is not a legal file name", currentValue->filePath.GetFullName().c_str()));
	return true;
	}
	
string NxsOutputCmdOption::GetDisplayType(bool includeIndefiniteArticle,  bool plural) const
	{
	return StandardNoun("file or stream name",includeIndefiniteArticle, plural);
	}


string  NxsOutFileCmdOption::GetHandlersChangesToCmd()
	{
	if (currentValue->userCorrection == NxsOutFilePath::kNoCorrection)
		return string();
	string retStr;
	retStr << prefixToOptions << "Replace = " << (currentValue->userCorrection == NxsOutFilePath::kCorrectedToReplace ? "true" : "false");
	retStr << ' ' << prefixToOptions << "Append = " << (currentValue->userCorrection == NxsOutFilePath::kCorrectedToAppend ? "true" : "false");
	currentValue->userCorrection = NxsOutFilePath::kNoCorrection;
	return retStr;
	}
			
NxsOutFileCmdOption::NxsOutFileCmdOption(
  const string & name, 
  const string & prefix, 
  NxsOutFilePath *manip, 
  NxsOutFilePath def,
  NxsCommand *cmd, 
  bool persis, 
  CmdPermissionLevel perm)
  	:SimpleCmdOptionInterface<NxsOutFilePath>(name, manip, def, true, persis, perm),
  	prefixToOptions(prefix)
  	{
  	ReturnValueToDefault();
  	if (cmd != NULL)
  		{
  		BoolCmdOption *replaceCmdOpt = new BoolCmdOption((prefixToOptions + "Replace"),
  										(bool *) &currentValue->replaceFile,
  										false, /*default */
  										false, /* not persistent*/
  										perm);
  		cmd->AddKeyword(NxsCmdOptionShPtr(replaceCmdOpt));
  		BoolCmdOption *appendCmdOpt = new BoolCmdOption((prefixToOptions + "Append"),
  										(bool *) &currentValue->appendToFile,
  										false,	/*default */
  										true,	/*persistent*/
  										perm);
  		cmd->AddKeyword(NxsCmdOptionShPtr(appendCmdOpt));
  		}
	
  	}
  
void AssignAllFieldsExceptAppAndReplace(NxsOutFilePath *leftOP, const NxsOutFilePath &rightOp)
 	{
 	bool prevAppend, prevReplace; // we don't want to return the append and replace fields to their defaults (they have their own cmdoptions)
	prevAppend = leftOP->appendToFile;
	prevReplace = leftOP->replaceFile;
	*leftOP = rightOp;
	leftOP->appendToFile = prevAppend;
	leftOP->replaceFile = prevReplace;
	}
  		
NxsInFileCmdOption::NxsInFileCmdOption(
  const string &n, 
  NxsInFilePath *manip, 
  NxsInFilePath def,
  bool persis, 
  CmdPermissionLevel perm)
  	:SimpleCmdOptionInterface<NxsInFilePath>(n, manip, def, true,  persis, perm)
  	{
  	ReturnValueToDefault();
  	}

bool NxsOutFileCmdOption::ReadValue(
  NxsToken &token, 
  bool equalsAlreadyRead)
  	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		NxsOutFilePath newlyRead(token.GetTokenReference());
		AssignAllFieldsExceptAppAndReplace(currentValue, newlyRead); 
		if (WasValidRead())
			{
			++token;
			currentValue->sourceOfName = NxsFilePath::kUserSupplied;
			return true;
			}
		}
	return false;
  	}

void NxsOutFileCmdOption::ReturnValueToDefault()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	AssignAllFieldsExceptAppAndReplace(currentValue, defaultValue);
	currentValue->sourceOfName = NxsFilePath::kDefaultOpt;
	}

void NxsOutFileCmdOption::StorePreviousValue()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	previousSourceOfSetting = (unsigned) currentValue->sourceOfName;
	valBefCmd = *currentValue;
	}

void NxsOutFileCmdOption::RevertValueBecauseCommandFailed() //need to use AssignAllFieldsExceptAppAndReplace 
	{
	AssignAllFieldsExceptAppAndReplace(currentValue, valBefCmd);
	if (previousSourceOfSetting == NxsFilePath::kDefaultOpt)
		currentValue->sourceOfName = NxsFilePath::kDefaultOpt;
	else
		currentValue->sourceOfName = NxsFilePath::kUserSupplied;
	}
  
bool NxsInFileCmdOption::ReadValue(
  NxsToken &token, 
  bool equalsAlreadyRead)
  	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		*currentValue = NxsInFilePath(token.GetTokenReference());
		if (WasValidRead())
			{
			++token;
			currentValue->sourceOfName = NxsFilePath::kUserSupplied;
			return true;
			}
		}
	return false;
  	}
  
bool NxsOutFileCmdOption::IsCurrentlyValid()
	{
	return (currentValue->IsLegalFileName() ? true : FlagError(query_for_error, MakeStrPrintF("The file %s is not a legal file name", currentValue->GetFullName().c_str())));
	}
	
bool NxsInFileCmdOption::IsCurrentlyValid()
	{ // checking for existence of a file anytime other than right before use is dicey.  we'll just verify that the name is legal here.
	return (currentValue->IsLegalFileName() ? true : FlagError(query_for_error, MakeStrPrintF("The file %s is not a legal file name", currentValue->GetFullName().c_str())));
	}
	
string NxsOutFileCmdOption::GetDisplayType(bool includeIndefiniteArticle,  bool plural) const
	{
	return StandardNoun("file name",includeIndefiniteArticle, plural);
	}
	
string NxsInFileCmdOption::GetDisplayType(bool includeIndefiniteArticle,  bool plural) const
	{
	return StandardNoun("file name",includeIndefiniteArticle, plural);
	}

void NxsInFileCmdOption::ReturnValueToDefault()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	*currentValue = defaultValue;
	currentValue->sourceOfName = NxsFilePath::kDefaultOpt;
	}

void NxsInFileCmdOption::StorePreviousValue()	//need to use AssignAllFieldsExceptAppAndReplace
	{
	previousSourceOfSetting = (unsigned) currentValue->sourceOfName;
	valBefCmd = *currentValue;
	}

void NxsInFileCmdOption::RevertValueBecauseCommandFailed() //need to use AssignAllFieldsExceptAppAndReplace 
	{
	*currentValue = valBefCmd;
	if (previousSourceOfSetting == NxsFilePath::kDefaultOpt)
		currentValue->sourceOfName = NxsFilePath::kDefaultOpt;
	else
		currentValue->sourceOfName = NxsFilePath::kUserSupplied;
	}
  
	
