#ifndef NCL_NXS_CMD_OPTION_H
#define NCL_NXS_CMD_OPTION_H
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/nxs_token.hpp"
#include "phypy/src/ncl/command/nxs_command_decl.hpp"
#include "phypy/src/ncl/misc/nxs_test.hpp"

std::string UserLevelTooLowString(int requiredUserLevel);
/*--------------------------------------------------------------------------------------------------------------------------
|	Abstract base class that encapsulates an option that may follow a command.  It stores the name of option (and its minimal 
|	abbreviation), info on when the option was read, its current state, and functions for reading the command froma a stream
|	of NexusTokens and interacting with objects of type Test so that command settings can be checked.
|	Each NxsAutomaticallyParsedCommand has a collection of CmdOptions and uses them to parse its command lines.
|	To work correctly the following order of function calls must be used (the NxsAutomaticallyParsedCommand class
|	implements this protocol for interacting with command Options).
|	For unnamed commands:
|	1	PrepareToRead()
|	2	if (HasPermission())
|	3		TryToRead(tokenStream, true)
|
|	For named commands:
|	1	PrepareToRead()
|	2	if (CanReadKeyword())
|	3		TryToRead(tokenStream, false)
|
|	NxsCmdOption::TryToRead makes sure the option hasn't been read already and then calls ReadValue.
|	
|	minimally a derived class must override:
|	ReadValue() for interpretting the stream of tokens from the user, 
|		ReadValue should call WasValidRead before exiting if the command is correct (whenever true is returned by read value, 
|		one should be able to assume that WasValidRead() returned true.
|	WasValidRead verifies that the requirements for the command were passed and then calls the virtual IsCurrentlyValid()
|		which should just check the "manipulated variable" field to make sure it is legal.
|	and GetCurrentValueAsString() and GetDisplayType() for providing error/help messages about the command
*/
class NxsCmdOption
	{
	public :
		enum ReadStatus 
			{
			at_factory, 				/* the option has not been changed from the factory default */
			set_before_last_execute, 	/* the user has changed the value, but the command has been executed since that change */
			old_setting, 				/* the user has changed the value since the last time the command was executed, but it wasn't set in the command that is currently being read/executed */
			new_setting					/* the value of the option was changed during the current command */
			};
		enum	CmdOptionErrState {
				no_err 				= 0x0000,
				miss_wd				= 0x0001,
				unrecognized 		= 0x0002,
				too_big				= 0x0004,
				too_small			= 0x0008,
				unrecog_labile 		= 0x0010,
				too_big_labile		= 0x0020,
				too_small_labile	= 0x0040,
				read_twice			= 0x0080,
				illegal_modulus		= 0x0100,
				reserved_word		= 0x0200,
				illegal_range		= 0x0400,
				parse_req_failed	= 0x0800,
				illegal_punctuation = 0x1000,
				unexpected_char	    = 0x2000,
				illegal_name	    = 0x4000,
				query_for_error		= 0x8000,
				cancelled_by_user	= 0x10000
				}; 
		//	Constructor/modifiers
		//
		NxsCmdOption(const std::string &n, bool validityIsLabile, bool persist, CmdPermissionLevel pLevel );	
		virtual ~NxsCmdOption();
		
		//	Modifiers used by AutoCommand
		//
		void				AddAvailabilityReq(NxsTestShPtr cr);
		void				AddRequirement(NxsTestShPtr req);
		virtual bool		AbbreviateChoices();
		void 				CommandWasExecuted();
		void				PrepareToRead();
		void				SetAbbreviation(const std::string &s);
		virtual void		SetConverseAbbreviation(const std::string &);
		virtual	bool		ReadValue(NxsToken &token, bool equalsAlreadyRead = false) = 0;
		bool 				ReadStringAsValue(std::string s);
		void				ReturnToDefault();
		void				RevertBecauseCommandFailed();
		void				SetDescription(const std::string &d);
		void				TryToRead(NxsToken &token, bool equalsAlreadyRead);
		
		//	Accessors
		//
		virtual bool		CanReadKeyword(const std::string &s, int permLevel) const;
		bool				CheckIfOld() const;
		const std::string & GetAbbreviation() const;
		const std::string & GetName() const;
		const std::string & GetDescription() const;
		const std::string & GetErrorSnippet();
		virtual std::string	GetConverseName() const;
		virtual std::string GetCurrentValueAsString() const = 0;
		virtual std::string	GetHandlersChangesToCmd() {return std::string();}
		virtual VecString   GetLegalChoices() const;
		virtual std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const = 0;
		int					GetErrorState() const;
		int 				GetPermissionLevel() const;
		ReadStatus			GetReadStatus() const;
		virtual VecString   GetValidArgument() = 0;
		bool				HadError() const;
		bool				HasPermission(int) const;
		bool				HasBeenRead();
		bool				WasValidRead();
		std::string			GetReasonForUnavailability() const;
		bool				IsAvailable() const;
		bool				IsPersistent() const;
		virtual bool		IsCurrentlyValid() = 0;
		bool				WasErrorFatal(bool commandContext) const;
		virtual bool		IsParserSpecificOption() const {return false;} //@temp hack for detecting RequiredToken command "parameters"
		virtual void		WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const = 0;
		bool				operator==(const std::string &) const;
		
		int					errState;
		
	protected:	
		
		virtual void		ReturnValueToDefault() = 0;
		virtual void		RevertValueBecauseCommandFailed() = 0;
		virtual void		StorePreviousValue() = 0;

		bool				AdvanceThenEatEquals(NxsToken &);
		bool				AdvanceThenEatEqualsThenAdvance(NxsToken &);
		bool				AdvanceThenEatWord(NxsToken &, const std::string &, ncl::StringCmpEnum );
		bool				AdvanceThenEatWordThenAdvance(NxsToken &, const std::string &, ncl::StringCmpEnum );
		void				ClearErrorState();
		bool				EatEquals(NxsToken &);
		bool				EatEqualsThenAdvance(NxsToken &);
		bool				EatWord(NxsToken &, const std::string &, ncl::StringCmpEnum );
		bool				EatWordThenAdvance(NxsToken &, const std::string &, ncl::StringCmpEnum );
		bool				FlagError(int );
		bool 				FlagError(int  newError, const std::string &erS);
		void				ResetCheckIfOld(bool);
		
		std::string			name;
		std::string			abbreviation;
		std::string			description;	
		std::string			errSnippet;	
		int					userRestriction;
		ReadStatus 			status;
		ReadStatus			statusBefCmd;
		mutable std::vector<NxsTestShPtr>	parseLevelRequirements;
		mutable std::vector<NxsTestShPtr>	availableRequirements;
	private :
		enum OptionFlagBits
			{
			kNotPersitent 	  = 0x00,
			kPersistent 	  = 0x01,
			kCheckIfOld 	  = 0x02
			};
		char				optionFlag;
		
		
		NxsCmdOption(const NxsCmdOption&);			//never use	- don't define.  Pointer to the manipulated values (in derived classes) precludes simple copy by value.
		NxsCmdOption & operator=(const NxsCmdOption&);			//never use	- don't define.  Pointer to the manipulated values (in derived classes) precludes simple copy by value.
	};
/*---------------------------------------------------------------------------------------------------
|	Adds a test that must be passed for the command option to be available for reading
*/
inline void  NxsCmdOption::AddAvailabilityReq(NxsTestShPtr cr)  
	{
	availableRequirements.push_back(cr);
	}
		
class NxsRequiredToken : public NxsCmdOption
	{
	public:
		NxsRequiredToken(const std::string &rt);
		
		VecString		GetValidArgument()
			{
			return VecString(1, reqToken);
			}
		bool 		ReadValue(NxsToken &token, bool equalsAlreadyRead = false);
		std::string GetCurrentValueAsString() const;
		std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
		void		ReturnValueToDefault(){}
		void		RevertValueBecauseCommandFailed(){}
		void		StorePreviousValue() {}
		bool		IsCurrentlyValid() {return true;}
		bool		IsParserSpecificOption() const {return true;}
		void		WriteTypeInfoStateElement(NxsHiddenQueryStream & ) const
				{ // nothing to be done, required tokens are not  command parameter from any perspective other than the parser
				} 
	protected:
		std::string reqToken;
	};
	
inline std::string NxsRequiredToken::GetCurrentValueAsString() const 
	{
	return reqToken;
	}
inline std::string NxsRequiredToken::GetDisplayType(bool , bool) const 
	{
	return reqToken;
	}



inline NxsCmdOption::ReadStatus NxsCmdOption::GetReadStatus() const
	{
	return status;
	}

inline bool NxsCmdOption::CheckIfOld() const
	{
	return (bool) ((optionFlag & kCheckIfOld) != 0);
	}

inline bool NxsCmdOption::IsPersistent() const
	{
	return (bool) ((optionFlag & kPersistent) != 0);
	}

inline void NxsCmdOption::SetDescription(
  const std::string &d)
	{
	description = d;
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template <typename T> class SimpleCmdOptionInterface : public NxsCmdOption
	{
	public :
		virtual VecString		GetValidArgument()
			{
			return VecString();
			}
		virtual std::string GetCurrentValueAsString() const = 0;
		virtual std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const = 0;
		T					GetValue() const {return *currentValue;}
		virtual bool		IsCurrentlyValid() {return true;} 
		virtual void 		StorePreviousValue();
		virtual void		RevertValueBecauseCommandFailed();
		virtual bool		ReadValue( NxsToken &token, bool equalsAlreadyRead) = 0;
  		virtual void		ReturnValueToDefault();
		void				SetValue(const T &newVal) {*currentValue = newVal;}
	
		SimpleCmdOptionInterface(const std::string &n, T *manipVal, T def, bool validityIsLabile, bool persist, CmdPermissionLevel pLevel );
	
	protected:
		
		T	   *currentValue;
		T		defaultValue;
		T		valBefCmd;
	};

template <typename T> class SimpleCmdOption : public SimpleCmdOptionInterface<T>
	{
	public :
		 std::string	GetCurrentValueAsString() const;
		 std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
		virtual bool	ReadValue( NxsToken &token, bool equalsAlreadyRead);
		void			WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
 		SimpleCmdOption(const std::string &n, T *manipVal, T def, bool validityIsLabile, bool persist, CmdPermissionLevel pLevel );
	};


/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<typename T> SimpleCmdOption<T>::SimpleCmdOption(
  const std::string &n, 
  T *manipVal, 
  T def, 
  bool validityIsLabile,
  bool persist, 
  CmdPermissionLevel pLevel)
	:SimpleCmdOptionInterface<T>(n, manipVal, def, validityIsLabile, persist, pLevel)
	{
	}

//	simple function to deal with commonality in GetDisplayType
//
std::string StandardNoun(std::string retStr, bool includeIndefiniteArticle, bool plural); 

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<typename T> SimpleCmdOptionInterface<T>::SimpleCmdOptionInterface(
  const std::string &n, 
  T *manipVal, 
  T def, 
  bool validityIsLabile, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:NxsCmdOption(n, validityIsLabile, persist, pLevel),
  	currentValue(manipVal),
  	defaultValue(def)
  	{
valBefCmd = *currentValue = defaultValue;
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<typename T> inline void SimpleCmdOptionInterface<T>::ReturnValueToDefault()
	{
	*currentValue = defaultValue;
	}
		
/*--------------------------------------------------------------------------------------------------------------------------
|		
*/
template<typename T> inline void SimpleCmdOptionInterface<T>::StorePreviousValue()
	{
	valBefCmd = *currentValue;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<typename T> 
inline void SimpleCmdOptionInterface<T>::RevertValueBecauseCommandFailed()
	{
	*currentValue = valBefCmd;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Adds a test to be checked every time the command is read.  If the test fails an error message will be generated 
|	assuming that the NxsCmdOption is being used by an NxsAutomaticallyParsedCommand)
*/
inline void NxsCmdOption::AddRequirement(NxsTestShPtr req) 	
  {
  parseLevelRequirements.push_back(req);
  }
  
/*--------------------------------------------------------------------------------------------------------------------------
|	returns abbreviation string 
*/
inline const std::string &NxsCmdOption::GetAbbreviation() const 	
  {
  return abbreviation;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline const std::string &NxsCmdOption::GetName() const 	
  {
  return name;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline const std::string &NxsCmdOption::GetDescription() const	
  {
  return description;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string NxsCmdOption::GetConverseName() const		
  {
  return std::string();
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline VecString NxsCmdOption::GetLegalChoices() const 	
	{
	return VecString();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline int NxsCmdOption::GetErrorState() const 	
  {
  return errState;
  }
inline void NxsCmdOption::ClearErrorState() 
	{
	errState = no_err;
	errSnippet.clear();
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::AbbreviateChoices() 	
  {
  return true;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::HasBeenRead()	
  {
  return (status == new_setting || status == old_setting);
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline void NxsCmdOption::SetAbbreviation(const std::string &s)	
  {
  abbreviation = s;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline const std::string &NxsCmdOption::GetErrorSnippet()	
  {
  return errSnippet;
  }

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::WasErrorFatal(
  bool isExecuting) const
  	{
  	if (isExecuting)
  		return HadError();
  	//	error codes marked as labile aren't fatal (unless the command will be executed immediately)
  	//
  	if (HadError())
		return !(GetErrorState() == unrecog_labile || GetErrorState() ==too_big_labile || GetErrorState() == too_small_labile);
	return false;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline int NxsCmdOption::GetPermissionLevel() const
	{
	return userRestriction;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::HadError() const
	{
	return (GetErrorState() != no_err);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Base class does nothing. This is just here so that Bool cmd option can override (pretty clumsy, but not time-critical)
*/
inline void NxsCmdOption::SetConverseAbbreviation(const std::string &)
	{
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::HasPermission(
  int usersLevel) const
  	{
  	return (GetPermissionLevel() <= usersLevel);
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline void NxsCmdOption::ReturnToDefault()
	{
	status = at_factory; //@@ need to detect if the user has changed the default 
	ReturnValueToDefault();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline void NxsCmdOption::RevertValueBecauseCommandFailed()
	{
	status = statusBefCmd;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline 	void NxsCmdOption::TryToRead(
  NxsToken &token,
  bool		 equalsRead)
 	{
 	if (status == NxsCmdOption::new_setting)
		FlagError(read_twice);
	else
		{
		if (ReadValue(token, equalsRead))
			status = new_setting;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	returns false so that one can write "return FlagError(x);"
*/	
inline 	bool NxsCmdOption::FlagError(int newError)
	{
	errState |= newError;
	return false; 
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	replaces previous errSnippet with erS and calls FlagError
|	returns false so that one can write "return FlagError(x, s);"
*/
inline 	bool NxsCmdOption::FlagError(int newError, const std::string &erS)
	{
	errSnippet = erS;
	return FlagError(newError);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::operator==(
  const std::string &r) const
  	{
  	return EqualsCaseInsensitive(GetCurrentValueAsString(), r);
  	}

inline bool NxsCmdOption::EatWord(
  NxsToken &token,
  const std::string &word,
  ncl::StringCmpEnum m)
	{
	if (!StrEquals(token.GetTokenReference(), word,m))
		{
		FlagError(NxsCmdOption::miss_wd,word);
		return false;
		}
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Case
*/
inline bool NxsCmdOption::EatWordThenAdvance(
  NxsToken &token,
  const std::string &word,
  ncl::StringCmpEnum m)
  	{
  	if (!EatWord(token, word, m))
  		return false;
	++token;
	return true;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Case
*/
inline bool NxsCmdOption::AdvanceThenEatWordThenAdvance(
  NxsToken &token,
  const std::string &word,
  ncl::StringCmpEnum m)
  	{
  	++token;
	return EatWordThenAdvance(token,word, m);
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Case
*/
inline bool NxsCmdOption::AdvanceThenEatWord(
  NxsToken &token,
  const std::string &word,
  ncl::StringCmpEnum m)
  	{
  	++token;
	return EatWord(token,word, m);
	}
	

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::EatEquals(
  NxsToken &token)
  	{
	return EatWord(token, "=", ncl::kStringRespectCase);
  	}
  
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::EatEqualsThenAdvance(
  NxsToken &token)
  	{
	return EatWordThenAdvance(token, "=", ncl::kStringRespectCase);
  	}
  
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::AdvanceThenEatEqualsThenAdvance(
  NxsToken &token)
  	{
	return AdvanceThenEatWordThenAdvance(token, "=", ncl::kStringRespectCase);
  	}
  
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsCmdOption::AdvanceThenEatEquals(
  NxsToken &token)
  	{
	return AdvanceThenEatWord(token, "=", ncl::kStringRespectCase);
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<> 
inline std::string SimpleCmdOption<std::string>::GetDisplayType(bool includeIndefiniteArticle,
  bool plural) const
	{
	return StandardNoun("string",includeIndefiniteArticle, plural);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<> 
inline std::string SimpleCmdOption<char>::GetDisplayType(bool includeIndefiniteArticle,
  bool plural) const
	{
	return StandardNoun("character",includeIndefiniteArticle, plural);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline void NxsCmdOption::RevertBecauseCommandFailed()
	{
	status = statusBefCmd ;
	RevertValueBecauseCommandFailed();
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<> 
inline std::string SimpleCmdOption<std::string>::GetCurrentValueAsString() const
	{
	std::string s;
	s.append(*currentValue);
	return s;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<> 
inline std::string SimpleCmdOption<char>::GetCurrentValueAsString() const
	{
	std::string s(1, *currentValue);
	return s;
	}



typedef SimpleCmdOption<std::string> NxsStringCmdOption; 
typedef SimpleCmdOption<char> NxsCharCmdOption; 


#endif
