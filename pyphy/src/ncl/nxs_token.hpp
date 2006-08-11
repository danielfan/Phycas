#ifndef NCL_NXSTOKEN_H
#define NCL_NXSTOKEN_H

#include <boost/noncopyable.hpp>
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/string_extensions.hpp"

// April 2003 changes:
//	Correctly parses embedded comments: test[embedded]comment as one token "testcomment"
//	Parsing 2X-3X faster (on metrowerks compiles, at least)
//	breaks token at ' :  test'string' is now two tokens "test" and "string" previously had been an error
//	comments that are contiguous (or embedded) within a token are now stored until the token is advanced.  This 
//		will be useful if one supports "marking up" tokens with formatting comments
//	To speed up ReadToken several interface changes had to be made:
//		SetLabileFlagBit() has been replaced with AlterTokenReading which is much more restrictive (only affects handling of newlines and hyphens)
//			Note that unlike SetLabileFlagBit, the fields changed by a call to AlterTokenReading remain changed until a subsequent call to function 
//			with the opposite argument turns them off.
//		Instead of setting the labile flag bit saveCommandComments now you must call ReadCommandCommentOrToken()
//		Instead of setting the labile flag bit parentheticalToken, curlyBracketedToken, doubleQuotedToken now you must call ReadAllTokensUntil(char tokenBreaker)
//		Instead of setting the labile flag bit singleCharacterToken now you must call ReadTokenAsSingleCharacter()
//		tildeIsPunctuation, useSpecialPunctuation, preserveUnderscores, ignorePunctuation flags have been deprecated
//	NxsToken is now const-correct
//	If SUPPORT_NXS_TOKEN_LABILE_FLAGS is defined in nxsdefs.h, then NxsToken offers legacy support 
//		to the old (slower, but more flexible) tokenizing that relies on GetNextToken and Labile flags
//		As with all deprecated functions, new code shouldn't use these functions as they may not be supported much longer the functions are:
//		GetNextToken()
//		SetSpecialPunctuationCharacter(char c); 
//		SetLabileFlagBit(int bit);
//
//	Deprecate functions:
//		Begins (this would be easy to bring back) IsAbbreviation is still available
//		BlanksToUnderscores
//		GetTokenAsCStr
//		IsPlusMinusToken
//		IsWhitespaceToken
//		StoppedOn(char ch)  (instead of StoppedOn(x) you can call use x == PeekAtNextChar() )
//		StripWhitespace
//		ToUpper()
//		Write(ostream &out);
//		Writeln(ostream &out);
//	Added GetNumComments, GetNextCommentIndex
//	Also added support for backing up the stream to allow for "reparsing" via Get/SeekTokenizerState
//	

//
/*-------------------------------------------------------------------------------------------------------------------------- 
|	Exception thrown by NxsToken::GetComment
*/
class NxsX_InvalidCommentIndex{};
/*-------------------------------------------------------------------------------------------------------------------------- 
|	Exception thrown by NxsToken::SeekTokenizerState
*/
class NxsX_InvalidTokenizerState{};

/*--------------------------------------------------------------------------------------------------------------------------
| 	This class stores the information about where the token is in its stream.
|	Note that the entire state of the token is not cached in this class, only the file postition fields.
|	To cache up in the token stream (so you can later back up) you need a NxsTokenizerState object.
|	
*/
class NxsTokenPosInfo 
	{
	public :
		
		inline unsigned		GetFileColumn() const;
		inline unsigned  	GetFileLine() const;
		inline file_pos		GetFilePosition() const;
		
	private :
		inline NxsTokenPosInfo(file_pos p, unsigned l, unsigned c, char nextChar, bool a, bool eol);
		
		file_pos 			filepos;	/* file position. the number of characters from the beginning of the stream */
		unsigned       		fileline;	/* the file line (number of newline characters since the beginning of the stream) */
		unsigned       		filecol;	/* the file column (number of characters since the last newline) */
		bool 				atEOF;		/* true if at the end of the file */
		bool				atEOL;		/* true if at the end of a line */
		char 				nextCharInStream;
		
		friend class NxsTokenizerState;
		friend class NxsToken;
		friend class NxsCommand;	
	};

class NxsToken;
/*--------------------------------------------------------------------------------------------------------------------------
| 	used to store NxsToken state (current token and it's place in the stream)
*/
class NxsTokenizerState
	{
	public :
		
		inline NxsTokenPosInfo		GetPosInfo() const;
	
	private :
		NxsTokenizerState(const NxsToken & t);
		NxsTokenizerState();

		NxsTokenPosInfo posInfo;
		std::string	 	token;
#		if defined(NXS_INTERNAL_TOKEN_BUFFER)
			file_pos		posBeforeLastRead;	/* file position. the number of characters from the beginning of the stream */
			unsigned		posInLocalBuffer;
			unsigned		nCharsInLocalBuffer;
			bool			lastReadHitEnd;
#		endif
		friend class NxsToken;	
		friend class NxsCommand;	
	};


/*--------------------------------------------------------------------------------------------------------------------------
| 	Stores the text of a comment, where in the token it was found, and any whether it is an output/command comment
*/
class NxsComment
	{
	public :
		
		inline const 	std::string  &GetCommentText() const;
		inline	 		unsigned	GetLocationInToken() const;
		inline 			bool 		IsCommandComment() const;
		inline 			bool		MatchesModifier(char modChar) const;
		
	private :
		inline NxsComment(const std::string &s, unsigned pos, char comChar, const std::string &progSpecifier);
		
		std::string 	body; 					/* everything that was between the [] except for command codes */
		unsigned	placeInToken;			/* location of the comment in the current token (number of characters from the beginning of the token) */
		char 		modifierChar;			/* the command comment character e.g. !, \ or & */
		std::string 	intendedProgramCode;	/* for && comments intended for only one program not implemented */
	
		friend class NxsToken;
	};

/*--------------------------------------------------------------------------------------------------------------------------
| 	Class used for all input in NCL.  
|	Takes a stream and breaks it into legal nexus tokens through calls to ReadToken() (the prefix ++ operator can also
|		be used to advance the token reader).
|	Access to the current token is through GetToken() or GetTokenReference() (the latter avoids copying the std::string that
|		holds the token and is preferable if a simple read-only check of the token is needed.  WARNING: because it returns the
|		token reference any function call that changes the NxsToken state could update alter this reference).
|	SimpleNexus comments that are encountered are simply stored (and will be flushed when the token position is advanced 
|		again) along with information about where in the current token they were encountered. For access to comments see
|		GetNumComments(), GetNextCommentIndex() and GetComment().
|	Overriding OutputComment(const std::string &) const allows you to display (or ignore) any [! ]- style output comments
|
|	In accordance with the Nexus standard:  
|		_ encountered (in non-quoted strings) are converted to ' '
|		comments are stored but do not affect tokenizing
|		Tokens are broken by whitespace or punctuation
| 		Whitespace is any char < 7, ' ', tab or any newline.
|		Punctuation consist of any of the following ;()]{}/\,:=*\"`+-<>' 
|		All punctuation characters except single-quotes are returned as individual tokens.  
|		When a token begins with a single quote everything between the closing quote is treated as the token (the token will 
|			not contain the surrounding quotes (Note adjacent '' within a single-quoted string are treated as one single-quote).
|			Thus 'test tok'  will be read as one 8 character string and would be identical to test_tok.
|		All three OS newline designations are returned as the character '\n'
|
|	Additional features:
|		To ease reading matrices of single characters use ReadSingleCharacter() instead of ReadToken().  
|		To return command comments as tokens use ReadCommandCommentOrToken() instead of ReadToken().
|		To read and concatenate all tokens until some punctuation mark use ReadAllTokensUntil(char).  For example to read a double-
|			quoted string, after seeing the opening " one could call ReadAllTokensUntil('\"')
|		 AlterTokenReading() can be used to disobey Nexus Tokenizing Rules: Reading interleaved matrices requires recognizing 
|			newline characters as tokens.  Reading numbers is easier if hyphens	are not treated as punctuation.  These can 
|			both be dealt with by calls to AlterTokenReading().
|		End of files are only acceptable in certain places (in between blocks).  Calling SetEOFAllowed(false) will result 
|			in an NxsException exception being thrown if eof is found.  This reduces the need to manually verify that you aren't at eof.
|		Fast reading of streams is implemented through unbuffered io.
|		Implementation of backing up in the stream would be inefficient (and rarely needed).  In those instances in which 
|			one might need to back up to a previous point, use GetTokenizerState() and then SeekTokenizerState().
|
|	For improved code readability/maintainability:
|		IsAbbreviation() and Equals() are provided to allow testing of the current token (these functions help you avoid code 
|			like:  IsCapAbbreviation(nxsToken.GetTokenReference(), "Test")
|		ThrowIfNot(std::string,bool) makes sure that the token equals the argument.  If not an NxsException with an appropriate message 
|			is thrown.
|
|	Implementation notes: at first glance the file position info seems wrong.  this is because the true file position is always one 
|		character ahead and the nextCharInStream is stored.  This make tokenizing much easier. 
|	Several functions that are inlined appear too long to inline.  Given how crucial speed is (all input occurs through NxsToken)
|	the code bloat seems worth ignoring the rules of thumb here (these functions are only calle by a couple of other functions
*/ 
class NxsToken NON_COPYABLE
	{
	public:
			STATELESS_FUNC VecString TokenizeString(const std::string &);
			NxsToken(std::istream & i);
			NxsToken(const std::string & s);
			
#	if defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
			enum NxsTokenFlags	/* For use with the variable labileFlags */
				{
				saveCommandComments		= 0x0001,	/* if set, command comments of the form [&X] are not ignored but are instead saved as regular tokens (without the square brackets, however) */
				parentheticalToken		= 0x0002,	/* if set, and if next character encountered is a left parenthesis, token will include everything up to the matching right parenthesis */
				curlyBracketedToken		= 0x0004,	/* if set, and if next character encountered is a left curly bracket, token will include everything up to the matching right curly bracket */
				doubleQuotedToken		= 0x0008,	/* if set, grabs entire phrase surrounded by double quotes */
				singleCharacterToken	= 0x0010,	/* if set, next non-whitespace character returned as token */
				newlineIsToken			= 0x0020,	/* if set, newline character treated as a token and atEOL set if newline encountered */
				tildeIsPunctuation		= 0x0040,	/* if set, tilde character treated as punctuation and returned as a separate token */
				useSpecialPunctuation	= 0x0080,	/* if set, character specified by the data member special is treated as punctuation and returned as a separate token */
				hyphenNotPunctuation	= 0x0100,	/* if set, the hyphen character is not treated as punctutation (it is normally returned as a separate token) */
				preserveUnderscores		= 0x0200,	/* if set, underscore characters inside tokens are not converted to blank spaces (normally, all underscores are automatically converted to blanks) */
				//	The following two flags are NOT LEGAL arguments to SetLabileFlags
				//
				newlineIsNotToken 		= 0x0400,
				hyphenIsPunctuation 	= 0x0800
				};
			const std::string & GetNextToken();
			void				SetSpecialPunctuationCharacter(char c);
			void				SetLabileFlagBit(int bit);
			std::string			GetErrorMessage() const {return errormsg;}
		private: 
			inline void		OLDAppendToComment(char);
			void			OLDReadComment();
			void			OLDReadParentheticalToken();
			void			OLDReadCurlyBracketedToken();
			void			OLDReadDoubleQuotedToken();
			void			OLDReadQuoted();
			bool			OLDIsPunctuation(char ch) const;
			std::string		comment;
			char			special;			/* ad hoc punctuation character; default value is '\0' */
			int				labileFlags;		/* storage for flags in the NxsTokenFlags enum */
			char			saved;
			std::string		errormsg;
#	else //defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
		enum NxsTokenFlags	/* arguments to AlterTokenReading used to change the standard tokenizing rules */
				{
				kNewlineIsToken        		= 0x0001, /* causes newline characters to be treated like punctuation (useful in reading data matrices) */
				kNewlineIsNotToken     		= 0x0002, /* causes newline to be treated as whitespace (this is the default state).  Used to undo  kNewlineIsToken */
				kHyphenNotPunctuation  		= 0x0004, /* causes hyphen to not be punctuation character (useful in reading in numbers) */
				kHyphenIsPunctuation   		= 0x0008  /* causes hyphen to  be punctuation character (this is the default state).  Used to undo  kHyphenNotPunctuation*/
				};
			const std::string & GetNextToken()
				{
				return ReadToken();
				}
#	endif //defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
		
			virtual ~NxsToken();
		//	Primary Services
		//
				const std::string & ReadToken();
				void				ReadSingleCharacter();
		  		bool				ReadCommandCommentOrToken();
				void				ReadAllTokensUntil(char tokenBreaker); // used to read the rest of double-quoted, parenthetical, or curly-bracketed token groups
				bool				ReadDoubleToken(double *dbl);
				NxsToken		   &operator++(); //equivalent to ReadToken();
		
		//	Accessors
		//
		bool      			AtEOF() const;
		bool       			AtEOL() const;
		bool       			Equals(const std::string &s, bool must_respect_case = false) const;
		NxsComment			GetComment(unsigned commentIndex) const;
		unsigned   			GetFileColumn() const;
		file_pos   			GetFilePosition() const;
		unsigned    	   	GetFileLine() const;
		unsigned			GetNumComments() const;
		unsigned			GetNextCommentIndex(char commandSpec = '\0',unsigned searchFrom = 0) const;
		std::string  		GetToken(bool must_respect_case = true) const;
		unsigned   			GetTokenLength() const;
		NxsTokenizerState	GetTokenizerState() const;
		const std::string & GetTokenReference() const;
		bool       			IsAbbreviation(const std::string &s) const;
		bool       			IsPunctuationToken() const;
		bool 				IsWhitespace(char ch) const;
		char				PeekAtNextChar() const;
		
		//	Modifiers
		//
		void       			AlterTokenReading(NxsTokenFlags bit);
		void       			ReplaceToken(const std::string &s);
		void       			ResetToken();
		void				SeekTokenizerState(const NxsTokenizerState &tpi);
		bool				SetEOFAllowed(bool b); //returns previous setting!
		void				ThrowIfNot(const std::string &s, bool respect_case = false) const;

		void				ThrowUnexpectedTokenNxsException(const std::string &) const;
	protected :
		//	
		//
		 
		void 		OutputComment(const std::string &msg) const;
	
	private :
		
		void 		AdvanceToNextCharInStream();
		void		AppendRestOfToken(char ch);
		void 		AppendToToken(char ch);
		bool		CharIsATokenByItself(char ch) const;
#		if defined(NXS_INTERNAL_TOKEN_BUFFER)
			void		FillInternalBuffer();
#		endif
		void		Initialize();
		bool		IsTokenBreaker(char) const;
		unsigned	ReadComment();
		char 		ReadFirstGraphChar();
		char 		ReadFirstGraphCharOrCommandComment();
		char 		ReadNextChar();
		void 		ReadQuoted();		
		std::string	ReadCommandCommentProgramIdentifier();
		class InternalStringAndStream
			{
			public:
				std::string 	nString;
				
#				if defined(NCL_USING_SSTREAM)
					std::stringstream 	strStream;
					
					InternalStringAndStream(const std::string &strCopy) 
						:nString(strCopy),
						strStream(nString.c_str())
						{}
#				else
					std::istrstream 	strStream;
					InternalStringAndStream(const std::string &strCopy) 
						:nString(strCopy),
						strStream(nString.c_str())
						{}
#				endif
			};
		
		class InternalPosInfo
			{
			public:
#				if defined(NXS_INTERNAL_TOKEN_BUFFER)
					InternalPosInfo();
					mutable file_pos 	posBeforeLastRead;	/* file position. the number of characters from the beginning of the stream */
					unsigned			posInLocalBuffer;
					unsigned			nCharsInLocalBuffer;
					bool				lastReadHitEnd;
#				else
					mutable file_pos 	filePos;	/* file position. the number of characters from the beginning of the stream */
#				endif
				unsigned       		fileLine;	/* the file line (number of newline characters since the beginning of the stream) */
				unsigned       		fileColumn;	/* the file column (number of characters since the last newline) */
				bool 				atEOF;		/* true if at the end of the file */
				bool				atEOL;		/* true if at the end of a line */
			};
		InternalStringAndStream *intStrStream;
#		if defined(NXS_INTERNAL_TOKEN_BUFFER)
			char * localBuffer;					/* the input stream that is being tokenized */
#		endif
		std::istream	   &inputStream;					/* the input stream that is being tokenized */
		char				nextCharInStream; 	/* the character that broke the current token */
		std::string			token;				/* the current token */
		InternalPosInfo     posInfo;			/* current position in stream info. Note:  posInfo refers the current character in the stream not the nextCharInStream */
		std::vector<NxsComment>	comments;			/* full information on all of the comments associated with the current token */
		bool				eofAllowed;			/* checked when eof is encountered.  if false, an NxsX_UnexpectedEOF exception is thrown.  Can be changed with SetEOFAllowed(bool) */
		char				currentSingleTokenStr[20];	/* list of characters that are single token characters (current punctuation) altered as necessary by AlterTokenReading() */
		char				currentTokenBreakerStr[23];	/* list of characters that are token breakers (current punctuation and whitespace) altered as necessary by AlterTokenReading() */
		char 				newLineCheck;				/* used to avoid extra slim down the IsWhitespace function.  If kNewlineIsToken, this char is set to ' ', so IsWhitespace('\n') won't return true */
		bool 				tokenStateHyphenIsPunc; /* Do NOT change this - use it to read the state ONLY (Use AlterTokenReading to change behaviour of the tokenizer) reports part of tokenizing state. */  

		friend class NxsTokenizerState;
	};

#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsToken NexusToken;
#endif
/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if ch is any nexus punctuation mark except [  ' 
|	[ and ' never get returned as punctuation, 
|	if kHyphenNotPunctuation the Caller must check to make sure the char isn't a -
|
*/
inline bool IsAlwaysSingleTokenPunctChar(char ch)
	{
	return (std::strchr("(){}\"-]/\\,;:=*`+<>", ch) != NULL);	
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	calls ReadToken().  Only prefix ++ is defined.
*/
inline NxsToken &NxsToken::operator++()
	{
	ReadToken(); 
	return *this;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the number of comments read between the end of the last token and the end of the current token.
|	Comments don't break tokens so a single token could contain several comments.  See GetNextCommentIndex()
*/
inline unsigned NxsToken::GetNumComments() const
	{
	return (unsigned)comments.size();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the current token has all of the letters in the capitalized prefix of s and doesn't disagree with other
|	characters in s.  
|	Written to avoid proliferation of ugly code like : IsCapAbbreviation(nxsToken.GetTokenReference(), s)
*/
inline bool NxsToken::IsAbbreviation(const std::string &s ) const
	{
	return IsCapAbbreviation(token, s);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the input stream is at its end.
*/
inline bool NxsToken::AtEOF() const
	{
	return posInfo.atEOF;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the current character (the last character in the current token) is an end of line
*/
inline bool NxsToken::AtEOL() const
	{
	return posInfo.atEOL;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	returns true if the current token equals s.
|	Written to avoid proliferation of ugly code like : nxsToken.GetTokenReference().Equal(s, ncl::kStringNoRespectCase)
*/
inline bool NxsToken::Equals(
  const std::string &s, 
  bool must_respect_case ) const	/* true to make the comparison case sensitive */
	{
	return StrEquals(token, s, (must_respect_case ? ncl::kStringRespectCase : ncl::kStringNoRespectCase) );
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns the column number (number of characters since last newline)
*/
inline unsigned  NxsToken::GetFileColumn() const
	{
	return posInfo.fileColumn;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns number of characters since the beginning of the file
*/
inline file_pos  NxsToken::GetFilePosition() const
	{
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		return file_pos(posInfo.posBeforeLastRead + file_pos(posInfo.posInLocalBuffer));
#	else
		posInfo.filePos = inputStream.rdbuf()->pubseekoff(0,std::ios::cur, std::ios::in);
		return posInfo.filePos;
#	endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	return the line number (number of newlines since the beginning of the file
*/
inline unsigned  NxsToken::GetFileLine() const
	{
	return posInfo.fileLine;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| 	return the state of the NxsToken object.  Used if the stream might need to be backed up to this point.  In which case
|	SeekTokenizerState is called with the stored NxsTokenizerState object as an argument.
*/
inline NxsTokenizerState NxsToken::GetTokenizerState() const
	{
	return NxsTokenizerState(*this);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns the next character in the stream. This is very fast because that character is always cached in the process of 
|	reading tokens so this doesn't involve a call to peek().
*/
inline char NxsToken::PeekAtNextChar() const
	{
	return nextCharInStream;
	}
			

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the current token is length one and is a punctuation char (){}\"-]/\\,;:=*`+<>
*/
inline bool NxsToken::IsPunctuationToken() const
	{
	return (GetTokenLength() == 1 && IsAlwaysSingleTokenPunctChar( token.at(0) ) );
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Clears the current token and all stored comments
*/
inline void NxsToken::ResetToken()
	{
	token.clear();
	comments.clear();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Used to control whether or not an exception is thrown when eof is encountered.  To throw exceptions (e.g. when entering 
|	a block) call SetEOFAllowed(false). then when returning to a portion of the file between blocks call SetEOFAllowed(true)
|	eof are allowed by default.
*/
inline bool NxsToken::SetEOFAllowed(
  bool b)
	{
	bool prev = eofAllowed;
	eofAllowed = b;
	return prev;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns true if ch is a command comment symbol or output comment symbol currently only recognizes ! \ and &
*/
inline bool IsCommentModifier(char ch)	//@need to add the other command comment ids
	{
	return (std::strchr("!\\&", ch) != NULL);	
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	returns true if ch is a single character token (punctuation other than single quote) (){}\"]/\\,;:=*`+<> always return 
|	true.  - might return true (it does by default, but AlterTokenReading allows this to be changed).
|	newline might return true (if AlterTokenReading(kNewlineIsToken) has been called)
*/
inline bool NxsToken::CharIsATokenByItself(char ch) const
	{
	return 	(std::strchr(currentSingleTokenStr, ch) != NULL);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if ch is any nexus punctuation mark or whitespace.
|	Whether or not hyphens are punctuation can be changed by AlterTokenReading()
*/
inline bool NxsToken::IsTokenBreaker(char ch) const
	{
	return 	(std::strchr(currentTokenBreakerStr, ch) != NULL);
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Reads the next character and then stops. To be used when character by character input is easier to handle.
|	NOTE 	 dangerous skips whitespace and comments, but this does not follow all nexus token rules. 
|			A single quote can be returned as a token.   single quotes will generate errors in the calling function!
*/
inline void NxsToken::ReadSingleCharacter()
	{
	token = ReadFirstGraphChar();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the character is a space, tab, newline or < 7.
|	Whether or not newlines are whitespace can be changed by AlterTokenReading()
*/
inline bool NxsToken::IsWhitespace(char ch) const 
	{
	return (ch == ' ' || (ch == newLineCheck) || ch == '\t' || ch < 7);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	
*/
inline NxsTokenPosInfo::NxsTokenPosInfo(
  file_pos p, 	/* current file position */
  unsigned l, 		/* current line number */
  unsigned c,		/* current column number */
  char nextC,
  bool	ateof,	/* true if at the EOF */
  bool  ateol)	/* true if at the end of a line */
	: filepos(p),
	fileline(l),
	filecol(c),
	atEOF(ateof),
	atEOL(ateol),
	nextCharInStream(nextC)
	{
	};

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns the position in file info (NxsTokenPosInfo object)
*/
inline NxsTokenPosInfo NxsTokenizerState::GetPosInfo() const 
	{
	return posInfo;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns the column number (number of characters since last newline)
*/
inline unsigned  NxsTokenPosInfo::GetFileColumn() const
	{
	return filecol;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns number of characters since the beginning of the file
*/
inline file_pos  NxsTokenPosInfo::GetFilePosition() const
	{
	return filepos;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	return the line number (number of newlines since the beginning of the file
*/
inline unsigned  NxsTokenPosInfo::GetFileLine() const
	{
	return fileline;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	returns the comment's location (number of characters from the beginning of the current token)
*/
inline unsigned NxsComment::GetLocationInToken() const 
	{
	return placeInToken;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	returns the comment, whatever was between [] (excluding the first character if it was a comment modifier).
*/
inline const std::string &NxsComment::GetCommentText() const 
	{
	return body;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns the comment, whatever was between [] (excluding the first character if it was a comment modifier).
*/
inline NxsComment::NxsComment(
  const std::string &s, 
  unsigned pos, 
  char comChar, 
  const std::string &progSpecifier)
  	: body(s),	/* the text of the comment */
  	placeInToken(pos),	/* the number of characters in the current token before the beginning of the comment was reached */
	modifierChar(comChar),	/* command/output comment specifier that applies to this commend, send '\0' if not applicable */
	intendedProgramCode(progSpecifier) /* program code (for use with command comments that are directed to only one program (used with [&&ProgramName ...] comments */
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns true if there is a modifier character and it isn't ! (the printing modifier).
*/
inline bool NxsComment::IsCommandComment() const
 	{
 	return ((modifierChar != '\\') && (modifierChar != '&'));
 	}


/*--------------------------------------------------------------------------------------------------------------------------
| 	returns true if comments modifier character (e.g. ! for output comment or & for some command comments) matches the 
|	argument.  
*/
inline bool NxsComment::MatchesModifier(
  char modChar) const
  	{
  	return (modifierChar == modChar);
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Function to improve code readability.
|	Throws an exception if the current token is not equal to the first argument.
*/
inline void NxsToken::ThrowIfNot(
  const std::string &s, 		/* string that should be identical to the current token. */
  bool respect_case ) const	/*true if the compare should be case sensitive */
	{
	if (!Equals(s, respect_case))
		ThrowUnexpectedTokenNxsException(s);
	}

#if defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
		/*----------------------------------------------------------------------------------------------------------------------
		|	Sets the bit specified in the variable `labileFlags'. The available bits are specified in the NxsTokenFlags enum.
		|	All bits in `labileFlags' are cleared after each token is read.
		*/
		inline void NxsToken::SetLabileFlagBit(
		  int bit)	/* the bit (see NxsTokenFlags enum) to set in `labileFlags' */
			{
			labileFlags |= bit;
			if (bit == newlineIsToken || bit == hyphenNotPunctuation)
				AlterTokenReading((NxsTokenFlags) bit);
			}
		/*----------------------------------------------------------------------------------------------------------------------
		|	Returns true if character supplied is considered a punctuation character. The following twenty characters are 
		|	considered punctuation characters:
		|>
		|	()[]{}/\,;:=*'"`+-<>
		|>
		|	Exceptions:
		|~
		|	o The tilde character ('~') is also considered punctuation if the tildeIsPunctuation labile flag is set
		|	o The special punctuation character (specified using the SetSpecialPunctuationCharacter) is also considered 
		|	  punctuation if the useSpecialPunctuation labile flag is set
		|	o The hyphen (i.e., minus sign) character ('-') is not considered punctuation if the hyphenNotPunctuation 
		|	  labile flag is set
		|~
		|	Use the SetLabileFlagBit method to set one or more NxsLabileFlags flags in `labileFlags'
		*/
		inline bool NxsToken::OLDIsPunctuation(
		  char ch) const	/* the character in question */
			{
			// PAUP 4.0b10 
			//  o allows ]`<> inside taxon names
			//  o allows `<> inside taxset names
			//
			bool is_punctuation = false;
			if (IsTokenBreaker(ch))
				return true;
			if (labileFlags & tildeIsPunctuation  && ch == '~')
				is_punctuation = true;
			if (labileFlags & useSpecialPunctuation  && ch == special)
				is_punctuation = true;
			if (labileFlags & hyphenNotPunctuation  && ch == '-')
				is_punctuation = false;

			return is_punctuation;
			}
		/*----------------------------------------------------------------------------------------------------------------------
		|	Adds `ch' to end of comment std::string.
		*/
		inline void NxsToken::OLDAppendToComment(
		  char ch)	/* character to be appended to comment */
			{
			comment += ch;  @@@-POL-@@@ convert to << operator
			}
#	endif // defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)

/*--------------------------------------------------------------------------------------------------------------------------
|	Appends the character ch to the current token. 
*/
inline void NxsToken::AppendToToken( char ch )
	{
   token << ch;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns a copy ofthe token as a std::string.
*/
inline std::string NxsToken::GetToken(
  bool must_respect_case) const
	{
	return (!must_respect_case ? GetCapitalized(token) : token);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns the size of the current token.
*/
inline unsigned NxsToken::GetTokenLength() const
	{
	return (unsigned) token.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the token for functions that only need read only access - faster than GetToken.
*/
inline const std::string &NxsToken::GetTokenReference() const
	{
	return token;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current token std::string with s.
*/
inline void NxsToken::ReplaceToken(
  const std::string &s)	/* std::string to replace current token std::string */
	{
	token = s;
	}

#if OLD_NXS_TOKEN || defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
	/*----------------------------------------------------------------------------------------------------------------------
	|	Sets the special punctuation character to `c'. If the labile bit useSpecialPunctuation is set, this character will 
	|	be added to the standard list of punctuation symbols, and will be returned as a separate token like the other 
	|	punctuation characters.
	*/
	inline void NxsToken::SetSpecialPunctuationCharacter(
	  char c)	/* the character to which `special' is set */
		{
		special = c;
		}

#endif




#endif
