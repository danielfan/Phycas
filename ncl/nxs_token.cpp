#include "phycas/force_include.h"
#if (MWERKS_LIB_BUILD)
#	pragma export on
#endif
#include "ncl/nxs_defs.hpp"
#include "ncl/nxs_token.hpp"
#include "ncl/nxs_exception.hpp"
#include "ncl/output/nxs_output.hpp"
#if (MWERKS_LIB_BUILD)
#	pragma export off
#endif

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::strcpy;
	using std::isgraph;
#endif
using std::string;

const unsigned kNxsTokenInternalBufferSize = 8192;

/*--------------------------------------------------------------------------------------------------------------------------
| 
*/
NxsToken::~NxsToken()
	{
	delete intStrStream;
#	if NXS_INTERNAL_TOKEN_BUFFER
		delete [] localBuffer;
#	endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	
*/
NxsTokenizerState::NxsTokenizerState(
  const NxsToken &t)	 /* current token */
  	: posInfo(t.GetFilePosition(), t.GetFileLine(), t.GetFileColumn(), t.nextCharInStream, t.AtEOF(), t.AtEOL()), 
  	token(t.GetTokenReference()) 
  	{
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		posBeforeLastRead  = t.posInfo.posBeforeLastRead;	/* file position. the number of characters from the beginning of the stream */
		posInLocalBuffer   = t.posInfo.posInLocalBuffer;
		nCharsInLocalBuffer = t.posInfo.posInLocalBuffer;
		lastReadHitEnd = t.posInfo.lastReadHitEnd;
#	endif
  	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	
*/
NxsTokenizerState::NxsTokenizerState()	 /* current token */
  	: posInfo(0,0,0,'\0',false,false), 
  	token(string()) 
  	{
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		posBeforeLastRead  = 0;	/* file position. the number of characters from the beginning of the stream */
		posInLocalBuffer   = 0;
		nCharsInLocalBuffer = 0;
		lastReadHitEnd = false;
#	endif
  	}



VecString NxsToken::TokenizeString(const std::string &s)
	{
	const string nstr(s);
	NxsToken token(nstr);
	VecString v;
	if (!token.AtEOF())
		{
		token.ReadToken();
		for (;;)
			{
			v.push_back(token.GetTokenReference());
			if (token.AtEOF())
				break;
			++token;
			}
		}
	return v;
	}
	
#if defined(NXS_INTERNAL_TOKEN_BUFFER)
	NxsToken::InternalPosInfo::InternalPosInfo()
		:posBeforeLastRead(0),
		posInLocalBuffer(0),
		nCharsInLocalBuffer(0),
		lastReadHitEnd(false),
		fileLine(0),
		fileColumn(0),
		atEOF(false),
		atEOL(false)
		{
		}
#endif

/*--------------------------------------------------------------------------------------------------------------------------
| 	Writes `s' to the output comment stream if the NxsToken::nxsOutCommentStream pointer is not NULL
*/
void NxsToken::OutputComment(const string &s) const
	{
	NxsOutput::GetNullCheckingStream<NxsOutput::kOutputComment>() << s.c_str() << ncl::endl;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Reads a token. treating - as if it were NOT a token breaker.
|	If the token can be treated as a double, dbl will be set to the double value and true will be returned
|	otherwise false will be returned (and dbl will not be changed).
*/
bool NxsToken::ReadDoubleToken(double *dbl)
	{
	const bool hyphenWasPuncOnEntry(tokenStateHyphenIsPunc);
	if (hyphenWasPuncOnEntry)
		AlterTokenReading(kHyphenNotPunctuation);
	try	{
		ReadToken();
		}
	catch (...)
		{
		if (hyphenWasPuncOnEntry)
			AlterTokenReading(kHyphenIsPunctuation);
		throw;
		}
	if (hyphenWasPuncOnEntry)
		AlterTokenReading(kHyphenIsPunctuation);
	return IsADouble(GetTokenReference(), dbl);
	}


#if defined(NXS_INTERNAL_TOKEN_BUFFER)
	void NxsToken::FillInternalBuffer()
		{
		posInfo.posBeforeLastRead = inputStream.tellg();
		posInfo.posInLocalBuffer = 0;
		if (posInfo.lastReadHitEnd)
			posInfo.nCharsInLocalBuffer = 0;
		else
			{
			inputStream.read(localBuffer, kNxsTokenInternalBufferSize);
			posInfo.nCharsInLocalBuffer = inputStream.gcount();
			posInfo.lastReadHitEnd = (posInfo.nCharsInLocalBuffer < kNxsTokenInternalBufferSize);
			}
		}
#endif

void NxsToken::Initialize()
	{
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		if (!localBuffer)
			localBuffer = new char [kNxsTokenInternalBufferSize];
		this->FillInternalBuffer();
#	else
		posInfo.filePos = inputStream.tellg();
#	endif
	posInfo.fileLine = 1L;
	posInfo.fileColumn = 1L;
	posInfo.atEOF = posInfo.atEOL = false;
 
#	if defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
		saved = '\0';
		special = '\0';
		labileFlags = 0;		
#	endif

	//	by default hyphen is a single token character, but newline isn't 
	//
	strcpy(currentSingleTokenStr,"(){}\"]/\\,;:=*`+<>-");
	currentSingleTokenStr[19]= '\0';
	
	//	by default hyphen- is a token breaker, (newline always is)
	//
	strcpy(currentTokenBreakerStr," \t\n;\'()]{}/\\,:=*\"`+<>-");
	
	//	read the first character into  nextCharInStream
	//
	nextCharInStream = 'a';	//anything other than EOF will work
	AdvanceToNextCharInStream();
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Creates a NxsToken to read the stream i.
*/
NxsToken::NxsToken(
  std::istream & i) 
  : intStrStream(NULL),
# if defined(NXS_INTERNAL_TOKEN_BUFFER)
	localBuffer(NULL), 
	inputStream(i),
# else
    inputStream(i),
# endif
   eofAllowed(true),
  newLineCheck('\n'),
  tokenStateHyphenIsPunc(true)
	{
	Initialize();
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Creates a NxsToken to read the stream i.
*/
NxsToken::NxsToken(const string & s)
  : intStrStream(new InternalStringAndStream(s)),
# if defined(NXS_INTERNAL_TOKEN_BUFFER)
    localBuffer(NULL),
	inputStream(intStrStream->strStream) , 
# else
    inputStream(intStrStream->strStream) , 
# endif
  eofAllowed(true),
  newLineCheck('\n'),
  tokenStateHyphenIsPunc(true)
	{
	Initialize();
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Advance the stream and store it in nextCharInStream.  Deal with the 3 ways of specifying return charaters 
|		(nextCharInStream will be set to \n if any of the return styles are found)
*/	
inline void NxsToken::AdvanceToNextCharInStream()
	{
	if (nextCharInStream == EOF)
		return;
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		if (posInfo.posInLocalBuffer == posInfo.nCharsInLocalBuffer)
			{
			if (posInfo.lastReadHitEnd == true) 
				{
				nextCharInStream = EOF;
				return;
				}
			this->FillInternalBuffer();
			}
		nextCharInStream  = this->localBuffer[posInfo.posInLocalBuffer++];
		if (nextCharInStream == 13 || nextCharInStream == 10)
			{
			nextCharInStream = '\n';
			if(nextCharInStream == 13)
				{
				if (posInfo.posInLocalBuffer == posInfo.nCharsInLocalBuffer)
					{
					this->FillInternalBuffer();
					if (posInfo.nCharsInLocalBuffer == 0)
						return;
					}
				if (this->localBuffer[posInfo.posInLocalBuffer] == 10)	//peeks at the next char
					posInfo.posInLocalBuffer++;
				}
			}
#	else
		nextCharInStream  = (char) (inputStream.rdbuf())->sbumpc();
		if (nextCharInStream == 13 || nextCharInStream == 10)
			{
			if(nextCharInStream == 13)
				{
				if ((inputStream.rdbuf())->sgetc() == 10)	//peeks at the next char
					(inputStream.rdbuf())->sbumpc();
				}
			nextCharInStream = '\n';
			}
#	endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	returns the character that had been stored in nextCharInStream, but also calls AdvanceToNextCharInStream() so 
|	nextCharInStream is advanced.
|	Does all of the fileposition bookkeeping.
|	Throws an NxsX_UnexpectedEOF exception if eof is found but eofAllowed is false.
*/
inline char NxsToken::ReadNextChar()
	{
	//	filepos = in.tellg();
	//
	// 	Why this was changed:  calls to tellg seem slow and unnecessary - we're storing filepos in terms of the
	//	number of times we call sbumpc().
	//	if we go back to getting the filepos via in.tellg(), remember to call it 
	//	twice after both sgetc() calls in the case of the \13\10 endline
	
	char ch = nextCharInStream;
	AdvanceToNextCharInStream();
	if(ch == EOF)
		{
		posInfo.atEOF = true;
		if (eofAllowed)
			return '\0';
		throw NxsX_UnexpectedEOF(*this);
		}
	if(ch == '\n')
		{
		posInfo.fileLine++;
		posInfo.fileColumn = 1L;
		posInfo.atEOL = true;
		return '\n';
		}
	if (ch == '\t')
		posInfo.fileColumn += 4 - ((posInfo.fileColumn - 1)%4);	//@assumes that tab will be 4 in the editor we use
	else
		posInfo.fileColumn++;
	posInfo.atEOL = false;
	return ch;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Reads in a comment.  Called after [ is found.  a NxsComment is created and pushed back onto the vector of comments in 
|	this token.  The index of the comment in the vector is returned.
*/
string	NxsToken::ReadCommandCommentProgramIdentifier()
	{
	string s;
	char ch = ' ';
	while (!isgraph(ch) && nextCharInStream != ']')
		ch = ReadNextChar();
	while (isgraph(ch) && nextCharInStream != ']')
		{
		s << ch;
		ch = ReadNextChar();
		}
	return s;
	}
	
unsigned NxsToken::ReadComment()
	{
	bool prevAllowEOF = eofAllowed;
	SetEOFAllowed(false);
	unsigned retVal = (unsigned)comments.size();
	try	{
		string progSpec;
		string commentBody;
		char modChar = '\0';
		// see if first character is the output comment symbol ('!')
		// or command comment symbol (&)
		char ch = ReadNextChar();
		if(IsCommentModifier(ch))
			{
			modChar = ch;
			if (ch == '&' && nextCharInStream == '&')
				{
				ReadNextChar();
				progSpec = ReadCommandCommentProgramIdentifier();
				}
			ch = ReadNextChar();
			}
		for(int level = 1;; ch = ReadNextChar())	
			{
			if(ch == ']')
				{
				--level;
				if(level == 0)
					break;
				}
			else if(ch == '[')
				++level;
			commentBody << ch;
			}
		comments.push_back(NxsComment(commentBody, GetTokenLength(), modChar, progSpec));
	  	if(modChar == '!') 
			OutputComment(commentBody);
		}
	catch (NxsX_UnexpectedEOF & x)
		{
		x.msg << " in a comment";
		throw x;
		}
	SetEOFAllowed(prevAllowEOF);
	return retVal;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Gets remainder of a quoted Nexus word (the first single quote character was read in already by ReadToken). This 
|	function reads characters until the next single quote is encountered.  An exception occurs if two single quotes occur 
|	one after the other, in which case the function continues to gather characters until an isolated single quote is found. 
|	The tandem quotes are stored as a single quote character in the token string.
*/
void NxsToken::ReadQuoted()
	{		// Note: within quoted tokens, underscores  should be preserved
	// always throw a NxsX_UnexpectedEOF if you reach the eof in the middle of a token
	bool prevEOFAllowed = eofAllowed;
	eofAllowed = false;
	for(;;)
		{
		char ch = ReadNextChar();
		if (ch == '\'')
			{
			if (nextCharInStream== '\'') 
				ReadNextChar();	//skip the second '
			else
				{
				eofAllowed = prevEOFAllowed;
				return;
				}
			}
		AppendToToken(ch);
		}
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Called to read the rest of a token (everything after the character ch).  This function assumes that the beginning of 
|	a token has been found (through GetNextGraphChar())
*/
void NxsToken::AppendRestOfToken(
  char ch)	/* the last character read before this function was called*/
	{	
	if (CharIsATokenByItself(ch))
		AppendToToken(ch);
	else
		{
		if(ch == '\'')
			{
			ReadQuoted();
#			if defined (NXS_THROW_IF_EMPTY_TOKENS)
				if (token.empty())
					throw NxsX_EmptyToken(*this);
#			else
				return;
#			endif
			}
		else
			{
			AppendToToken(ch);
			for(;!IsTokenBreaker(nextCharInStream) && !AtEOF();)
				{
				if (nextCharInStream == '[')
					{
					ReadNextChar();
					ReadComment();
					}
				else if (nextCharInStream == '_')
					{
					AppendToToken(' ');
					ReadNextChar();
					}
				else 
					{
					ch = ReadNextChar();
					if (ch != '\0')
						AppendToToken(ch);
					else
						return;
					}
				}
			}
		}
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Reads characters from in until a complete token has been read and stored in token.  Performs a number of useful 
|	operations in the process of retrieving tokens:
|		underscore characters encountered are stored as blank spaces (except in single quoted tokens).
|		strings inside single quotes are returned without the surrounding quotes.  Paired single quotes are converted to '
|		comments are handled automatically (normal comments are stored and output comments are passed to the function 
|			the virtual OutputComment(string).
|		leading whitespace is automatically skipped
|		if the end of the file is reached on reading this token, the atEOF flag is set and may be queried using the AtEOF()
|		punctuation characters are always returned as individual tokens.
*/
const std::string & NxsToken::ReadToken()
	{
	ResetToken();
	char ch = ReadFirstGraphChar();
	if (ch != '\0')
		AppendRestOfToken(ch);
	return GetTokenReference();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns True if command comment was the only thing found (in which case the token length will still be zero
|	it is possible to return false but still have command comments (if they are embedded in a token)
*/
bool NxsToken::ReadCommandCommentOrToken()
	{
	ResetToken();
	
	char ch = ReadFirstGraphCharOrCommandComment();
	//	ch will be 0 if we read a command comment only if this is the case, we're done
	//
	if (ch != '\0')
		{
		AppendRestOfToken(ch);
		return false;
		}
	return !AtEOF();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This function is very limited in usefulness.
|	Squashes together (with no whitespace, unless the newline is token is on) all of the tokens until the character 
|	tokenBreaker is encountered.  
|	NOTE: tokenBreaker MUST already be a punctuation character for this function to work.
|	NOTE: NONE of the tokens read can be single-quoted tokens
|	ALSO NOTE:  the terminating tokenBreaker is NOT added to the token.
|	this function is intended for reading the contents of double-quoted strings, () and {}
*/
void NxsToken::ReadAllTokensUntil(
  char tokenBreaker)
  	{
  	assert(IsTokenBreaker(tokenBreaker));
  	ResetToken();
	char ch = ReadFirstGraphChar();
	while (ch != tokenBreaker)
		{
		if (AtEOF())
			{
			string errormsg = "Unexpected end of file while searching for a closing ";
			errormsg << tokenBreaker << " character";
			throw NxsException(errormsg, *this);
			}
		if (ch == '\'')
			{
			string errormsg = "Unexpected single-quoted word while searching for a closing";
			errormsg << tokenBreaker << " character";
			throw NxsException(errormsg, *this);
			}
		if (ch != '\0')	// ch could be \0 if a "bare" command comment is found
			AppendRestOfToken(ch);
		ch = ReadFirstGraphChar();	
		} 
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Returns next darkspace or 0 for the eof.
|	If first darkspace is [ then the comment is read, and the search for the first graphical character continues
*/
char NxsToken::ReadFirstGraphChar()
	{
	if (AtEOF())
		return '\0';
	char ch = ' ';
	do	{
		if (nextCharInStream == '[')
			{
			ReadNextChar();
			ReadComment();
			}
		else 
			{
			ch = ReadNextChar();
			if (AtEOF())
				return '\0';
			}
		}
	while (IsWhitespace(ch));
	if (ch == '_')
		return ' ';
	return ch;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Returns next darkspace or 0 (if the eof or a command comment is found)
|	Just like ReadFirstGraphChar, except that it returns '\0' if a command comment is found.
*/
char NxsToken::ReadFirstGraphCharOrCommandComment()
	{
	if (AtEOF())
		return '\0';
	char ch = ' ';
	do	{
		if (nextCharInStream == '[')
			{
			ReadNextChar();
			unsigned comInd = ReadComment();
			if (comments[comInd].IsCommandComment())
				return '\0';
			if (IsWhitespace(nextCharInStream))
				{
				comments.clear();//flush comments that aren't adjacent to a token;
				ReadNextChar();
				}
			else
				ch = ReadNextChar();
			}
		else 
			{
			ch = ReadNextChar();
			if (AtEOF())
				return '\0';
			}
		}
	while (IsWhitespace(ch));
	if (ch == '_')
		return ' ';
	return ch;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Sets the NxsToken to the file position and the token described in the argument.  The tokenizer state should have been
|	obtained via a call of GetTokenizerState() with this NxsToken object.
|	Note reported file positions in the NxsToken are actually one character off (hence AdvanceToNextCharInStream() is 
|	called after the file info is set).
*/
void NxsToken::SeekTokenizerState(
  const NxsTokenizerState &tpi)
	{
	using std::ios;
	token = tpi.token;
	posInfo.fileColumn = tpi.posInfo.filecol;
	posInfo.fileLine = tpi.posInfo.fileline;
	posInfo.atEOF = tpi.posInfo.atEOF;
	posInfo.atEOL = tpi.posInfo.atEOL;
#	if defined(NXS_INTERNAL_TOKEN_BUFFER)
		if (tpi.posBeforeLastRead != posInfo.posBeforeLastRead)
			{
				//need to make the local buffer's match.
			inputStream.seekg(posInfo.posBeforeLastRead);
			FillInternalBuffer();
			if (tpi.posInLocalBuffer > posInfo.nCharsInLocalBuffer)
				throw NxsX_InvalidTokenizerState(); //this shouldn't happen.  Should try to figure out if these cases are IO errors, or bad arguments to this function
			}
		posInfo.lastReadHitEnd = tpi.lastReadHitEnd;
		posInfo.posInLocalBuffer = tpi.posInLocalBuffer;
		nextCharInStream = tpi.posInfo.nextCharInStream;
#	else
		posInfo.filePos = tpi.posInfo.filepos;
		inputStream.seekg(posInfo.filePos);
		inputStream.rdbuf()->pubseekoff(-1, ios::cur, ios::in);
		nextCharInStream = tpi.posInfo.nextCharInStream;
		AdvanceToNextCharInStream();
#	endif

	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	returns the NxsComment with commentIndex.  Throws XInvalidCommentIndex if the commentIndex is greater than the number of
|	comments stored.
*/
NxsComment NxsToken::GetComment(
  unsigned commentIndex) const
	{
	if (commentIndex >= comments.size())
		throw NxsX_InvalidCommentIndex();
	return comments[commentIndex];
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the index of the first command comment that has the same commandSpecifier (commandSpec) and with an 
|	index >= searchFrom.
|	If commandSpec is '\0' ALL comments are returned
*/
unsigned NxsToken::GetNextCommentIndex(
  char commandSpec,
  unsigned searchFrom) const
  	{
  	for (; searchFrom < comments.size(); ++searchFrom)
  		{
  		if (commandSpec == '\0' || comments[searchFrom].MatchesModifier(commandSpec))
  			return searchFrom;
  		}
  	return UINT_MAX;
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Used to modify the standard tokenizing rules to make it easier to read some parts of the Nexus format.
|	Should be called with the argument kNewlineIsToken if newlines need to be read as tokens (for instance when reading 
|		interleaved matrices).  This setting will continue to affect the NxsToken until AlterTokenReading(kNewlineIsNotToken)
|		is called.
|	Use AlterTokenReading(kHyphenNotPunctuation) to read in tokens that might contain hyphens as negative signs (to stop - from
|		breaking tokens).  AlterTokenReading(kHyphenIsPunctuation) can then be called to return to the default state of 
|		treating hyphens as punctuation characters.
|
|	Implementation note:  	this code is pretty odd looking.  the character arrays currentSingleTokenStr and currentTokenBreakerStr
|		are used as argument to strchr when NxsToken is deciding if the next character constitutes a token by itself or should
|		break the current token.  Calls to AlterTokenReading modify these character arrays. 
|		IsWhitespace checks if the character is ==' ' or =='\t' or < 7 or == newLineCheck.  So when AlterTokenReading(kNewlineIsToken)
|		is called, newLineCheck is set to ' '.  This means that \n will no longer trigger yes to the question IsWhitespace.
|		AlterTokenReading(kNewlineIsNotToken) returns newLineCheck to the more intuitive value of '\n'
|
 */
 void NxsToken::AlterTokenReading(
   NxsTokenFlags bit ) /* specifies what tokenizing change is being requested */
	{
	switch (bit)
		{
   		case (kNewlineIsToken) :
    		if (currentSingleTokenStr[17] =='-')
    			currentSingleTokenStr[18] = '\n';
    		else
    			currentSingleTokenStr[17] = '\n';
    		//currentWhitespaceStr[3] = '\0';
    		newLineCheck = ' ';
    		break;
   		
    	case (kNewlineIsNotToken) :
    		if (currentSingleTokenStr[17] == '\n')
    			currentSingleTokenStr[17] = currentSingleTokenStr[18];
    		currentSingleTokenStr[18] = '\0';
    		//currentWhitespaceStr[3] = '\n';
    		newLineCheck = '\n';
    		break;
   		
    	case (kHyphenNotPunctuation) :
    		if (currentSingleTokenStr[17] == '-')
    			currentSingleTokenStr[17] = currentSingleTokenStr[18];
    		currentSingleTokenStr[18] = '\0';
    		currentTokenBreakerStr[21] = '\0';
    		tokenStateHyphenIsPunc = false;
    		break;
    		
   		case (kHyphenIsPunctuation) :
    		if (currentSingleTokenStr[17] =='\n')
    			currentSingleTokenStr[18] = '-';
    		else
    			currentSingleTokenStr[17] = '-';
    		currentTokenBreakerStr[21] = '-';
    		tokenStateHyphenIsPunc = true;
    		break;
   		
    	}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Throws a NxsException with the message "Expecing s but found current_token" where s is the argument, and current_token
|	is the string currently in the NxsToken.
|	Written so that ThrowIfNot() could be inlined with minimal code-bloat (we don't need to inline this function because
|	exception throwing is inherently slow).
*/
void NxsToken::ThrowUnexpectedTokenNxsException(
  const string &s) const	/* string that should be identical to the current token. */
	{
	string e;
	e << s << " was exepected, but " << token << " was entered.";
	throw NxsException(e, *this);
	}


#	if defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)
		/*----------------------------------------------------------------------------------------------------------------------
		|	Reads characters from in until a complete token has been read and stored in token. GetNextToken performs a number 
		|	of useful operations in the process of retrieving tokens:
		|~
		|	o any underscore characters encountered are stored as blank spaces (unless the labile flag bit preserveUnderscores
		|	  is set)
		|	o if the first character of the next token is an isolated single quote, then the entire quoted string is saved 
		|	  as the next token
		|	o paired single quotes are automatically converted to single quotes before being stored
		|	o comments are handled automatically (normal comments are treated as whitespace and output comments are passed to 
		|	  the function OutputComment which does nothing in the NxsToken class but can be overridden in a derived class to 
		|	  handle these in an appropriate fashion)
		|	o leading whitespace (including comments) is automatically skipped
		|	o if the end of the file is reached on reading this token, the atEOF flag is set and may be queried using the AtEOF 
		|	  member function
		|	o punctuation characters are always returned as individual tokens (see the Maddison, Swofford, and Maddison paper 
		|	  for the definition of punctuation characters) unless the flag ignorePunctuation is set in labileFlags,
		|	  in which case the normal punctuation symbols are treated just like any other darkspace character.
		|~
		|	The behavior of GetNextToken may be altered by using labile flags. For example, the labile flag saveCommandComments 
		|	can be set using the member function SetLabileFlagBit. This will cause comments of the form [&X] to be saved as 
		|	tokens (without the square brackets), but only for the aquisition of the next token. Labile flags are cleared after 
		|	each application.
		*/
		const std::string & NxsToken::GetNextToken()
			{
			ResetToken();

			char ch = ' ';
			if (saved == '\0' || IsWhitespace(saved))
				{
				// Skip leading whitespace
				//
				while( IsWhitespace(ch) && !AtEOF())
					ch = ReadNextChar();
				saved = ch;
				}

			for(;;)
				{
				// Break now if singleCharacterToken mode on and token length > 0.
				//
				if (labileFlags & singleCharacterToken && token.size() > 0)
					break;

				// Get next character either from saved or from input stream.
				//
				if (saved != '\0')
					{
					ch = saved;
					saved = '\0';
					}
				else
					ch = ReadNextChar();

				// Break now if we've hit EOF.
				//
				if (AtEOF())
					break;

				if (ch == '\n' && labileFlags & kNewlineIsToken)
					{
					if (token.size() > 0)
						{
						// Newline came after token, save newline until next time when it will be 
						// reported as a separate token.
						//
						posInfo.atEOL = 0;
						saved = ch;
						}
					else
						{
						posInfo.atEOL = 1;
						AppendToToken(ch);
						}
					break;
					}

				else if (IsWhitespace(ch))
					{
					// Break only if we've begun adding to token (remember, if we hit a comment before a token,
					// there might be further white space between the comment and the next token).
					//
					if (token.size() > 0)
						break;
					}

				else if (ch == '_')
					{
					// If underscores are discovered in unquoted tokens, they should be 
					// automatically converted to spaces.
					//
					if (!(labileFlags & preserveUnderscores))
						ch = ' ';
					AppendToToken(ch);
					}

				else if (ch == '[')
					{
					// Get rest of comment and deal with it, but notice that we only break if the comment ends a token,
					// not if it starts one (comment counts as whitespace). In the case of command comments 
					// (if saveCommandComment) GetComment will add to the token string, causing us to break because
					// token.size() will be greater than 0.
					//
					OLDReadComment();
					if (token.size() > 0)
					break;
					}

				else if (ch == '(' && labileFlags & parentheticalToken)
					{
					AppendToToken(ch);

					// Get rest of parenthetical token.
					//
					OLDReadParentheticalToken();
					break;
					}

				else if (ch == '{' && labileFlags & curlyBracketedToken)
					{
					AppendToToken(ch);

					// Get rest of curly-bracketed token.
					//
					OLDReadCurlyBracketedToken();
					break;
					}

				else if (ch == '\"' && labileFlags & doubleQuotedToken)
					{
					// Get rest of double-quoted token.
					//
					OLDReadDoubleQuotedToken();
					break;
					}

				else if (ch == '\'')
					{
					if (token.size() > 0)
						{
						// We've encountered a single quote after a token has
						// already begun to be read; should be another tandem
						// single quote character immediately following.
						//
						ch = ReadNextChar();
						if (ch == '\'')
							AppendToToken(ch);
						else
							{
							errormsg = "Expecting second single quote character";
							throw NxsException( errormsg, GetFilePosition(), GetFileLine(), GetFileColumn());
							}
						}
					else
						{
						// Get rest of quoted NEXUS word and break, since
						// we will have eaten one token after calling GetQuoted.
						//
						OLDReadQuoted();
						}
					break;
					}

				else if (OLDIsPunctuation(ch))
					{
					if (token.size() > 0)
						{
						// If we've already begun reading the token, encountering
						// a punctuation character means we should stop, saving
						// the punctuation character for the next token.
						//
						saved = ch;
						break;
						}
					else
						{
						// If we haven't already begun reading the token, encountering
						// a punctuation character means we should stop and return
						// the punctuation character as this token (i.e., the token
						// is just the single punctuation character.
						//
						AppendToToken(ch);
						break;
						}
					}

				else
					{
					AppendToToken(ch);
					}

				}

			if (labileFlags & kNewlineIsToken)
			 	AlterTokenReading(kNewlineIsNotToken);
			if (labileFlags & kHyphenNotPunctuation)
				AlterTokenReading(kHyphenIsPunctuation);
			labileFlags = 0;
			return GetTokenReference();
			}
	/*----------------------------------------------------------------------------------------------------------------------
	|	Reads rest of comment (starting '[' already input) and acts accordingly. If comment is an output comment, and if 
	|	an output stream has been attached, writes the output comment to the output stream. Otherwise, output comments are 
	|	simply ignored like regular comments. If the labileFlag bit saveCommandComments is in effect, the comment (without 
	|	the square brackets) will be stored in token. 
	*/
	void NxsToken::OLDReadComment()
		{
		// Set comment level to 1 initially.  Every ']' encountered reduces
		// level by one, so that we know we can stop when level becomes 0.
		//
		int level = 1;

		// Get first character
		//
		char ch = ReadNextChar();
		if (AtEOF())
			{
			errormsg = "Unexpected end of file inside comment";
			throw NxsException( errormsg, GetFilePosition(), GetFileLine(), GetFileColumn());
			}

		// See if first character is the output comment symbol ('!')
		// or command comment symbol (&)
		//
		int printing = 0;
		int command = 0;
		if (ch == '!')
			printing = 1;
		else if (ch == '&' && labileFlags & saveCommandComments)
			{
			command = 1;
			AppendToToken(ch);
			}
		else if (ch == ']')
			return;

		// Now read the rest of the comment
		//
		for(;;)
			{
			ch = ReadNextChar();
			if (AtEOF())
				break;

			if (ch == ']')
				--level;
			else if (ch == '[')
				++level;

			if (level == 0)
				break;

			if (printing)
				OLDAppendToComment(ch);
			else if (command)
				AppendToToken(ch);
			}

		if (printing)
			{
			// Allow output comment to be printed or displayed in most appropriate
			// manner for target operating system
			//
			OutputComment(comment);
			}
		
		// Now that we are done with it, free the memory used to store the comment
		//
		comment.clear();
		}


	/*----------------------------------------------------------------------------------------------------------------------
	|	Reads rest of a token surrounded with curly brackets (the starting '{' has already been input) up to and including
	|	the matching '}' character. All nested curly-bracketed phrases will be included.
	*/
	void NxsToken::OLDReadCurlyBracketedToken()
		{
		// Set level to 1 initially.  Every '}' encountered reduces
		// level by one, so that we know we can stop when level becomes 0.
		//
		int level = 1;

		char ch;
		for(;;)
			{
			ch = ReadNextChar();
			if (AtEOF())
				break;

			if (ch == '}')
				--level;
			else if (ch == '{')
				++level;

			AppendToToken(ch);

			if (level == 0)
				break;
			}
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Gets remainder of a double-quoted NEXUS word (the first double quote character was read in already by GetNextToken).
	|	This function reads characters until the next double quote is encountered. Tandem double quotes within a 
	|	double-quoted NEXUS word are not allowed and will be treated as the end of the first word and the beginning of the 
	|	next double-quoted NEXUS word. Tandem single quotes inside a double-quoted NEXUS word are saved as two separate 
	|	single quote characters; to embed a single quote inside a double-quoted NEXUS word, simply use the single quote by 
	|	itself (not paired with another tandem single quote).
	*/
	void NxsToken::OLDReadDoubleQuotedToken()
		{
		char ch;

		for(;;)
			{
			ch = ReadNextChar();
			if (AtEOF())
				break;

			if (ch == '\"')
				break;
			else
				AppendToToken(ch);
			}
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Gets remainder of a quoted NEXUS word (the first single quote character was read in already by GetNextToken). This
	|	function reads characters until the next single quote is encountered. An exception occurs if two single quotes occur
	|	one after the other, in which case the function continues to gather characters until an isolated single quote is
	|	found. The tandem quotes are stored as a single quote character in the token string.
	*/
	void NxsToken::OLDReadQuoted()
		{
		char ch;

		for(;;)
			{
			ch = ReadNextChar();
			if (AtEOF())
				break;

			if (ch == '\'' && saved == '\'')
				{
				// Paired single quotes, save as one single quote
				//
				AppendToToken(ch);
				saved = '\0';
				}
			else if (ch == '\'' && saved == '\0')
				{
				// Save the single quote to see if it is followed by another
				//
				saved = '\'';
				}
			else if (saved == '\'')
				{
				// Previously read character was single quote but this is something else, save current character so that it will
				// be the first character in the next token read
				//
				saved = ch;
				break;
				}
			else
				AppendToToken(ch);
			}
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Reads rest of parenthetical token (starting '(' already input) up to and including the matching ')' character.  All
	|	nested parenthetical phrases will be included.
	*/
	void NxsToken::OLDReadParentheticalToken()
		{
		// Set level to 1 initially.  Every ')' encountered reduces
		// level by one, so that we know we can stop when level becomes 0.
		//
		int level = 1;

		char ch;
		for(;;)
			{
			ch = ReadNextChar();
			if (AtEOF())
				break;

			if (ch == ')')
				--level;
			else if (ch == '(')
				++level;

			AppendToToken(ch);

			if (level == 0)
				break;
			}
		}

#	endif	//defined (SUPPORT_NXS_TOKEN_LABILE_FLAGS)



