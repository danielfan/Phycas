// This file is part of Oak, a program for likelihood and bayesian
#ifndef NCL_NXS_CMD_OPTION_PRIMITIVES_H
#define NCL_NXS_CMD_OPTION_PRIMITIVES_H
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_cmd_param.hpp"
#include <boost/function.hpp>

/*----------------------------------------------------------------------------------------------------------------------
|	reads boolean values.
|	Note because this class has to recognized NoName as well as Name,  CanReadKeyword sets the mutable recognizedAsConverse
|	flag to indicate which form of the command was encountered.  This means that the ReadValue function only behaves
|	correctly if CanReadKeyword has been called immediately before.
*/
class BoolCmdOption : public SimpleCmdOptionInterface<bool>
	{
	public :
		bool		CanReadKeyword(const std::string &s, int permLevel) const; 
		std::string GetConverseName() const;
		std::string	GetCurrentValueAsString() const;
		std::string GetDisplayType(bool includeIndefiniteArticle,  bool plural) const;
		VecString   GetValidArgument()
			{
			VecString s(2, "YES");
			s[1] = "NO";
			return s;
			}
		bool		IsCurrentlyValid() {return true;}
  		void 		StorePreviousValue();
		bool		ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		void		SetConverseAbbreviation(const std::string &);
		void		WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
		
		bool operator==(const std::string &) const;
		BoolCmdOption( 	const std::string &n, 
						bool *manipVal, 
						bool def, 
						bool persist, 
						CmdPermissionLevel pLevel );
	protected:	
		
		std::string 		converseAbbrev;
		mutable bool	recognizedAsConverse;
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets recognizedAsConverse to false and calls base class StorePreviousValue
*/
inline void BoolCmdOption::StorePreviousValue()
	{
	recognizedAsConverse = false;
	SimpleCmdOptionInterface<bool>::StorePreviousValue();
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Reads double quotes strings as a Single string.  WHITESPACE IS REMOVED !! 
*/
class QuotedStringCmdOption : public NxsStringCmdOption
	{
	public :
		bool	ReadValue( NxsToken &token, bool equalsAlreadyRead);
		QuotedStringCmdOption(const std::string &n, std::string *manipVal, std::string def, bool persist, CmdPermissionLevel pLevel );
	};

/*----------------------------------------------------------------------------------------------------------------------
|	returns "boolean"
*/
inline std::string BoolCmdOption::GetDisplayType(
  bool includeIndefiniteArticle,
  bool plural) const
	{
	return StandardNoun("boolean",includeIndefiniteArticle, plural);
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string BoolCmdOption::GetCurrentValueAsString() const
	{
	return (*(this->currentValue) ? "true" : "false");
	}
	
template <typename T>
class NamedNumberSource
	{
	public:
		typedef boost::function0<T> CallbackPtr;
		NamedNumberSource()
			:ptr(NULL)
			{}
		NamedNumberSource(const std::string &n, T *p)
			:ptr(p),
			name(n)
			{
			NXS_ASSERT(p);
			}
		NamedNumberSource(const std::string &n, CallbackPtr p)
			:callback(p),
			ptr(NULL),
			name(n)
			{
			NXS_ASSERT(p);
			}
		operator bool () const // implicit cast to bool
			{
			return (ptr != NULL) || callback;
			}
		const std::string & GetName() const
			{
			return name;
			}
		std::string GetDescription(const T & v) const
			{
			const std::string & n = GetName();
			if (n.empty())
				{
				std::string s;
				s << v;
				return s;
				}
			return n;
			}
		T operator()() const
			{
			return (ptr != NULL ? *ptr : (callback ? callback() : T(0)));
			}
	protected:
		CallbackPtr		callback;
		const T		  * ptr;	
		std::string		name;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	template class for the standard numerical types
*/
template <typename T> class NumberCmdOption: public SimpleCmdOptionInterface<T>
	{
	public :
		typedef NamedNumberSource<T> NumberSource;
						NumberCmdOption(const std::string &n, T *manipVal, T def, NumberSource minV, NumberSource maxV, bool persist, CmdPermissionLevel pLevel );
						NumberCmdOption(const std::string &n, T *manipVal, T def,			 T minV, NumberSource maxV, bool persist, CmdPermissionLevel pLevel );
						NumberCmdOption(const std::string &n, T *manipVal, T def, NumberSource minV,			T maxV, bool persist, CmdPermissionLevel pLevel );
						NumberCmdOption(const std::string &n, T *manipVal, T def,			 T minV,			T maxV, bool persist, CmdPermissionLevel pLevel );
		std::string 	GetDisplayType(bool includeIndefiniteArticle,  bool plural) const;
  		bool			IsCurrentlyValid();
		bool			ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		std::string		GetCurrentValueAsString() const;
		VecString		GetValidArgument();
		void			WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
		void			WriteTypeInfoStateElementPrivate(NxsHiddenQueryStream & outStream, const std::string & tagName) const;

	protected:	
		enum ValueDescriptionIdentifier
			{
			kMinValueDescription,  //kDefault... should be added when the default is allowed to be labile
			kMaxValueDescription
			};
		std::string DescribeValue(typename NumberCmdOption<T>::ValueDescriptionIdentifier d) const
			{
			std::string s;
			RefreshBounds();
			if (d == kMinValueDescription)
				{
				if (minimumProvider)
					return minimumProvider.GetDescription(minVal);
				s << minVal;
				}
			else if (d == kMaxValueDescription)
				{
				if (maximumProvider)
					return maximumProvider.GetDescription(maxVal);
				s << maxVal;
				}
			return s;
			}

		void 		RefreshBounds() const;
		T			ReadValueAfterEquals(NxsToken &);
		void		FlagTooSmall();
		void		FlagTooLarge();
		
		NumberSource maximumProvider;
		NumberSource minimumProvider;

		mutable T   maxVal;
		mutable T   minVal;
		
	};
template<typename T>
inline std::string NumberCmdOption<T>::GetCurrentValueAsString() const
	{
	std::string s;
	return AppendNumber<T>(s, *(this->currentValue));
	}

		
  		
/*----------------------------------------------------------------------------------------------------------------------
|	Flag's the errState as too_big or too_big_labile
*/
template<typename T> 
inline void NumberCmdOption<T>::FlagTooLarge()
	{
	std::string s;
	AppendNumber<T>(s, maxVal);
	if (maximumProvider)
		NxsCmdOption::FlagError(NxsCmdOption::too_big_labile, s);
	else
		NxsCmdOption::FlagError(NxsCmdOption::too_big, s);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Flag's the errState as too_small or too_small_labile
*/
template<typename T> inline void NumberCmdOption<T>::FlagTooSmall()
	{
	std::string s;
	AppendNumber<T>(s, minVal);
	if (minimumProvider)
		NxsCmdOption::FlagError(NxsCmdOption::too_small_labile);
	else
		NxsCmdOption::FlagError(NxsCmdOption::too_small, s);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Checks to make sure the currentValue is within the bounds
*/
template<typename T> inline bool NumberCmdOption<T>::IsCurrentlyValid()
	{
	RefreshBounds();
	if (*this->currentValue < minVal)
		{
		FlagTooSmall();
		return false;
		}
	if (*this->currentValue > maxVal)
		{
		FlagTooLarge();
		return false;
		}
	return true;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Checks to make sure the currentValue is within the bounds
*/
template<typename T> inline VecString NumberCmdOption<T>::GetValidArgument()
	{
	RefreshBounds();
	VecString s;
	std::string strForm;
	AppendNumber<T>(strForm, minVal);
	s.push_back(strForm);
	strForm.clear();
	AppendNumber<T>(strForm, maxVal);
	s.push_back(strForm);
	return s;
	}
	
template<typename T> inline void NumberCmdOption<T>::RefreshBounds() const
	{
	if (minimumProvider)
		minVal = minimumProvider();
	if (maximumProvider)
		maxVal = maximumProvider();
	}

//	template specializations for direct interactions with std::string
//
template <>
inline double NumberCmdOption<double>::ReadValueAfterEquals(
  NxsToken&token)
	{
	return ConvertToDoubleOrThrow(token.GetTokenReference());
	}
	
template <>
inline int NumberCmdOption<int>::ReadValueAfterEquals(
  NxsToken&token)
	{
	return ConvertToIntOrThrow(token.GetTokenReference());
	}

template <>
inline unsigned NumberCmdOption<unsigned>::ReadValueAfterEquals(
  NxsToken&token)
	{
	return ConvertToUnsignedOrThrow(token.GetTokenReference());
	}
	
template <>
inline long NumberCmdOption<long>::ReadValueAfterEquals(
  NxsToken&token)
	{
	return ConvertToLongOrThrow(token.GetTokenReference());
	}

template<typename T>
NumberCmdOption<T>::NumberCmdOption(
  const std::string &n, 
  T *manipVal, 
  T def, 
  T minV, 
  NumberSource maxV, 
bool persist, 
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<T>(n, manipVal, def, true, persist, pLevel),
	maximumProvider(maxV),
	minVal(minV)
	{
	}

template<typename T> 
NumberCmdOption<T>::NumberCmdOption(
  const std::string &n, 
  T *manipVal, 
  T def, 
  NumberSource minV, 
  T maxV, 
 bool persist, 
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<T>(n, manipVal, def, true, persist, pLevel),
	minimumProvider(minV),
	maxVal(maxV)
	{
	}

template<typename T> 
NumberCmdOption<T>::NumberCmdOption(
  const std::string &n, 
  T *manipVal, 
  T def, 
  T minV, 
  T maxV, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<T>(n, manipVal, def, false, persist, pLevel),
	maxVal(maxV),
	minVal(minV)
	{
	}

template<typename T> 
NumberCmdOption<T>::NumberCmdOption(
  const std::string &n, 
  T *manipVal, 
  T def, 
  NumberSource minV, 
  NumberSource maxV, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<T>(n, manipVal, def, false, persist, pLevel),
	maximumProvider(maxV),
	minimumProvider(minV)
	{
	}


template <>
inline std::string NumberCmdOption<unsigned>::GetDisplayType(
  bool includeIndefiniteArticle,
  bool plural) const
	{
	return StandardNoun("positive integer",includeIndefiniteArticle, plural);
	}

template <>
inline std::string NumberCmdOption<int>::GetDisplayType(bool includeIndefiniteArticle,
  bool plural) const
	{
	if (plural)
		return "integers";
	if (includeIndefiniteArticle)
		return "an integer";
	return "integer";
	}

template <>
inline std::string NumberCmdOption<long>::GetDisplayType(bool includeIndefiniteArticle,
  bool plural) const
	{
	if (plural)
		return "integers";
	if (includeIndefiniteArticle)
		return "an integer";
	return "integer";
	}

template <>
inline std::string NumberCmdOption<double>::GetDisplayType(bool includeIndefiniteArticle,
  bool plural) const
	{
	return StandardNoun("number",includeIndefiniteArticle, plural);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Reads number and makes sure that it is within the bounds of the setting
*/
template<typename T>
bool NumberCmdOption<T>::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	bool negateNumber = (equalsAlreadyRead && token.GetTokenReference() == '-');
	token.AlterTokenReading(NxsToken::kHyphenNotPunctuation);
	if (equalsAlreadyRead || NxsCmdOption::EatEqualsThenAdvance(token))
		{
		try	
			{
			*(this->currentValue) = ReadValueAfterEquals(token);
			if (negateNumber)
				*(this->currentValue) = -(*this->currentValue);
			}
		catch (NxsX_NumberIsTooLarge &)
			{
			FlagTooLarge();
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch (NxsX_NumberIsTooSmall &)
			{
			FlagTooSmall();
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch (NxsX_NotANumber &)
			{
			NxsCmdOption::FlagError(NxsCmdOption::unrecognized);
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch(...)
			{
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
		if (NxsCmdOption::WasValidRead())
			{
			++token;
			return true;
			}
		}
	else
		token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
	return false;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Reads number and makes sure that it is within the bounds of the setting
*/
template <>
inline bool NumberCmdOption<unsigned>::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	bool negateNumber = (equalsAlreadyRead && token.GetTokenReference() == '-');
	if (negateNumber)
		{
		NxsCmdOption::FlagError(NxsCmdOption::unrecognized, "-");
		return false;
		}
	token.AlterTokenReading(NxsToken::kHyphenNotPunctuation);
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		try	
			{
			*this->currentValue = ReadValueAfterEquals(token);
			}
		catch (NxsX_NumberIsTooLarge &)
			{
			FlagTooLarge();
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch (NxsX_NumberIsTooSmall &)
			{
			FlagTooSmall();
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch (NxsX_NotANumber &)
			{
			NxsCmdOption::FlagError(NxsCmdOption::unrecognized);
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		catch(...)
			{
			token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
			return false;
			}
		token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
		if (WasValidRead())
			{
			++token;
			return true;
			}
		}
	else
		token.AlterTokenReading(NxsToken::kHyphenIsPunctuation);
	return false;
	}
	


typedef NumberCmdOption<unsigned> 	UIntCmdOption; 
typedef NumberCmdOption<int> 		IntCmdOption; 
typedef NumberCmdOption<double> 	DblCmdOption; 
typedef NumberCmdOption<long> 		LongCmdOption; 

#endif
