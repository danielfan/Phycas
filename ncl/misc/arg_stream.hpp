#if !defined(ARG_STREAM_H)
#define ARG_STREAM_H
#include <string>
#include <vector>
class ArgStream
	{
	public:
		ArgStream(int argc, char *argv[]);
		bool ReadNextToken(std::string *s);
		
	private:
		typedef std::vector<std::string> VecToken;
		typedef VecToken::const_iterator VecToken_ConIt;
		
		VecToken 	Split(const std::string &s) const;

		VecToken	 	tokenStream;	
		VecToken_ConIt	tokenIt;	
	};
	
inline bool ArgStream::ReadNextToken(std::string *s)
	{
	if (tokenIt == tokenStream.end())
		return false;
	*s = *tokenIt++;
	return true;
	}
	
#endif
