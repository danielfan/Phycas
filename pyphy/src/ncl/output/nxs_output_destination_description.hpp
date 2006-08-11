#if ! defined(NCL_NXS_OUTPUT_DESTINATION_DESCRIPTION_HPP)
#define NCL_NXS_OUTPUT_DESTINATION_DESCRIPTION_HPP
#include "pyphy/src/ncl/misc/nxs_file_path.hpp"

class NxsOutputDestinationDescription
	{
	public:
		enum OutputDestinationType
			{
			kSuppressOutput,
			kDirectOutputToNCLStream,
			kDirectOutputToFile,
			kDirectOutputToBoth
			};
		enum ValidNCLOutputRedirect
			{
			kRedirectToOutputStream		= 0,
			kRedirectToCommentStream	= 1,
			kRedirectToErrorStream		= 2,
			kRedirectToPlotStream		= 3
			};
		STATIC_CONST_ACCESSOR std::string TranslateRedirectionStreamName(ValidNCLOutputRedirect);
		STATIC_CONST_ACCESSOR const VecString 	& getLegalRedirectionStreamsNames();
		
		NxsOutputDestinationDescription();
		NxsOutputDestinationDescription(ValidNCLOutputRedirect o);
		NxsOutputDestinationDescription(const std::string &fn, bool appendIfPresent = false, bool replaceIfPresent = false);
		NxsOutputDestinationDescription(ValidNCLOutputRedirect o, const std::string &fn, bool appendIfPresent = false, bool replaceIfPresent = false);
		NxsOutputDestinationDescription &operator=(const NxsOutputDestinationDescription &other);
		bool IsSuppressOutput() const
			{
			return outDestinationType == kSuppressOutput;
			}
		bool IsRedirectToNCLStream() const
			{
			return outDestinationType == kDirectOutputToNCLStream ||  outDestinationType == kDirectOutputToBoth;
			}
		bool IsRedirectToFile() const
			{
			return outDestinationType == kDirectOutputToFile ||  outDestinationType == kDirectOutputToBoth;
			}
		std::string GetRedirectionStreamName(unsigned i) const
			{
			return TranslateRedirectionStreamName(redirectStreams[i]);
			}
		std::string GetNexusDescription() const;
		
		typedef std::vector<ValidNCLOutputRedirect> VecRedirection;
		
		OutputDestinationType   outDestinationType;
		VecRedirection			redirectStreams;
		mutable NxsOutFilePath  filePath; 
	};

#endif
