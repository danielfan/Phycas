#if !defined(CIPRES_NEXUS_READER_HPP)
#define CIPRES_NEXUS_READER_HPP

#include <string>
#include <stack>
#include "boost/shared_ptr.hpp"
#include "pyphy/src/ncl/nxs_reader.hpp"
#include "pyphy/src/cipres/CipresNativeC.h"

class PhoTaxaManager;
class PhoTreesManager;
class PhoCharactersManager;

#if defined(POL_PHYCAS)
struct my_exception // : std::exception
	{
	my_exception(std::string s) : msg(s) {}
	char const* what() const { return msg.c_str(); }
	std::string msg;
	};
#endif

class FullTreeDescription;
class NxsToken;
class CipresNexusReader: public NxsReader
	{
	public: 
		CipresNexusReader(const int blocksToRead);

		void 	NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult resCode, NxsBlock* b = NULL);
		void	ReadNxsFilePath(NxsInFilePath & filePath);
		void	ReadStringAsNexus(const std::string & s);
		void	ReadFilePath(const std::string &filePath)
			{
			NxsInFilePath fp(filePath);
			ReadNxsFilePath(fp);
			}

#		if defined(POL_PHYCAS)
		int		GetNChar();
		std::string GetErrorMessage()
			{
			return errorMsg;
			}
#		endif

		boost::shared_ptr<const PhoTaxaManager> GetTaxaManager() const
			{
			return phoTaxaMgr;
			}
		boost::shared_ptr<const PhoTreesManager> GetTreesManager() const
			{
			return phoTreesMgr;
			}
		boost::shared_ptr<const PhoCharactersManager> GetCharactersManager() const
			{
			return phoCharactersMgr;
			}
			
		const std::vector<FullTreeDescription> & GetTrees() const;

#	if defined(POL_PHYCAS) //@POL 10-Nov-2005 added
		std::vector<std::string> GetTaxLabels() const;
#	endif

	private:

		void	ReadNxsTokenStream(NxsToken & tokenStream);
		
		boost::shared_ptr<PhoTaxaManager> phoTaxaMgr;
		boost::shared_ptr<PhoTreesManager> phoTreesMgr;
		boost::shared_ptr<PhoCharactersManager> phoCharactersMgr;
		std::string errorMsg;
		file_pos	filePositionOfError;
		long		lineOfError;
		long		columnOfError;
		std::stack<std::string> fileStack;
	};

namespace ncl
	{
	class NxsDiscreteMatrix;
	}
namespace CipresNative
	{
	class DiscreteMatrix;
	}
	
CIPR_Matrix 	toCIPRMatrix(const ncl::NxsDiscreteMatrix & inMatrix, const bool verbose);
void 			freeCIPRMatrix(CIPR_Matrix & mat);
CipresNative::DiscreteMatrix *createNativeDiscreteMatrix(const CipresNexusReader & nexusReader, unsigned int charBlockIndex);
std::string		incrementTreesTaxaIndices(const std::string & zeroBasedNewick);

#endif
