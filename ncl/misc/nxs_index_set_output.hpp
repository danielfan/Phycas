#if !defined(NXS_INDEX_SET_OUTPUT_HPP)
#define NXS_INDEX_SET_OUTPUT_HPP

#include "ncl/output/temp_basic_output_operators.hpp"
class NxsBasicListManager;
class NxsTaxaManager;
class NxsCharactersManager;
class NxsTreesManager;
class NxsIndexSet;

// could use an NxsIndexSetWrapper class for this
enum  NCLSetOutputType
		{
		kTaxaSetOutput,
		kCharactersSetOutput,
		kTreesSetOutput,
		kGenericSetOutput
		};

	
class NxsIndexSetOutput
	{
	public:
		NxsIndexSetOutput(const NxsIndexSet &s, const NxsBasicListManager & m, NCLSetOutputType typeOfSet)
			:indexSet(s),
			labelManager(m),
			setType(typeOfSet)
			{}
		const NxsIndexSet	  & indexSet;
		const NxsBasicListManager & labelManager;
		const NCLSetOutputType	setType;
	};

template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, NxsIndexSetOutput, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const NxsIndexSetOutput & iso)
			{
			const unsigned mgrSize = iso.labelManager.GetSize();
			const std::string & name = iso.indexSet.GetName();
			if (name.empty())
				outStream << "(unnamed set)";
			else
				outStream << "The set " << GetAsNexusToken(name);
			if (mgrSize == 0 || iso.indexSet.empty())
				{
				outStream <<" is empty";
				return;
				}
			outStream << " =";
			for (NxsIndexSet::const_iterator sIt = iso.indexSet.begin(); sIt != iso.indexSet.end(); ++sIt)
				{
				if (*sIt < mgrSize)
					outStream << ' ' << iso.labelManager.GetLabel(*sIt);
				else
					{
					outStream << " [and other out-of-range indices]";
					break;
					}
				}
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kConciseOutStyle, NxsIndexSetOutput, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const NxsIndexSetOutput & iso)
			{
			const unsigned mgrSize = iso.labelManager.GetSize();
			const std::string & name = iso.indexSet.GetName();
			if (name.empty())
				outStream << "(unnamed set)";
			else
				outStream << "The set " << GetAsNexusToken(name);
			if (mgrSize == 0 || iso.indexSet.empty())
				{
				outStream <<" is empty";
				return;
				}
			outStream << " =";
			outStream << iso.indexSet.GetNexusDescription();
			}
	};

	
	

#endif
