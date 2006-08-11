#if defined (CORBA_CLIENT_PHYCAS)

#include "phycas/corba/put_tree_wrapper.hpp"
#include "phycas/corba/corba_wrapper.hpp"
#include "phycas/trees/tree.hpp"
#include "ncl/misc/string_extensions.hpp"
using cipresCORBA::PutTreeWrapper;

///@move to cipresCORBA/utils.cpp
///@duplication necessary only because typedef is repeated in our crude IDL
cipresCORBA::TreeDrawer::stringArray_var toCORBA_Drawer(const StrVec & strVec)
	{
	cipresCORBA::TreeDrawer::stringArray_var labelArr = new cipresCORBA::TreeDrawer::stringArray();
	const CORBA::ULong len = strVec.size();
	labelArr->length(len);
	for (CORBA::ULong i = 0U; i < len; ++i)
		labelArr[i] = CORBA::string_dup(strVec[i].c_str());
	return labelArr._retn();
	}

// duplication necessary only because typedef is repeated in our crude IDL
cipresCORBA::TreeDatabase::stringArray_var toCORBA_Database(const StrVec & strVec)
	{
	cipresCORBA::TreeDatabase::stringArray_var labelArr = new cipresCORBA::TreeDatabase::stringArray();
	const CORBA::ULong len = strVec.size();
	labelArr->length(len);
	for (CORBA::ULong i = 0U; i < len; ++i)
		labelArr[i] = CORBA::string_dup(strVec[i].c_str());
	return labelArr._retn();
	}
	

template<typename TVAR>
class RemoteTreeAcceptorFunctor
	{
	public:
		RemoteTreeAcceptorFunctor()
			:useCORBAService(false)
			{}
		RemoteTreeAcceptorFunctor(cipresCORBA::CorbaWrapper &, const VecString &taxLabels);
		void operator()(const Tree &);
	protected:
		TVAR treeAcceptor;
		bool useCORBAService;
		std::string taxNameContext;
	};
typedef RemoteTreeAcceptorFunctor<cipresCORBA::TreeDrawer_var> TreeDrawerWrapper;
typedef RemoteTreeAcceptorFunctor<cipresCORBA::TreeDatabase_var> TreeDatabaseWrapper;

template<>
void RemoteTreeAcceptorFunctor<cipresCORBA::TreeDrawer_var>::operator()(const Tree & tree) 
	{
	std::string newick;
	tree.AppendNewickRepresentation(newick, true, true); // use tax labels and edge lengths
	newick << ';';
	try	{
		if (useCORBAService)
			treeAcceptor->drawTree(newick.c_str()); // removed CORBA::string_dup
		}
	catch (CORBA::SystemException & x)
		{
		std::cerr << "CORBA::SystemException from drawTree(), disabling corbaTreeDrawer" << std::endl;
		useCORBAService = false;
		}
	catch (cipresCORBA::TreeDrawer::X_InvalidTree & y)
		{
		std::cerr << "cipresCORBA::X_InvalidTree from drawTree(), disabling corbaTreeDrawer (and tripping assert in debug mode)" << std::endl;
		useCORBAService = false;
		assert(false);
		}
#	if defined (DEBUG_PUT_TREE_OUTPUT)
		std::cerr << newick << std::endl; //@temp, debugging
#	endif
	}

template<>
RemoteTreeAcceptorFunctor<cipresCORBA::TreeDrawer_var>:: RemoteTreeAcceptorFunctor(
  cipresCORBA::CorbaWrapper & corbaWrapper, 
  std::vector<std::string> const &taxLabels)
	:useCORBAService(true)
	{
#   if defined(USING_FACTORY)
		try {
			treeAcceptor = corbaWrapper.treeDrawerFactory->build();
			}
		catch (...)
			{
			std::cerr << "Factory did not build a treedrawer"<< std::endl;
			throw;
			}
#   endif
	cipresCORBA::TreeDrawer::stringArray_var taxLabelArr = toCORBA_Drawer(taxLabels);
	taxNameContext = CORBA::string_dup(treeAcceptor->findOrSetTaxa(*taxLabelArr));
	}
		
template<>
void  RemoteTreeAcceptorFunctor<cipresCORBA::TreeDatabase_var>::operator()(const Tree & tree) 
	{
	std::string newick;
	tree.AppendNewickRepresentation(newick, true, true); // use tax labels but no edge lengths
	newick << ';';
	try	{
		if (useCORBAService)
			treeAcceptor->putTree(taxNameContext.c_str(), newick.c_str()); // removed CORBA::string_dup
		}
	catch (CORBA::SystemException & x)
		{
		std::cerr << "CORBA::SystemException from putTree(), disabling corbaTreeDatabase" << std::endl;
		useCORBAService = false;
		}
#	if defined (DEBUG_PUT_TREE_OUTPUT)
		std::cerr << newick << std::endl; //@temp, debugging
#	endif
	}

template<>	
RemoteTreeAcceptorFunctor<cipresCORBA::TreeDatabase_var>::RemoteTreeAcceptorFunctor(
  cipresCORBA::CorbaWrapper & corbaWrapper, 
  std::vector<std::string> const &taxLabels)
	:useCORBAService(true)
	{
#   if defined(USING_FACTORY)
		try {
			treeAcceptor = corbaWrapper.treeDBFactory->build();
			}
		catch (...)
			{
			std::cerr << "Factory did not build a treeDatabase"<< std::endl;
			throw;
			}
#   endif
	cipresCORBA::TreeDatabase::stringArray_var taxLabelArr = toCORBA_Database(taxLabels);
	taxNameContext = CORBA::string_dup(treeAcceptor->findOrSetTaxa(*taxLabelArr));
	}
	
PutTreeWrapper::PutTreeWrapper(const VecString &taxLabels)
	{
	cipresCORBA::CorbaWrapper * corbaWrapper = NULL;
	try
		{
		corbaWrapper = cipresCORBA::CorbaWrapper::GetInstance();
		}
	catch (CORBA::SystemException &x)
		{
		std::cerr << "CORBA::SystemException from CorbaWrapper::GetInstance(), CORBA calls disabled" << std::endl;
		}
	if (corbaWrapper)
		{
		try
			{
			vecTreeAcceptorFunc.push_back(TreeAcceptorFunc(TreeDrawerWrapper(*corbaWrapper, taxLabels)));
			}
		catch (CORBA::SystemException &x)
			{
			std::cerr << "CORBA::SystemException from TreeDrawer Initialization.  CORBA TreeDrawing disabled" << std::endl;
			}
		try
			{
			vecTreeAcceptorFunc.push_back(TreeAcceptorFunc(TreeDatabaseWrapper(*corbaWrapper, taxLabels)));
			}
		catch (CORBA::SystemException &x)
			{
			std::cerr << "CORBA::SystemException from TreeDatabase Initialization.  CORBA treeDatabase storage disabled" << std::endl;
			}
		}
	else
		std::cerr << "No CorbaWrapper instance is available, trees will not be drawn or deposited using CORBA." << std::endl;
	}

void PutTreeWrapper::PutTree(const Tree & tree)
	{
	for (VecTreeAcceptorFunc::iterator tfIt = vecTreeAcceptorFunc.begin(); tfIt != vecTreeAcceptorFunc.end(); ++tfIt)
		{
		try	{
			(*tfIt)(tree);
			}
		catch (CORBA::SystemException & )
			{
			}
		}
	}


	


#endif


