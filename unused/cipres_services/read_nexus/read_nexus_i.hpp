#ifndef READNEXUSI_H_
#define READNEXUSI_H_

#include "boost/shared_ptr.hpp"
#include "CipresIDL/ReadNexusS.h"
#include "CipresCommlib/LifeCycle_i.h"

#include "CipresCommlib/CipresDataMatrixHelper.h"
//namespace ncl
//	{
//	class NxsDiscreteMatrix;
//	}
class CipresNexusReader;
	
//CIPR_Matrix toCIPRMatrix(const ncl::NxsDiscreteMatrix & inMatrix, const bool verbose);
//CipresNative::DiscreteMatrix * createNativeDiscreteMatrix(const CipresNexusReader & nexusReader, unsigned int charBlockIndex);


class  CipresIDL_ReadNexus_i :
  public virtual POA_CipresIDL::ReadNexus, 
  public virtual CipresIDL_LifeCycle_i 
	{
	public:
		CipresIDL_ReadNexus_i (const bool verbose);
	  
		virtual ~CipresIDL_ReadNexus_i (void);
	  
		virtual ::CipresIDL::ReadNexus::NumBlockReadSequence * readNexusFile (const char * filePath, CORBA::Short blocksToRead)
		  ACE_THROW_SPEC ((CORBA::SystemException, CipresIDL::NexusException));

		virtual ::CipresIDL::ReadNexus::NumBlockReadSequence * readStringAsNexus (const char * nexusContent, CORBA::Short blocksToRead)
		  ACE_THROW_SPEC ((CORBA::SystemException, CipresIDL::NexusException));
	
		virtual ::CipresIDL::DataMatrix * getCharacters (CORBA::Short charactersBlockIndex)
		  ACE_THROW_SPEC ((CORBA::SystemException));
		virtual ::CipresIDL::TaxonInfoSeq * getTaxa (CORBA::Short taxaBlockIndex)
		  ACE_THROW_SPEC ((CORBA::SystemException));
		virtual ::CipresIDL::TreeSeq * getTrees (CORBA::Short treesBlockIndex)
		  ACE_THROW_SPEC ((CORBA::SystemException));
		virtual char * getUIXml ()
		  ACE_THROW_SPEC ((CORBA::SystemException));
		virtual CORBA::Boolean execute (const char * command, CORBA::String_out display)
		  ACE_THROW_SPEC ((CORBA::SystemException));

	protected:
		typedef boost::shared_ptr<CipresNexusReader> CipresNexusReaderShPtr;
		
		CipresNexusReaderShPtr newNexusReader(CORBA::Short blocksToRead) const;
		::CipresIDL::ReadNexus::NumBlockReadSequence * postNexusReading(CipresNexusReaderShPtr, CORBA::Short blocksToRead);
		
		CipresNexusReaderShPtr	nexusReader;
		bool					verboseMode;
	};


#endif /* READNEXUSI_H_  */
// -*- C++ -*-
//
// $Id: read_nexus_i.hpp,v 1.7 2005/09/07 20:22:14 mholder Exp $

// ****  Code generated by the The ACE ORB (TAO) IDL Compiler ****
// TAO and the TAO IDL Compiler have been developed by:
//       Center for Distributed Object Computing
//       Washington University
//       St. Louis, MO
//       USA
//       http://www.cs.wustl.edu/~schmidt/doc-center.html
// and
//       Distributed Object Computing Laboratory
//       University of California at Irvine
//       Irvine, CA
//       USA
//       http://doc.ece.uci.edu/
// and
//       Institute for Software Integrated Systems
//       Vanderbilt University
//       Nashville, TN
//       USA
//       http://www.isis.vanderbilt.edu/
//
// Information about TAO is available at:
//     http://www.cs.wustl.edu/~schmidt/TAO.html