#ifndef PHYCAS_TREE_MERGE_I_HPP
#define PHYCAS_TREE_MERGE_I_HPP

#include "CipresIDL/TreeMergeS.h"
#include "CipresCommlib/LifeCycle_i.h"

#if !defined (ACE_LACKS_PRAGMA_ONCE)
#	pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

class  CipresIDL_TreeMerge_i : 
	public virtual POA_CipresIDL::TreeMerge,
	public virtual CipresIDL_LifeCycle_i 
	{
	public:
		CipresIDL_TreeMerge_i(const bool verbose)
			:verboseMode(verbose)
			{}
		virtual ~CipresIDL_TreeMerge_i (void)
			{}
		
		virtual ::CipresIDL::Tree * mergeTrees (const CipresIDL::TreeSeq & trees)
			ACE_THROW_SPEC ((CORBA::SystemException));
		virtual char * getUIXml ()
			ACE_THROW_SPEC ((CORBA::SystemException));
		virtual CORBA::Boolean execute (const char * command, CORBA::String_out display)
			ACE_THROW_SPEC ((CORBA::SystemException));
	protected:
		bool verboseMode;
	};

#endif 
// -*- C++ -*-
//
// $Id: tree_merge_i.hpp,v 1.2 2005/06/07 02:52:10 mholder Exp $

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