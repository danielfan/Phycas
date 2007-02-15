/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*
These constants have been segregated here because they are used not
only in ReadNexus.idl but also in code that needs to be independent
of CORBA
*/

/* bits used in the blocksToRead arg of readNexusFile */
const short NEXUS_TAXA_BLOCK_BIT = 0x01;
const short NEXUS_TREES_BLOCK_BIT = 0x02;		
const short NEXUS_CHARACTERS_BLOCK_BIT = 0x04;	
const short NEXUS_ASSUMPTIONS_BLOCK_BIT = 0x08;	
const short NEXUS_SETS_BLOCK_BIT = 0x10;	

/* block to index mapping used in the NumBlockReadSequence sequence returned by readNexusFile */
const short NEXUS_TAXA_BLOCK_INDEX = 0;
const short NEXUS_TREES_BLOCK_INDEX = 1;		
const short NEXUS_CHARACTERS_BLOCK_INDEX = 2;	
const short NEXUS_ASSUMPTIONS_BLOCK_INDEX = 3;	
const short NEXUS_SETS_BLOCK_INDEX = 4;	
