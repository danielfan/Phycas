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
