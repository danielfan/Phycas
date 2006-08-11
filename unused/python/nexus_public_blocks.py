#!/usr/bin/python
from nexus_primitives import *
from nexus_parser import *
from nexus_taxa_block import *
from nexus_trees_block import *
from nexus_characters_block import *
import string
from sets import Set

_npb_identityTrans = string.maketrans('', '')
ALL_PUBLIC_BLOCKS = {}
def initializePublicBlocks():
	global ALL_PUBLIC_BLOCKS
	initializeCharDataBlock()
	ALL_PUBLIC_BLOCKS = {	'TAXA'	:	NexusTaxaBlock,
						'CHARACTERS':	NexusCharactersBlock,
						'DATA'		:	NexusDataBlock,
						'TREES'		:	NexusTreesBlock,
					}
initializePublicBlocks()
			