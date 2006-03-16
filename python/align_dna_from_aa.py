#!/usr/bin/python
import sys
from nexus_public_blocks import *

usingBBEdit = True # True to use bbedit to highlight the lines of Nexus Errors

class NoCorrespondingKeyError(ValueError): pass

def reorderList(listToChange, currKeyOrder, desiredKeyOrder):
	for i in range(len(desiredKeyOrder)):
		nextKey = desiredKeyOrder[i]
		try:
			j = currKeyOrder.index(nextKey)
		except ValueError: raise NoCorrespondingKeyError, nextKey
		if i != j:
			listToChange[i], listToChange[j] = listToChange[j], listToChange[i]
			currKeyOrder[i], currKeyOrder[j] = currKeyOrder[j], currKeyOrder[i]

def findFrameAndCrop(dna, prot, geneticCodeTranslator):
	'''canBeSuperSubSequence => 1 bit for "translation in prot", 2 bit for "prot in translation", neither for translation.isMoreResolvedSubset(prot)'''
	for i in range(3):
		translation = geneticCodeTranslator.translate(dna, False)
		sp = None
		sp = translation.findSubSequenceStartPoint(prot)
		if sp != None:
			firstNuc = 3 * sp
			if firstNuc > 0:
				del dna[:firstNuc]
			toCropInd = 3*len(prot) - len(dna)
			if toCropInd < 0:
				del dna[toCropInd:]
			return True
		dna.insertCode(0, DNASeq.anyNonGapCode)
	del dna[0:3] # get rid of the 2 N's that we added
	return False
def printNoFrame(out, dna, prot, translator):
	out.write('using the %s genetic code\n' % translator.name)
	translations = []
	for i in range(3):
		translations.append(str(translator.translate(dna)))
		dna.insertCode(0, DNASeq.anyNonGapCode)
	out.write('desired aa   %s\n' % '  '.join([i for i in str(prot)]))
	out.write('frame 1       %s\n' % '  '.join([i for i in translations[0]]))
	out.write('frame 2      %s\n' % '  '.join([i for i in translations[1]]))
	out.write('frame 3     %s\n' % '  '.join([i for i in translations[2]]))
	out.write('dna          %s\n' % ''.join([i for i in str(dna)[2:]]))
	
if __name__ == '__main__':
	if len(sys.argv) != 4:
		print >>sys.stderr, 'Usage: %s <aligned protien nexus file> <dna nexus file> <Standard|VertMito>' % sys.argv[0]
		sys.exit(1)
	import genetic_code
	mostLikelyTrans = 0
	requestedCode = sys.argv[3].upper()
	if requestedCode == 'STD' or requestedCode.startswith('STAND'):
		requestedCode = 'STANDARD'
	elif requestedCode.startswith('VERTMIT') or requestedCode.startswith('VMIT'):
		requestedCode = 'VERTMITO'
	translators = copy.copy(genetic_code.genetic_codes) 
	for t in translators:
		if t.name.upper() == requestedCode:
			translators = [t]
			break
	NexusParsing.statusStream = sys.stderr
	if len(translators) != 1:
		print >> sys.stderr, 'The genetic code %s was not found.  Valid choices are:\n\t%s\n' %(requestedCode, '\n\t'.join([i.name for i in translators]))
		sys.exit(10)
	else:
		NexusParsing.statusMessage('Using the %s genetic code' % translators[0].name)
	
	matrices = []
	for inFilename in sys.argv[1:-1]:
		try:	
			NexusParsing.pushCurrentFile(inFilename)
			inF = open(inFilename)
			NexusParsing.statusMessage('Reading %s ...' % inFilename)
			for b in  NexusBlockStream(inF, ALL_PUBLIC_BLOCKS, True):
				if isinstance(b, NexusCharactersBlock):
					matrices.append(b)
					break
			NexusParsing.popCurrentFile()
			inF.close()
		except NexusError, x:
			import os
			print >>sys.stderr,  x
			if usingBBEdit and NexusParsing.currentFile != '':
				if x.startPos != None:
					s = 'bbedit +%d %s' % (x.startPos, NexusParsing.currentFile)
					#s = 'bbedit +%d %s' % (x.startPos.line, NexusParsing.currentFile)
				else:
					s = 'bbedit %s' % (NexusParsing.currentFile)
				os.system(s)
				sys.exit(1)
			NexusParsing.statusStream = None
	
	assert(len(matrices) == 2)
	dnaMat, protMat = matrices[0].dataType == 'DNA' and (matrices[0], matrices[1]) or (matrices[1], matrices[0])
	if dnaMat.dataType != 'DNA': raise ValueError, 'Expecting at least one matrix of DNA Sequences'
	if protMat.dataType != 'PROTEIN': raise ValueError, 'Expecting at least one matrix of PROTEIN Sequences'
	
	protSeqList = protMat.matrix
	dnaSeqList = dnaMat.matrix
	taxList = dnaMat.taxLabels
	 # make sure that the matrix lists are in same order.  We'll use the dna matrix order because that is the matrix we'll print
	try:
		reorderList(protSeqList, protMat.taxLabels, taxList)
	except NoCorrespondingKeyError, k:
		print >>sys.stderr, 'No protein sequence was found for the taxon %s' % k
		sys.exit(1)
	ntaxa = len(dnaSeqList)
	if ntaxa == 0: NexusParsing.statusMessage('Found no DNA Sequences in the Data/Characters Block.\nDone.')
	
	 # remove gaps from dna sequences, and determine their frame 
	 # we just try codes until we get a match (though we assume that all of the sequences should be translatable 
	 # using the same genetic code
	
	inFrameDNASeqs = []
	problems = []
	nDropped = 0
	for i in range(len(dnaSeqList)):
		NexusParsing.statusMessage('...Finding frame for %s...\n' % taxList[i - nDropped])
		gaplessAA = copy.copy(protSeqList[i - nDropped])
		gaplessAA.removeGaps()
		dnaWithGaps = copy.copy(dnaSeqList[i])
		gaplessDNA = dnaSeqList[i]
		gaplessDNA.removeGaps()
		framed = None
		for t in translators:
			stPhase = findFrameAndCrop(gaplessDNA, gaplessAA, t) # find the reading frame allowing for either sequence to be a subset of the other
			if stPhase == True:
				framed = gaplessDNA
				break
				# print 'Rejecting the genetic code %s from further consideration...' % translatorToTry.name
				# translators.pop(0)
		if framed == None:
			problems.append((taxList[i - nDropped], gaplessDNA, gaplessAA, dnaWithGaps))
			print >> sys.stderr, 'Could not find frame for which the DNA sequence for %s gave the specified Amino acid sequence!\n' % taxList[i - nDropped],
			taxList.pop(i - nDropped)
			protSeqList.pop(i - nDropped)
			nDropped += 1
		else:
			inFrameDNASeqs.append(framed)
		del dnaWithGaps
	if len(inFrameDNASeqs) == 0:
		print >>sys.stderr, 'No dna sequences were successfully matched to the amino acid sequences.\nExiting.'
		sys.exit(5)
	NexusParsing.statusMessage('Sequences are consistent with the %s genetic code...' % t.name)
	nAlignedDNCChar = 3 * len(protSeqList[0])
	 # transfer --- for every - in the Amino acid alignment (and crop the DNA sequences)
	 # BUG We don't crop the starting part of the matrix yet !!!
	for i in range(len(inFrameDNASeqs)):
		NexusParsing.statusMessage('...Adding gaps to the DNA sequence for %s...\n' % taxList[i])
		protSeq = protSeqList[i]
		currNtPos = 0
		dnaSeq = inFrameDNASeqs[i]
		for c in protSeq:
			if currNtPos >= len(dnaSeq):
				break
			if c == AASeq.gapCode:
				dnaSeq.insertCode(currNtPos, DNASeq.gapCode)
				dnaSeq.insertCode(currNtPos, DNASeq.gapCode)
				dnaSeq.insertCode(currNtPos, DNASeq.gapCode)
			currNtPos += 3
		nTerminalGaps = nAlignedDNCChar - len(dnaSeq)
		if nTerminalGaps > 0: 
			for j in range(nTerminalGaps):
				dnaSeq.appendCode(DNASeq.missingCode)	# pad matrix with missingCode (?)
		elif nTerminalGaps < 0: 
			del dnaSeq[nTerminalGaps:]					# crop sites for which we don't have an amino acid alignment
		
		#replace the dna matrix from the read file, and update nchar
	dnaMat.matrix = inFrameDNASeqs	
	dnaMat.nchar = len(dnaMat.matrix[0])
	dnaMat.ntax = len(taxList)
	
	outStream = sys.stdout
	outStream.write('#NEXUS\n')
	dnaMat.writeNexusBlock(sys.stdout)
	probOut = open('problems.txt', 'w')
	if len(problems) > 0:
		for p in problems:
			print >>probOut, '\nProblems with the sequences for %s:' % p[0],
			printNoFrame(probOut, p[1], p[2], translators[mostLikelyTrans])
			print >> probOut, '\n'
				
		dnaMat.taxLabels = [i[0] for i in problems]
		dnaMat.matrix = [i[3] for i in problems]
		dnaMat.nchar = len(dnaMat.matrix[0])
		dnaMat.ntax = len(dnaMat.taxLabels)
		probOut.write('\n\n\n#NEXUS\n')
		dnaMat.writeNexusBlock(probOut)
	