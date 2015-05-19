#!/usr/bin/env python

from __future__ import print_function

import re
import os
import os.path
import optparse
import pysam

# sort bam file with pysam.sort function
# the bam file should be sorted by qname

def SortSam(inBam, outBam):
	pysam.sort("-n", inBam, outBam)

# get tagmeth of single read

def TagMethOfSingleRead(read):

	return (pos, len, meth, unmeth, undt)

# get tagmeth of paired read

def TagMethOfPairedReads(read1, read2):

	return (pos, len, meth, unmeth, undt)

def main():

	# parse the command line options

	usage = 'usage: %prog [options] input.bam refseq.fa -o output.csv'
	parser = optparse.OptionParser(usage=usage, version='%prog 0.1.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input BAM file before data processing')

	(options, args) = parser.parse_args()
	if(len(args) != 2):
		parser.print_help()
		sys.exit(0)
	
	inputBamFileName = args[0]
	refSeqFileName = args[1]
	bamFileName = inputBamFileName
	baseFileName = os.path.splitext(os.path.basename(inputBamFileName))
	outputFileName =  baseFileName + '.tagmeth.csv'
	logFileName = baseFileName + '.log'

	# sort the input bam file if -s option is set

	rmTemp = False
	if(options.sort):
		print('[*] Sorting by QNAME...')
		bamFileName = 'sorted.' + os.path.basename(inputBamFileName)
		SortBam(inputBamFileName, bamFileName)
		rmTemp = True

	# load input files

	print('[*] Initializing...')

	if(not os.path.exists(bamFileName)):
		print('error: Failed to open file "', bamFileName, '"')
		sys.exit(-1)
	bamFile = pysam.AlignmentFile(bamFileName, "rb")
	
	if(not os.path.exists(refSeqFileName)):
		print('error: Reference sequence file "', refSeqFileName, '"', ' doest not exist.')
		sys.exit(-1)

	refSeq = {}
	with open(refSeqFileName, 'r') as refSeqFile :
		chrname = ''
		seq = ''
		for line in refSeqFile:
			if(line[0] == '>'):
				
				# save current seq for current chr

				if(chrname != ''):
					refSeq[chrname] = seq
				
				# new chrname & seq

				chrname = line[1:].strip()
				seq = ''
			else:
				seq += line.strip().upper()
	refSeqFile.close()

	# prepare output files

	if(options.outputfile):
		outputFileName = options.outputfile
	try:
		outFile = open(outputFileName, 'w')
	except IOError:
		print('error: Write to output file failed!')
		sys.exit(-1)

	# analyse algnments

	print('[*] Analyzing...')
	
	dictUnpaired = {}
	readCount = 0
	for read in bamFile.fetch():
		pairKey = ''.join(pairKeyRe.findall(read.qname)[0])
		chrname = bamFile.getrname(read.rname)

		# check if it is properly paired

		if(read.is_proper_pair):

			if(pairKey in dictUnpaired):
				
				# find paired read, handle paired reads

				tagmeth = TagMethOfPairedReads(read, dictUnpaired[pairKey])
				dictUnpaired.pop(pairKey, None)

			else:

				# not paired yet

				dictUnpaired[pairKey] = read
				tagmeth = None
		else:

			# handle single read

			tagmeth = TagMethOfSingleRead(read)

		if(tagmeth):
			outFile.write("%s\t%s\t%ld\t%d\t%s\t%s\t%s\n" % (
				pairKey, chrname, 
				tagmeth[0], tagmeth[1],
				tagmeth[2], tagmeth[3], tagmeth[4]))

		# progress

		readCount += 1
		sys.stdout.write('\r    read: #%ld' % (readCount))
		sys.stdout.flush()
	
	# release resources

	bamFile.close()
	outFile.close()
	
	if(rmTemp):
		os.unlink(bamFileName)

	print('[*] Complete')

if __name__ == '__main__':
	main()
