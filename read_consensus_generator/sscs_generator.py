#!/usr/bin/env python

import re
from collections import defaultdict
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument("-i", action="store", dest="inSam", help="input *_PFFC_sorted_Part_of_family_of_SizeX.sam", required=True)
parser.add_argument("-outSin",  action="store", dest="outSin", help="output consensus SAM file for families of size 1 (singletons only)", required=True)
parser.add_argument("-out2OM",  action="store", dest="out2OM", help="output consensus SAM file for families of size 2 or more", required=True)
parser.add_argument('-c', type=float, default=1, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [1]")
parser.add_argument('--readlength', type=int, default=123, dest='readlength', help="Length of the input read that is being used. [123]") #THE CODE ASSUMES THAT ALL THE READS ARE AT THE SAME LENGTH
parser.add_argument('-n', type=float, default=0.05, dest='Nthreshold', help="Percentage of N nucleotides in SSCS")

args = parser.parse_args()
inSam = open( args.inSam , "r")
outSin = open( args.outSin , "w")
out2OM = open( args.out2OM , "w")

class data:
        def __init__(self, **kwargs):
               self.__dict__.update(kwargs)

###############
### FUNCTIONS
###############
def printRead(readIn):
        parts=readIn.split()
        read = data(family=parts[0], exTag=parts[1], qname=parts[2], flag=parts[3], tid=parts[4], pos=parts[5], mapq=parts[6], cigar=parts[7], mrnm=parts[8], mpos=parts[9], isize=parts[10], seq=parts[11], qual=parts[12])
        return read

def consensusMaker (groupedReadsList, groupedQualsList,  cutoff,  readLength) :
        nucIdentityList=[0, 0, 0, 0, 0, 0] # In the order of T, C, G, A, N, Total
        nucKeyDict = {0:'T', 1:'C', 2:'G', 3:'A', 4:'N'}
        consensusRead = ''
        consensusQual = ''
        qualAscii=[0, 0, 0, 0, 0]
        for i in xrange(readLength) :
                for j in xrange(len(groupedReadsList)):
                        if groupedReadsList[j][i] == 'T' :
                                nucIdentityList[0] += 1
                                if ord(groupedQualsList[j][i])>qualAscii[0]:
                                        qualAscii[0] = ord(groupedQualsList[j][i])
                        elif groupedReadsList[j][i] == 'C':
                                nucIdentityList[1] += 1
                                if ord(groupedQualsList[j][i])>qualAscii[1]:
                                        qualAscii[1]=ord(groupedQualsList[j][i])
                        elif groupedReadsList[j][i] == 'G':
                                nucIdentityList[2] += 1
                                if ord(groupedQualsList[j][i])>qualAscii[2]:
                                        qualAscii[2]=ord(groupedQualsList[j][i])
                        elif groupedReadsList[j][i] == 'A':
                                nucIdentityList[3] += 1
                                if ord(groupedQualsList[j][i])>qualAscii[3]:
                                        qualAscii[3]=ord(groupedQualsList[j][i])
                        elif groupedReadsList[j][i] == 'N':
                                nucIdentityList[4] += 1
                                if ord(groupedQualsList[j][i])>qualAscii[4]:
                                        qualAscii[4]=ord(groupedQualsList[j][i])
                        else:
                                nucIdentityList[4] += 1
                        nucIdentityList[5] += 1
                for k in [0, 1, 2, 3, 4]:
                        if float(nucIdentityList[k])/float(nucIdentityList[5]) >= cutoff:
                                consensusRead += nucKeyDict[k]
                                consensusQual += chr(qualAscii[k])
                                break
                        elif k==4:
                                consensusRead += 'N'
                                consensusQual += "!"
                nucIdentityList=[0, 0, 0, 0, 0, 0] # Reset for the next nucleotide position
                qualAscii=[0, 0, 0, 0, 0, 0]

        return consensusRead, consensusQual;


##########
### MAIN
##########
FamilySize = 0
lineNO = 0
groupSEQ = []
groupQUAL = []
for line in inSam:
        lineNO = lineNO + 1
        read = printRead(line)
        if lineNO == 1:
                FamilySize = int(read.family)
                groupEnd = lineNO + FamilySize
        if lineNO == groupEnd:
                FamilySize = int(read.family)
                groupEnd = lineNO + FamilySize
                if preFamilySize == 1 :
                	outSin.write("%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("S_",preRead.qname, preRead.flag, preRead.tid, preRead.pos, preRead.mapq, preRead.cigar, preRead.mrnm, preRead.mpos, preRead.isize, preRead.seq, preRead.qual))
                elif preFamilySize >= 2 :
                	consensus, newQual = consensusMaker(groupSEQ, groupQUAL, args.cutoff, args.readlength)
			if float(consensus.count('N'))/args.readlength<args.Nthreshold :
                        	out2OM.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (preRead.qname, preRead.flag, preRead.tid, preRead.pos, preRead.mapq, preRead.cigar, preRead.mrnm, preRead.mpos, preRead.isize, consensus, newQual))
                groupSEQ = []
                groupQUAL = []
        if lineNO < groupEnd:        
		groupSEQ.append(read.seq)
                groupQUAL.append(read.qual)
        preRead=read
	preFamilySize=FamilySize

##for the last consensus
if preFamilySize == 1 :
	outSin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (preRead.qname, preRead.flag, preRead.tid, preRead.pos, preRead.mapq, preRead.cigar, preRead.mrnm, preRead.mpos, preRead.isize, preRead.seq, preRead.qual))
elif preFamilySize >= 2 :
	consensus, newQual = consensusMaker(groupSEQ, groupQUAL, args.cutoff, args.readlength)
	if float(consensus.count('N'))/args.readlength<args.Nthreshold :
	        out2OM.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (preRead.qname, preRead.flag, preRead.tid, preRead.pos, preRead.mapq, preRead.cigar, preRead.mrnm, preRead.mpos, preRead.isize, consensus, newQual))

inSam.close()
outSin.close()
out2OM.close()
