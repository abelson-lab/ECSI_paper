#!/usr/bin/env python

import re
from collections import defaultdict
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument("-i", action="store", dest="inSam", help="input sam with concensus reads i.e *PFDM_sorted_preDuplex.sam", required=True)
parser.add_argument("-outDuplex",  action="store", dest="outDuplex", help="output consensus duplex SAM file", required=True)
parser.add_argument('-c', type=float, default=.7, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
parser.add_argument('-n', action="store", dest='usedForDuplex', help="Sample Prefix")
parser.add_argument('--readlength', type=int, default=123, dest='readlength', help="Length of the input read that is being used. [123]")

args = parser.parse_args()
inSam = open( args.inSam , "r")
outDuplex = open( args.outDuplex , "w")
usedForDuplex = open( args.usedForDuplex+"_usedForDuplex" , "w")

class data:
        def __init__(self, **kwargs):
               self.__dict__.update(kwargs)


###############
### FUNCTIONS
###############
def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]

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
                        if float(nucIdentityList[k])/float(nucIdentityList[5]) > cutoff:
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
groupTag = []
preRead = []
groupFlag = []
for line in inSam:
        lineNO = lineNO + 1
        read = printRead(line)
	lastgroupTag = groupTag
	lastgroupFlag = groupFlag
        if lineNO == 1:
                FamilySize = int(read.family)
                groupEnd = lineNO + FamilySize
        if lineNO == groupEnd:
                FamilySize = int(read.family)
                groupEnd = lineNO + FamilySize
		for i in xrange(len(groupTag)-1):
			for j in range(i+1,len(groupTag),1):
				if split(groupTag[i],2)[0]==split(groupTag[j],2)[1] and split(groupTag[i],2)[1]==split(groupTag[j],2)[0] and abs(int(groupFlag[i]) - int(groupFlag[j]))==64 :
					consensus, newQual = consensusMaker([preRead[i].seq,preRead[j].seq], [preRead[i].qual,preRead[j].qual], args.cutoff, args.readlength)
					outDuplex.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (preRead[j].qname, preRead[j].flag, preRead[j].tid, preRead[j].pos, preRead[j].mapq, preRead[j].cigar, preRead[j].mrnm, preRead[j].mpos, preRead[j].isize, consensus, newQual))
					usedForDuplex.write("%s:%s:%s:%s:%s\n" % (preRead[i].qname,preRead[i].tid, preRead[i].pos,preRead[i].mpos, preRead[i].cigar))
					usedForDuplex.write("%s:%s:%s:%s:%s\n" % (preRead[j].qname,preRead[j].tid, preRead[j].pos,preRead[j].mpos, preRead[j].cigar))
		groupTag = []
		preRead = []
		groupFlag = []
        if lineNO < groupEnd:        
		groupTag.append(read.qname.split(":")[7].split("]")[1])
		groupFlag.append(read.flag)
        preRead.append(read)
	preFamilySize=FamilySize

##for the last consensus
for i in xrange(len(groupTag)-1):
                        for j in range(i+1,len(groupTag),1):
                                if split(groupTag[i],2)[0]==split(groupTag[j],2)[1] and split(groupTag[i],2)[1]==split(groupTag[j],2)[0] and abs(int(groupFlag[i]) - int(groupFlag[j]))==64 :
                                        consensus, newQual = consensusMaker([preRead[i].seq,preRead[j].seq], [preRead[i].qual,preRead[j].qual], args.cutoff, args.readlength)
                                        outDuplex.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (preRead[j].qname, preRead[j].flag, preRead[j].tid, preRead[j].pos, preRead[j].mapq, preRead[j].cigar, preRead[j].mrnm, preRead[j].mpos, preRead[j].isize, consensus, newQual))
					usedForDuplex.write("%s:%s:%s:%s:%s\n" % (preRead[i].qname,preRead[i].tid, preRead[i].pos,preRead[i].mpos, preRead[i].cigar))
					usedForDuplex.write("%s:%s:%s:%s:%s\n" % (preRead[j].qname,preRead[j].tid, preRead[j].pos,preRead[j].mpos, preRead[j].cigar))
					

inSam.close()
outDuplex.close()
usedForDuplex.close()

