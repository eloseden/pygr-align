# Author Eden Elos
# June 29,08
# create NLMSA from Blat alignment file
# ! /usr/bin/env python2.5

import os
import glob
from pygr import cnestedlist, seqdb


#
# BlatLocalAlignment
#

class BlatLocalAlignment:
    """
    A blat gapped local alignment, consisting of multiple
    ungapped blocks.
    """
    def __init__(self, qStart, qEnd, tStart, tEnd, qSeqName, tSeqName,
                 orient, blocks):
        
        self.qStart = qStart
        self.qEnd = qEnd
        self.tStart = tStart
        self.tEnd = tEnd

        self.qSeqName = qSeqName
        self.tSeqName = tSeqName
        
        self.orient = orient

        self.blocks = blocks
    

#
# BlatUngappedBlock
#

class BlatUngappedBlock:
    """
    A single ungapped block in a blat alignment.
    """
    def __init__(self, qStart, qEnd, tStart, tEnd):
        self.qStart = qStart
        self.qEnd = qEnd
        self.tStart = tStart
        self.tEnd = tEnd



    def __len__(self):
        return self.tEnd - self.tStart

    def convert_to_text(self, seq1, seq2):
        """
        Return the ungapped alignment query and target as sequence.
        """
        q = seq1[self.qStart:self.qEnd]
        t = seq2[self.tStart:self.tEnd]

        return (q, t)



def calculate_end(Starts, blockSize):
    """
    Calculate the end coordinates of ungapped blocks
    """
    
    Ends = []
          
    for i in range(0, len(blockSize)):
        Ends.append(int(Starts[i]) + int(blockSize[i]))

    return Ends
       
        

def parse_blat(buf):
    """
    Takes a blat alignment buffer and returns a list of BlastLocalAlignments
    and names of the genomes.
    """

    assert buf[0:8] == 'psLayout', " This is not a blat alignment file"
    
    records = buf.strip().split('\n')
    records = records[5:]
    
    genome_names = set()
    matches = []

    # Each line/record in a blat alignment file contains a single
    # alignment to be represented by a BlatLocalAlignment object.
    # The elements/fields in a record are tab-separated
    # and the particular indices used is according to the specification
    # of the default blat alignment output.
    
    records = [ i.strip().split('\t') for i in records ]
    for record in records:
       # orientation information is not yet used in the NLMSA
       if record[8] == '+':
           orient = 1
       else:
           orient = -1
           
       qName = record[9]
       tName = record[13]
       qStart = int(record[11])
       tStart = int(record[15])
       qEnd = int(record[12])
       tEnd = int(record[16])
       blockCount = int(record[17])
       blockSize = map(int, record[18].strip(',').split(','))
       qStarts = map(int, record[19].strip(',').split(','))       
       tStarts = map(int, record[20].strip(',').split(','))
       qEnds = map(int, calculate_end(qStarts, blockSize))
       tEnds = map(int, calculate_end(tStarts, blockSize))
       genome_names = genome_names.union(set([qName, tName]))

       # construct a list of tuples with each tuple containing
       # the i^th ungapped blocks's coords
       # i.e. qStarts[i], tStarts[i], qEnds[i], tEnds[i]
       blocks = []  
       for i in range(0,len(qStarts)):
           blocks.append((qStarts[i], tStarts[i], qEnds[i], tEnds[i]))

       blocks = [ BlatUngappedBlock(a, c, b, d) \
               for (a, b, c, d) in blocks ]


       blatLocalAln = BlatLocalAlignment(qStart, qEnd, tStart, tEnd,
                                       qName, tName, orient, blocks)
       
       matches.append(blatLocalAln)       

    return matches, list(genome_names)

def create_NLMSA_blat(buf, seqDb, al):
    """
    takes blat file buffer as input and creates and returns NLMSA
    """
    blataln_list, genome_names = parse_blat(buf)
    
    #feed the genomes
    for i in range(0,len(genome_names)):
        al += seqDb[genome_names[i]]


    for blt_al in blataln_list:
        genome_name1 = getattr(blt_al, "qSeqName")
        genome_name2= getattr(blt_al, "tSeqName")
            
        block = getattr(blt_al, "blocks")
        for ungapped in block:
            
            a = getattr(ungapped, "qStart")
            b = getattr(ungapped, "qEnd")
            
            x = getattr(ungapped, "tStart")
            y = getattr(ungapped, "tEnd")

            ival1 = seqDb[genome_name1][a:b]
            ival2 = seqDb[genome_name2][x:y]

            
            al[ival1] += ival2

    # build alignment
    al.build()
    return al
            
        


