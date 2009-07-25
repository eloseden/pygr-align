
# Author: Eden T. Elos, <eloseden@msu.edu>
# June 29,08
# ! /usr/bin/env python2.5

"""
BLAT_NLMSA MODULE
===================
A module that parses blat output from blat alignment program, and
builds pygr NLMSAs with it. The module defines the following class:

- `BlatLocalAlignment`, a blat gapped local alignment, consisting of
  multiple ungapped blocks.

Functions:

- `parse_blat()`: takes a blat alignment buffer and returns a list of
  BlastLocalAlignments and names of the sequences
- `build_blat_ivals():`: takes blat file buffer and sequence db
  as input and builds the ivals
- `create_NLMSA_blat()`: takes blat alignment file buffer, sequence db and NLMSA
  as input and returns a modified/built NLMSA

 

How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import blat_NLMSA``.
   You will also need to ``from pygr import cnestedlist, seqdb``.

2. Obtain the NLMSA using create_NLMSA_blat(buf, seqDb, al)
   function. One needs to pass blat output file object (buf), sequence
   database(seqDb) and the NLMSA object (al) to the function and the
   function returns the modified/built NLMSA.
   ``nlmsa_aln = create_NLMSA_blat(buf, seqDb, al)``

"""

__docformat__ = 'restructuredtext'

from pygr import cnestedlist, nlmsa_utils, seqdb, translationDB

# BlatLocalAlignment

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
    
# BlatUngappedBlock

class BlatUngappedBlock:
    """
    A single ungapped block in a blat alignment.
    """
    def __init__(self, qStart, qEnd, tStart, tEnd, orient):
        self.qStart = qStart
        self.qEnd = qEnd
        self.tStart = tStart
        self.tEnd = tEnd
        self.orient = orient

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
    and names of the sequences.
    """

    assert buf[0:8] == 'psLayout', " This is not a blat alignment file"
    
    records = buf.strip().split('\n')
    records = records[5:]
    
    seqs_names = set()
    matches = []

    # Each line/record in a blat alignment file contains a single
    # alignment to be represented by a BlatLocalAlignment object.
    # The elements/fields in a record are tab-separated
    # and the particular indices used is according to the specification
    # of the default blat alignment output.
    
    records = [ i.strip().split('\t') for i in records ]
    for record in records:
       # orientation information is not yet used in the NLMSA
       if len(record[8]) == 1:
           orient = record[8]+record[8]
       else:
           orient = record[8]
           
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
       seqs_names = seqs_names.union(set([qName, tName]))

       # construct a list of tuples with each tuple containing
       # the i^th ungapped blocks's coords
       # i.e. qStarts[i], tStarts[i], qEnds[i], tEnds[i]
       
       blocks = []  
       for i in range(0,len(qStarts)):
           blocks.append((qStarts[i], tStarts[i], qEnds[i], tEnds[i], orient))

       blocks = [ BlatUngappedBlock(a, c, b, d, ori) \
               for (a, b, c, d, ori) in blocks ]


       blatLocalAln = BlatLocalAlignment(qStart, qEnd, tStart, tEnd,
                                       qName, tName, orient, blocks)
       
       matches.append(blatLocalAln)       

    return matches, list(seqs_names)

def build_blat_ivals(buf):
    """
    Takes a blat file buffer as input and builds the ivals
    """
    blataln_list, seqs_names = parse_blat(buf)
    
    for blt_al in blataln_list:
        seqs_name1 = getattr(blt_al, "qSeqName")
        seqs_name2 = getattr(blt_al, "tSeqName")
        ivals = []   
        block = getattr(blt_al, "blocks")
        for ungapped in block:
            
            a = getattr(ungapped, "qStart")
            b = getattr(ungapped, "qEnd")
            
            x = getattr(ungapped, "tStart")
            y = getattr(ungapped, "tEnd")

            orient = getattr(ungapped, "orient")
            
            if orient[0] == '+':
                orient1 = 1
            else:
                orient1 = -1
            
            if orient[1] == '+':
                orient2 = 1
            else:
                orient2 = -1
            
            ival1 = (seqs_name1, a, b, orient1)
            ival2 = (seqs_name2, x, y, orient2)
            
            ivals.append((ival1, ival2))

        yield ivals

def create_NLMSA_blat(buf, seqDb,al):
    """
    Takes a blat alignment file buffer, sequence db and NLMSA (al) as input
    and returns a built NLMSA
    """
    for ivals in build_blat_ivals(buf):
        alignedIvalsAttrs = dict(id=0, start=1, stop=2, idDest=0, startDest=1,
                                 stopDest=2, ori=3, oriDest=3)
        cti = nlmsa_utils.CoordsToIntervals(seqDb, seqDb,
                                            alignedIvalsAttrs)
        al.add_aligned_intervals(cti(ivals))
        
    #build alignment
    al.build()
    return al
 
def create_NLMSA_tblat(buf, srcDB, destDB, al):
    """
    Takes a blat alignment file buffer, srcDB, destDB and NLMSA (al) as input
    and returns a built NLMSA
    """
    #srcDB = translationDB.get_translation_db(srcDB)
    destDB = translationDB.get_translation_db(destDB)
    for ivals in build_blat_ivals(buf):
        alignedIvalsAttrs = dict(id=0, start=1, stop=2, idDest=0, startDest=1,
                                 stopDest=2, ori=3, oriDest=3)        
        
        cti = nlmsa_utils.CoordsToIntervals(srcDB, destDB,
                                            alignedIvalsAttrs)
        al.add_aligned_intervals(cti(ivals))
        
    #build alignment
    al.build()
    return al