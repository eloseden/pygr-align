
# Author: Eden T. Elos, <eloseden@msu.edu>
# July,01, 08
# ! /usr/bin/env python2.5

"""
MLAGAN_NLMSA MODULE
===================
A module that parses mlagan output from mlagan alignment program, and
builds pygr NLMSAs with it. The module does not define any class.

Functions:

- `read_mlagan()`: read aligned sequences from a mlagan alignment file
  buffer
- `find_markers()`: find the > markers in the mlagan alignment buffer and
  return a list of their location
- `build_interval_list()`: extract all ungapped aligned subintervals from
  a pair of aligned sequences
- `build_mlagan_ivals()`: takes a mlagan alignment file buffer as input and
  builds the ivals
- `create_NLMSA_mlagan()`: takes buffer of a mlagan alignment file,
  sequence db and NLMSA  as input and returns NLMSA
  

How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import mlagan_NLMSA``.
   You will also need to ``from pygr import cnestedlist, seqdb``.

2. Obtain the NLMSA using create_NLMSA_mlagan(buf, seqDb, al)
   function. One needs to pass blat output file object (buf), sequence
   database(seqDb) and the NLMSA object (al) to the function and the
   function returns the modified/built NLMSA.
   ``nlmsa_aln = create_NLMSA_mlagan(buf, seqDb, al)``

"""

__docformat__ = 'restructuredtext'

from pygr import cnestedlist, nlmsa_utils, seqdb

def read_mlagan(buf):
    """
    Read aligned sequences from a mlagan alignment file buffer
    """
    
    assert buf[0].startswith('>'), "This doesn't look like fasta file"

    # Extracting the sequence and names from mlagan output
    
    seqMarker_list, seqInfo_list = find_markers(buf)
    buf = buf[1:]
    seqName = []
    for seqinf in seqInfo_list:
        buf = buf.replace(seqinf,'')
        seqName.append(seqinf.split()[0])
        
    buf = buf.replace('\n','')
    seq_List = buf.split('>')
   
    return seq_List, seqName

def find_markers(buf):
    """
    Find the > markers and return a list of their locations and
    the seq info on those lines in a list
    """
    
    seqMarker_list = []
    seqInfo_list = []

    index1 = buf.find('>')
    while index1 != -1:
        seqMarker_list.append(index1)
        index2 = buf.find('\n',index1)
        seqInfo_list.append(buf[index1+1:index2])
        index1 = buf.find('>',index1+1)

    return seqMarker_list, seqInfo_list

def build_interval_list(a, b):
    """
    Hacky code to extract all ungapped aligned subintervals from a
    pair of aligned sequences.
    """
    interval_list = []

    a_start = None
    b_start = None

    a_count = b_count = 0
    for i in range(0, len(a)):
        if a[i] == '-' or b[i] == '-':
            if a_start is not None:           # want to end at i-1
                interval_list.append((a_start, a_count, b_start, b_count))

                a_start = b_start = None
        else:
            if a_start is None:
                a_start = a_count
                b_start = b_count

        if a[i] != '-':
            a_count += 1
        if b[i] != '-':
            b_count += 1

    if a_start is not None:
        interval_list.append((a_start, a_count, b_start, b_count))

    assert a_count == len(a.replace('-', ''))
    assert b_count == len(b.replace('-', ''))
        
    return interval_list

def build_mlagan_ivals(buf, seqDb):
    """
    Takes a lagan alignment file buffer as input and builds the
    ivals
    """
    seqList, seqNames = read_mlagan(buf)
    
    # a multiple alignment considered as a collection of
    # pairwise alignments so double iteration
    for i in range(0, len(seqList)):
        ivals = []
        seqs1_ival = seqDb[seqNames[i]]
        seqs1_ival_str = seqList[i]
        for j in range(i+1, len(seqList)):
            seqs2_ival = seqDb[seqNames[j]]
            seqs2_ival_str = seqList[j]
            interval_list = build_interval_list(seqs1_ival_str,
                                                    seqs2_ival_str)
            for (a, b, x, y) in interval_list:
                ival1 = (seqNames[i], a, b)
                ival2 = (seqNames[j], x, y)
                ivals.append((ival1, ival2))
        yield ivals
            
def create_NLMSA_mlagan(buf, seqDb,al):
    """
    Takes mlagan alignment file buffer as input and creates and
    returns NLMSA
    """
    for ivals in build_mlagan_ivals(buf, seqDb):
        alignedIvalsAttrs = dict(id=0, start=1, stop=2, idDest=0, 
                                 startDest=1, stopDest=2)
        cti = nlmsa_utils.CoordsToIntervals(seqDb, seqDb,
                                            alignedIvalsAttrs)
        al.add_aligned_intervals(cti(ivals))


    # build alignment
    al.build()
    return al
