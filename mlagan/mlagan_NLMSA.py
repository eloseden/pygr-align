# Author Eden Elos
# July,01
# create NLMSA from mlagan alignment file

# ! /usr/bin/env python2.5

import os
import glob
from pygr import cnestedlist, seqdb

def read_mlagan(buf):
    """
    Read aligned sequences from a lagan alignment file buffer
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
    # find the > markers and return a list of their locations and
    # also return the seq info on those lines in a list
    
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


def create_NLMSA_mlagan(buf, seqDb, al):
    """
        takes buffer of a mlagan alignment file as input and creates and
        returns NLMSA
    """
    seqList, seqNames = read_mlagan(buf)

    # feed the alignment
    al += seqDb[seqNames[0]]

    # a multiple alignment considered as a collection of
    # pairwise alignments so double iteration
    for i in range(0, len(seqList)):
        genome1_ival = seqDb[seqNames[i]]
        genome1_ival_str = seqList[i]
        for j in range(i+1, len(seqList)):
            genome2_ival = seqDb[seqNames[j]]
            genome2_ival_str = seqList[j]
            interval_list = build_interval_list(genome1_ival_str,
                                                    genome2_ival_str)
            for (a, b, x, y) in interval_list:
                ival1 = genome1_ival[a:b]
                ival2 = genome2_ival[x:y]
                al[ival1] += ival2

    al.build()
    return al
    

