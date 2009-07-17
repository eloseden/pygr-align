# Author Eden Elos
# July,01
# create NLMSA from lagan alignment file

# ! /usr/bin/env python2.5



"""
LAGAN_NLMSA MODULE
===================
A module that parses lagan output from lagan alignment program, and
builds pygr NLMSAs with it. The module does not define any class.



Functions:
- `read_lagan()`: read aligned sequences from a lagan alignment file
  buffer
- `build_interval_list()`: extract all ungapped aligned subintervals from
  a pair of aligned sequences
- `build_lagan_ivals()`: takes a lagan alignment file buffer as input and
  builds the ivals
- `create_NLMSA_lagan()`: takes buffer of a lagan alignment file,
  sequence db and NLMSA  as input and returns NLMSA
  
How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import lagan_NLMSA``.
   You will also need to ``from pygr import cnestedlist, seqdb``.

2. Obtain the NLMSA using create_NLMSA_lagan(buf, seqDb, al)
   function. One needs to pass blat output file object (buf), sequence
   database(seqDb) and the NLMSA object (al) to the function and the
   function returns the modified/built NLMSA.
   ``nlmsa_aln = create_NLMSA_lagan(buf, seqDb, al)``

"""

__docformat__ = 'restructuredtext'


from pygr import cnestedlist, seqdb

def read_lagan(buf):
    """
    Read aligned sequences from a lagan alignment file buffer
    """
    
    assert buf[0].startswith('>'), "This doesn't look like fasta file"

    # Extracting the sequence and names from lagan output
    index1 = buf.find('>')
    index2 = buf.find('\n')
    index3 = buf.find('>', index1+1)
    index4 = buf.find('\n', index3)
    seqInfo1 = buf[index1+1:index2]
    seqInfo2 = buf[index3+1:index4]
    buf = buf[1:]
    buf = buf.replace(seqInfo1,'')
    buf = buf.replace(seqInfo2,'')
    buf = buf.replace('\n','')
    seq_List = buf.split('>')
    seqName1 = seqInfo1.split()[0]
    seqName2 = seqInfo2.split()[0]
    
    return seq_List, (seqName1, seqName2)


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


def build_lagan_ivals(buf, seqDb):
    """
    Takes a lagan alignment file buffer and sequence db as input and
    builds the ivals
    """
    
    seqList, seqNames = read_lagan(buf)

    # Extract ungapped intervals
    interval_list = build_interval_list(seqList[0], seqList[1])
    ivals = []
    for (a, b, x, y) in interval_list:
             ival1 = seqDb[seqNames[0]][a:b]
             ival2 = seqDb[seqNames[1]][x:y]
             ivals.append((ival1,ival2))

    yield ivals

def create_NLMSA_lagan(buf, seqDb,al):
    """
    Takes a lagan alignment file buffer as input and creates and
    returns NLMSA
    """
    for ivals in build_lagan_ivals(buf, seqDb):
       al.add_aligned_intervals(ivals)

    # build alignment
    al.build()
    return al


    
