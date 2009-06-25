# Author Eden Elos
# June 17,08

# create NLMSA from Blastz alignment file

# ! /usr/bin/env python2.5

import os
import glob
from pygr import cnestedlist, seqdb


#
# BlastzLocalAlignment
#

class BlastzLocalAlignment:
    """
    A blastz gapped local alignment, consisting of multiple
    ungapped blocks.
    """
    def __init__(self, score, start_top, end_top, start_bot, end_bot,
                 genome_name1, genome_name2, orient, blocks):
        self.score = score
        self.start_top = start_top
        self.end_top = end_top
        self.start_bot = start_bot
        self.end_bot = end_bot

        self.genome_name1 = genome_name1
        self.genome_name2 = genome_name2
        
        assert abs(orient) == 1
        self.orient = orient

        self.blocks = blocks
    

#
# BlastzUngappedBlock
#

class BlastzUngappedBlock:
    """
    A single ungapped block in a blastz alignment.
    """
    def __init__(self, start_top, end_top, start_bot, end_bot, ident):
        self.start_top = start_top
        self.end_top = end_top
        self.start_bot = start_bot
        self.end_bot = end_bot

        self.ident = ident

    def __len__(self):
        return self.end_top - self.start_top

    def convert_to_text(self, seq1, seq2):
        """
        Return the ungapped alignment, as sequence.
        """
        top = seq1[self.start_top:self.end_top]
        bot = seq2[self.start_bot:self.end_bot]

        return (top, bot)



def find_lavmarkers(buf):
    """
    """
    lav_marker_list = []
    next_block = 1

    while True:
        lav_marker = buf.find('#:lav',next_block)
        next_block = lav_marker + 1

        if lav_marker != -1:
            lav_marker_list.append(lav_marker)
        elif len(lav_marker_list) >= 1:
            lav_marker_list.append(len(buf)-1)
            break
        else:
            break
    return lav_marker_list

def construct_coord(lav_marker_list):
    """
    """
    coords = []
    for i in range(0,len(lav_marker_list) - 1):
        if i != len(lav_marker_list) - 2:
            coords.append((lav_marker_list[i], lav_marker_list[i+1]-1))
        else:
            coords.append((lav_marker_list[i], lav_marker_list[i+1]))
            
    return coords
        

def get_orient(records):
    """
    This is the orientation to be obtained from each lav block
    The question is how do you make use of the orientation information in
    NLMSA (if at all) ?
    """
    orient = 1
    return orient

def get_names(records):
    """
    """
    names = []
    tempstr = records[1]

    pos = tempstr.find('h')
    index1 = tempstr.find('>',pos)
    index2 = tempstr.find('"',index1)

    index3 = tempstr.find('>',index2)
    index4 = tempstr.find('"',index3)

    if index1 != -1 or index2 != -1:
        names.append(tempstr[index1+1:index2])
    if index3 != -1 or index4 != -1:
        names.append(tempstr[index3+1:index4])

    return names

def parse_blastz(buf):
    """
    
    """
    assert buf[0:5] == '#:lav'," This does not look like a blastz file"
    lav_marker_list=find_lavmarkers(buf)

    if lav_marker_list:
        lav_coords=construct_coord(lav_marker_list)
    else:
        return [], []

    genome_names = set()
    matches = []
    for coord in lav_coords:
        start_index = int(coord[0])
        end_index = int(coord[1])
        lav_block = buf[start_index:end_index]
        records = lav_block.split('\r\n}\r\n')
        orient = get_orient(records)
        names = get_names(records)
        genome_names = genome_names.union(set(names))
        matches.extend(_parse_blastz_record_block(records, orient, names[0],
                                                  names[1]))

    return matches, list(genome_names)

def _parse_blastz_record_block(records, orient, genome_name1, genome_name2):
    """
    Run through each alignment block, parsing out the score and
    overall begin/end coords, as well as the individual ungapped
    blocks.
    """
    matches = []
    for record in records:
        # convert by linebreaks
        lines = record.split('\r\n')

        # get rid of comments
        lines = [ i for i in lines if len(i) and i[0] != '#' ]

        # get rid of now-empty records (should only be last line of file)
        if not lines:
            continue

        # double check format!
        assert lines[0][2] == '{'
        
        record_type = lines[0][0]

        if record_type == 'a':
            lines = lines[1:]
            matches.append(_parse_record(lines, orient, genome_name1,
                                         genome_name2))
        else:
            continue

    return matches


def _parse_record(record, orient, genome_name1, genome_name2):
    """
    Parse individual lines in an "a {" record block, and return a
    BlastzLocalAlignment.
    """
    record = [ i.strip().split() for i in record ]

    blocks = []
    score = begin_coords = end_coords = None
    for i in record:
        if i[0] == 's':                 # score line
            assert score is None
            score = int(i[1])
        elif i[0] == 'b':               # begin coords line
            assert begin_coords is None
            begin_coords = map(int, i[1:])
        elif i[0] == 'e':               # end coords line
            assert end_coords is None
            end_coords = map(int, i[1:])
        elif i[0] == 'l':               # ungapped block line
            blocks.append(map(int, i[1:]))
    
    start_top, start_bot = begin_coords
    end_top, end_bot = end_coords

    blocks = [ BlastzUngappedBlock(a - 1, c, b - 1, d, e) \
               for (a, b, c, d, e) in blocks ]
    
    return BlastzLocalAlignment(score,
                                start_top - 1, end_top,
                                start_bot - 1, end_bot,genome_name1,genome_name2,
                                orient, blocks)

        
def create_NLMSA_blastz(buf, seqDb, al):
    """
        takes blastz file buffer as input and creates and returns NLMSA
    """
    blastzaln_list, genome_names = parse_blastz(buf)
    
    #feed the genomes
    for i in range(0,len(genome_names)):
        al += seqDb[genome_names[i]]


    for blz_al in blastzaln_list:

        start1 = getattr(blz_al, "start_top")
        stop1 = getattr(blz_al, "end_top")
        genome_name1 = getattr(blz_al, "genome_name1")

        start2 = getattr(blz_al, "start_bot")
        stop2 = getattr(blz_al, "end_bot")       
        genome_name2= getattr(blz_al, "genome_name2")
            
        block = getattr(blz_al, "blocks")
        for ungapped in block:
            
            a = getattr(ungapped, "start_top")
            b = getattr(ungapped, "end_top")
            
            x = getattr(ungapped, "start_bot")
            y = getattr(ungapped, "end_bot")

            ival1 = seqDb[genome_name1][a:b]
            ival2 = seqDb[genome_name2][x:y]

            
            al[ival1] += ival2

    al.build()
    return al
            
        


