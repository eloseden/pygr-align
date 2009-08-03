# protein/DNA alignment test
# Author Eden Elos

import os
import unittest
from pygr import cnestedlist, seqdb, translationDB
import blat_NLMSA

class Blat_test(unittest.TestCase):
    """
    Blat test class, to test protein-DNA alignment
    """
    def setUp(self):
        self.buf = open('data/ProtDNA.psl').read()
        self.buf = self.buf.replace("\r\n","\n")
        self.aln_type = 1  # if protein-dna-1,else-0      
        
    def test_parse_blat(self) :
        matches, genome_names = blat_NLMSA.parse_blat(self.buf, self.aln_type)
        names_set = set(['gi|171854975|dbj|AB364477.1|',
                                                'hbb1_mouse_RC',
                                                'hbb1_mouse_RC_2',
                                                'hbb1_mouse_RC_3',
                                                 'HBB0_PAGBO','HBB1_ANAMI',
                                                 'HBB1_CYGMA','HBB1_IGUIG',
                                                 'HBB1_MOUSE','HBB1_ONCMY',
                                                 'HBB1_PAGBO','HBB1_RAT',
                                                 'HBB1_SPHPU','HBB1_TAPTE',
                                                 'HBB1_TORMA','HBB1_TRICR',
                                                 'HBB1_UROHA','HBB1_VAREX',
                                                 'HBB1_XENBO','HBB1_XENLA',
                                                 'HBB1_XENTR','MYG_DIDMA',
                                                 'MYG_ELEMA','MYG_ERIEU',
                                                 'MYG_ESCGI','MYG_GALCR',
                                                 'PRCA_ANASP','PRCA_ANAVA'])
        
        self.assertEqual(set(genome_names).difference(names_set),set())
        
        blat_aln = matches[0] # Extracting the first alignment 
        
        qStart = getattr(blat_aln, "qStart")
        tStart = getattr(blat_aln, "tStart")
        qEnd = getattr(blat_aln, "qEnd")
        tEnd = getattr(blat_aln, "tEnd")
        qName = getattr(blat_aln, "qSeqName")
        tName = getattr(blat_aln, "tSeqName")
        orient = getattr(blat_aln, "orient")
        blocks = getattr(blat_aln, "blocks")
                
        last_ungapped = blocks[-1] # Extracting the last ungapped block
        last_ungapped_qStart = getattr(last_ungapped, "qStart")
        last_ungapped_tStart = getattr(last_ungapped, "tStart")
        last_ungapped_qEnd = getattr(last_ungapped, "qEnd")
        last_ungapped_tEnd = getattr(last_ungapped, "tEnd")
        
        self.assertEqual(qStart, 20)
        self.assertEqual(tStart, 63)
        self.assertEqual(qEnd, 106)
        self.assertEqual(tEnd, 321)
        self.assertEqual(qName, 'HBB0_PAGBO')
        self.assertEqual(tName, 'gi|171854975|dbj|AB364477.1|')
        self.assertEqual(orient, '++')
        self.assertEqual(len(blocks), 2)
        
        self.assertEqual(last_ungapped_qStart, 81)
        self.assertEqual(last_ungapped_tStart, 246)
        self.assertEqual(last_ungapped_qEnd, 81+25)
        self.assertEqual(last_ungapped_tEnd, 246+3*25)


class Blat_NLMSA_test(unittest.TestCase):

    def setUp(self):
        
        self.buf = open('data/ProtDNA.psl').read()
        self.buf = self.buf.replace("\r\n","\n")

        thisdir = os.path.abspath(os.path.dirname(__file__))
        def thisfile(name):
            return os.path.join(thisdir, name)
        
        
        self.srcDB = seqdb.SequenceFileDB(thisfile('data/test_prot.fa'))
        self.destDB = seqdb.SequenceFileDB(thisfile('data/test_dna.fa'))
        self.db = seqdb.SequenceFileDB(thisfile('data/translatedDB.fa'))
        
        self.aln_type = 1 #protein-dna alignmet = 1, else - 0 
        
        self.matches, genome_names = blat_NLMSA.parse_blat(self.buf, 
                                                           self.aln_type)

        self.db = translationDB.get_translation_db(self.destDB)

        alignment = cnestedlist.NLMSA('test', mode='memory', pairwiseMode=True,
                                      bidirectional=False)

        self.temp_nlmsa = blat_NLMSA.create_NLMSA_blat(self.buf, alignment,
                                                       aln_type = self.aln_type,        
                                                       srcDB = self.srcDB, 
                                                       destDB = self.db)
    
    def test_align_manual1(self):
        """
        In this test, alignments from Protein/DNA blat output file are
        read and tested against the alignments read and built into the NLMSA
        """
        s1 = self.srcDB['HBB0_PAGBO']
        s2 = self.srcDB['HBB1_VAREX']
        s3 = self.db.annodb['gi|171854975|dbj|AB364477.1|:0']
        
        ival = s1[20:30]        
        temp_lst = []
        for s in self.temp_nlmsa[ival]:
            temp_lst.append(str(s))
        
        self.assertEqual(temp_lst, [str(s3[21:31])])
    
        temp_lst=[]
        ival = s2[13:23]
        
        for s in self.temp_nlmsa[ival]:
          temp_lst.append(str(s))
        
        self.assertEqual(temp_lst, [str(s3[14:24]) ])
  
    # additional manual tests
    def test_align_auto1(self):
        """
        In this test, alignments from Protein/DNA blat output file are
        read and tested against the alignments read and built into the NLMSA
        """
        s1_name = 'HBB1_VAREX'
        s1 = self.srcDB[s1_name]
        s1_blatalns = []
        
        for blataln in self.matches:
            qName = getattr(blataln, "qSeqName")
            tName = getattr(blataln, "tSeqName")
            if qName == s1_name or tName == s1_name:
                s1_blatalns.append(blataln)
        
        s1_blck_indices = [] # Extracting the ungapped aligned blocks
        s2_group = []    # 
        
        for blataln in s1_blatalns:
            qName = getattr(blataln, "qSeqName")
            tName = getattr(blataln, "tSeqName")
            blcks = getattr(blataln, "blocks")
            for blck in blcks:
                if qName == s1_name:
                    s1_blck_indices.append((getattr(blck,"qStart"), 
                                                        getattr(blck,"qEnd")))
                    s2_group.append((tName, getattr(blck,"tStart"), 
                                                        getattr(blck,"tEnd")))
                else:
                    s1_blck_indices.append((getattr(blck,"tStart"), 
                                                        getattr(blck,"tEnd")))
                    s2_group.append((qName, getattr(blck,"qStart"), 
                                                        getattr(blck,"qEnd")))
        
        temp_lst1 = [] # list of alignments from nlmsa
        for se_indices in s1_blck_indices:
            ival = s1[se_indices[0]:se_indices[1]]
            for s in self.temp_nlmsa[ival]:
                temp_lst1.append(str(s))
                       
        temp_lst2 = [] # list of alignments from blat output
        for s2grp in s2_group:
            start = s2grp[1]/3
            end = s2grp[2]/3
            frame_no = s2grp[1]%3            
            s2 = self.db.annodb[s2grp[0]+":"+str(frame_no)]
            temp_lst2.append(str(s2[start:end]))
        
        self.assertEqual(temp_lst1, temp_lst2)
    
 
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Blat_test))
    suite.addTest(unittest.makeSuite(Blat_NLMSA_test))
    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
