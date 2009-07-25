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
        self.buf = open('ProtDNA.psl').read()
        self.buf = self.buf.replace("\r\n","\n")
                
        
    def test_parse_blat(self) :
        matches, genome_names = blat_NLMSA.parse_blat(self.buf)
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
        #test ++ orientation
        blat_aln = matches[0]
        
        qStart = getattr(blat_aln, "qStart")
        tStart = getattr(blat_aln, "tStart")
        qEnd = getattr(blat_aln, "qEnd")
        tEnd = getattr(blat_aln, "tEnd")
        
        qName = getattr(blat_aln, "qSeqName")
        tName = getattr(blat_aln, "tSeqName")
        orient = getattr(blat_aln, "orient")
        blocks = getattr(blat_aln, "blocks")
        
        
        last_ungapped = blocks[-1]
        last_ungapped_qStart = getattr(last_ungapped, "qStart")
        last_ungapped_tStart = getattr(last_ungapped, "tStart")
        last_ungapped_qEnd = getattr(last_ungapped, "qEnd")
        last_ungapped_tEnd = getattr(last_ungapped, "tEnd")
        
        #self.assertEqual(qStart, 20)
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
        self.assertEqual(last_ungapped_tEnd, 246+25)

        #test +- orientation
        blat_aln = matches[3]
        
        qStart = getattr(blat_aln, "qStart")
        tStart = getattr(blat_aln, "tStart")
        qEnd = getattr(blat_aln, "qEnd")
        tEnd = getattr(blat_aln, "tEnd")
        
        qName = getattr(blat_aln, "qSeqName")
        tName = getattr(blat_aln, "tSeqName")
        orient = getattr(blat_aln, "orient")
        blocks = getattr(blat_aln, "blocks")
        
        
        last_ungapped = blocks[-1]
        last_ungapped_qStart = getattr(last_ungapped, "qStart")
        last_ungapped_tStart = getattr(last_ungapped, "tStart")
        last_ungapped_qEnd = getattr(last_ungapped, "qEnd")
        last_ungapped_tEnd = getattr(last_ungapped, "tEnd")
        
def translator(prot_db, dna_db, of):
    """
    Iterate over the frames and write them to a single output 
    fasta format file
    """
    t_dna_db = translationDB.get_translation_db(dna_db)
    for name in dna_db.keys():
        dna_seq = dna_db[name]
        for i in range(0,6):
            cur_frame = t_dna_db.annodb[dna_seq.id + ':'+ str(i)]
            of.write(">" + name + ":" + str(i) + "\n")
            of.write(str(cur_frame) + "\n")
    for name in prot_db.keys():
        prot_seq = prot_db[name]
        of.write(">" + name + "\n")
        of.write(str(prot_db[name])+ "\n")
        
class Blat_NLMSA_test(unittest.TestCase):

    def setUp(self):
        
        self.buf = open('ProtDNA.psl').read()
        self.buf = self.buf.replace("\r\n","\n")
        
        thisdir = os.path.abspath(os.path.dirname(__file__))
        
        self.srcDB = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_prot.fa'))
        self.destDB = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_dna.fa'))
        of  = open('translatedDB.fa','w')
        translator(self.srcDB, self.destDB, of)
        of.close()
        
        self.db = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'translatedDB.fa'))
        
        matches, genome_names = blat_NLMSA.parse_blat(self.buf)
        
        alignment = cnestedlist.NLMSA('test', mode='memory', pairwiseMode=True,
                               bidirectional=False, seqDict=self.db)


        self.temp_nlmsa = blat_NLMSA.create_NLMSA_tblat(self.buf, self.srcDB, 
                                                         self.destDB, alignment)
    
   
    def test_align_manual1(self):
        """
        in this test, alignments from Protein/DNA blat output file are
        read and tested against the alignments read and built into the NLMSA
        """
        s1 = self.db['HBB0_PAGBO']
        s2 = self.db['HBB1_VAREX']
        s3 = self.db['gi|171854975|dbj|AB364477.1|:0']
    
        ival = s1[81:91]        
        temp_lst = []
        for s in self.temp_nlmsa[ival]:
            temp_lst.append(str(s))
         
        self.assertEqual(temp_lst, [str(s3[82:92])])
    
        temp_lst=[]
        ival = s2[19:33]
        for s in self.temp_nlmsa[ival]:
          temp_lst.append(str(s))
        
        self.assertEqual(temp_lst, [str(s3[20:34])])
  
        # can add additional manual tests
  
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Blat_test))
    suite.addTest(unittest.makeSuite(Blat_NLMSA_test))
    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
