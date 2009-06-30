# Author Eden Elos
import os
import unittest
from pygr import cnestedlist
from pygr import seqdb
import blat_NLMSA

class Blat_test(unittest.TestCase):
    """
    """
    def setUp(self):
        self.buf = open('output.psl').read()
        self.buf = self.buf.replace("\r\n","\n")
        

    def test_calculate_end(self):
        ends = blat_NLMSA.calculate_end([10, 20], [12, 13])
        self.assertEqual(ends, [22, 33])
        
        
    def test_parse_blat(self) :
        matches, genome_names = blat_NLMSA.parse_blat(self.buf)
        self.assertEqual(set(genome_names), set(['testgenome1', 'testgenome2',
                                                 'testgenome3', 'testgenome4']))
        
        blat_aln = matches[-1]
        
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
        


        self.assertEqual(qStart, 64)
        self.assertEqual(tStart, 64)
        self.assertEqual(qEnd, 420)
        self.assertEqual(tEnd, 350)
        self.assertEqual(qName, 'testgenome4')
        self.assertEqual(tName, 'testgenome1')
        self.assertEqual(orient, 1)
        self.assertEqual(len(blocks), 2)
        
        self.assertEqual(last_ungapped_qStart, 351)
        self.assertEqual(last_ungapped_tStart, 281)
        self.assertEqual(last_ungapped_qEnd, 351+69)
        self.assertEqual(last_ungapped_tEnd, 281+69)




class Blat_NLMSA_test(unittest.TestCase):

    def setUp(self):
        self.buf = open('output.psl').read()
        self.buf = self.buf.replace("\r\n","\n")
        
        thisdir = os.path.abspath(os.path.dirname(__file__))
        self.db = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_genomes.fna'))
        matches, genome_names = blat_NLMSA.parse_blat(self.buf)
        
        alignment = cnestedlist.NLMSA('test', mode='memory', seqDict=self.db,
                                      use_virtual_lpo=True)
        alignment += self.db[genome_names[0]]



        self.temp_nlmsa = blat_NLMSA.create_NLMSA_blat(self.buf,
                                                   self.db, alignment)
    
        
        

    def test_align_manual1(self):
        # in this test, manual alignments from the blat alignment file are
        # read and tested against the alignments read and built into the NLMSA
        
        s1 = self.db['testgenome1']
        s2 = self.db['testgenome2']
        s3 = self.db['testgenome3']
        s4 = self.db['testgenome4']
        
        temp_lst=[]
        
        for s in self.temp_nlmsa[s1[281:300]]:
            temp_lst.append(str(s))
        self.assertEqual(temp_lst,[str(s2[281:300]),str(s3[351:370]),
                                   str(s4[351:370])])


        temp_lst=[]
        
        for s in self.temp_nlmsa[s1[:10]]:
            temp_lst.append(str(s))
        self.assertEqual(temp_lst,[])

        # can add additional manual tests
    
        

    
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Blat_test))
    suite.addTest(unittest.makeSuite(Blat_NLMSA_test))
    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
