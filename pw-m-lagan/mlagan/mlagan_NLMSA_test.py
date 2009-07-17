# Author Eden Elos
import os
import unittest
from pygr import cnestedlist
from pygr import seqdb
import mlagan_NLMSA



class mlagan_NLMSA_test(unittest.TestCase):

    
    def setUp(self):
        """
        initilize/create the NLMSA using mlagan alignment file
        the file containing the genomes created with in memory mode,
        read in the sequences from the genomes file
        """

        self.buf = open('output').read()
        self.buf = self.buf.replace("\r\n","\n")
        
        thisdir = os.path.abspath(os.path.dirname(__file__))
        self.db = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_genomes.fna'))
        
        
        seqList, seq_names = mlagan_NLMSA.read_mlagan(self.buf)
        
        alignment = cnestedlist.NLMSA('test', mode='memory', seqDict=self.db,
                                      use_virtual_lpo=True)


        
        

        
        self.temp_nlmsa = mlagan_NLMSA.create_NLMSA_mlagan(self.buf, self.db,
                                                               alignment)

        
    def test_align_manual1(self):
        """
        in this test, alignments from the mlagan alignment
        file are read and tested against the
        alignemnts built into the NLMSA
        """

        s1 = self.db['testgenome1']
        temp_lst = []
        
        for s2 in self.temp_nlmsa[s1[71:86]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst,['GCTTTTCATTCTGAC', 'GCTTTTCATTCTGAC'])

       
        temp_lst = []
        
        for s2 in self.temp_nlmsa[s1[0:10]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst, [])

        # can add additional manual tests
        
            
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(mlagan_NLMSA_test))

    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
