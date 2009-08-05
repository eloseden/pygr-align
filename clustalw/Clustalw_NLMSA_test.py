# Author Eden Elos
import os
import unittest
from pygr import cnestedlist
from pygr import seqdb
import Clustalw_NLMSA

# Test for the clustalwresidues class
class ClustalwResidues_test(unittest.TestCase):

    def setUp(self):
        """
        Initialize the residues from a test file
        """
        
        thisdir = os.path.abspath(os.path.dirname(__file__))
        self.db = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_genomes_file'))
        
        lines = open(os.path.join(thisdir, 'test_clustalw_alignment.aln')).readlines()
        
        self.clustalw_resds = Clustalw_NLMSA.read_clustalw(lines)
    
    def test_num_ClustalwResidues(self):
        """
        Test the number of blocks
        """
        
        self.assertEqual(len(self.clustalw_resds),6)

    def test_num_seq_in_all(self):
        """
        Test the number of seqs in each block
        """
        
        for clus_res in self.clustalw_resds:
            self.assertEqual(clus_res.get_no_seqs(),4)

    def test_seqs_names_in_all(self):
        """
        Test the names of the sequences in each block
        hast to be the same (for each sequence) in the different blocks
        """
        
        for clus_res in self.clustalw_resds:
            self.assertEqual(clus_res.get_names(),['query','P15522',
                                                   'AAB85326.1','NP009141'])
     
    def test_start_indices(self):
        """
        Test the start indices of the sequences in each block
        has to be the same (for each sequence) in the different blocks
        """
        
        start_indices = [[0,0,0,0],[41,13,60,45],[55,27,110,69],
                       [90,60,147,129],[121,91,182,189],[141,111,201,249]]
        for i in  range(0, len(self.clustalw_resds)):
            self.assertEqual(self.clustalw_resds[i].get_start_indices(),
                             start_indices[i])

    def test_end_indices(self):
        """
        Test the end indices of the sequences in each block
        must correpond to values read from the test file
        """
        
        end_indices = [[40,12,59,44],[54,26,109,68],[89,59,146,128],
                     [120,90,181,188],[141,111,201,248],[141,111,201,279]]
        for i in  range(0, len(self.clustalw_resds)):
            self.assertEqual(self.clustalw_resds[i].get_end_indices(),
                             end_indices[i])
            
    
    def test_ungapped_count(self):
        """
        A sample test...just from the last block
        """
        
        ungapped_count_last = [0,0,0,31]
        self.assertEqual(self.clustalw_resds[5].ungapped_count(),
                         ungapped_count_last)
    
# Test for the clustalwresidues class
class Clustalw_NLMSA_test(unittest.TestCase):

    
    def setUp(self):
        """
        Initilize/create the NLMSA using : .aln alignment file,
        the file containing the genomes created with in memory mode,
        read in the sequences from the genomes file
        """

        thisdir = os.path.abspath(os.path.dirname(__file__))
        
        def thisfile(name):
            return os.path.join(thisdir, name)
        
        self.db = seqdb.SequenceFileDB(thisfile('test_genomes_file'))
        
        buf = open(thisfile('test_clustalw_alignment.aln'), "r").read()
        
        clustal_res_list = Clustalw_NLMSA.read_clustalw(buf.split("\n"))
        
       

        alignment = cnestedlist.NLMSA('test', mode='memory', seqDict=self.db,
                                      use_virtual_lpo=True)
    
        self.temp_nlmsa = Clustalw_NLMSA.create_NLMSA_clustalw(buf, self.db,
                                                               alignment)

        
    def test_align_manual1(self):
        """
        In this test, alignments from the clustalw alignment
        file are read and tested against the
        alignemnts built into the NLMSA
        perhaps, a better systematic testing can be designed
        """
        
        s1 = self.db['query']
        temp_lst = []
        for s2 in self.temp_nlmsa[s1[:10]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst,['GSFRVLKSRT','RRRHMPLRLA'])

        s2 = self.db['NP009141']
        temp_lst = []
        
        for s2 in self.temp_nlmsa[s2[-10:]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst, [])

        # can add additional manual tests
        
            
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ClustalwResidues_test))
    suite.addTest(unittest.makeSuite(Clustalw_NLMSA_test))
    # suite.addTest(unittest.makeSuite(Clustalw_MAF_NLMSA_test))
    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
