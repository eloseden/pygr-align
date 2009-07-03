# Author Eden Elos
import os
import unittest
from pygr import cnestedlist
from pygr import seqdb
import blastz_NLMSA

class Blastz_test(unittest.TestCase):
    """
    blastz test class to test the various functions 
    """
    
    def setUp(self):
        self.buf = open('output').read()
        self.buf = self.buf.replace("\r\n","\n")
        
    def test_find_lav_marker(self):
        self.assertEqual(len(blastz_NLMSA.find_lavmarkers(self.buf)), 2)
        
    def test_construct_coord(self):
        lavmarkers_list = blastz_NLMSA.find_lavmarkers(self.buf)
        self.assertEqual(blastz_NLMSA.construct_coord(lavmarkers_list),
                         [(218, 490)])
        
    def test_get_names(self):
        lav_marker_list = blastz_NLMSA.find_lavmarkers(self.buf)
        lav_coords = blastz_NLMSA.construct_coord(lav_marker_list)
        for coord in lav_coords:
           start_index = int(coord[0])
           end_index = int(coord[1])
           lav_block = self.buf[start_index:end_index]
           records = lav_block.split('\n}\n')
        self.assertEqual(blastz_NLMSA.get_names(records),
                         ['testgenome1', 'testgenome2'])
        
    def test_parse_blastz(self) :
        matches, genome_names = blastz_NLMSA.parse_blastz(self.buf)
        self.assertEqual(set(genome_names), set(['testgenome1', 'testgenome2']))
        
        blastz_aln = matches[-1]
        score = getattr(blastz_aln, "score")
        start_top = getattr(blastz_aln, "start_top")
        start_bot = getattr(blastz_aln, "start_bot")
        end_top = getattr(blastz_aln, "end_top")
        end_bot = getattr(blastz_aln, "end_bot")
        genome_name1 = getattr(blastz_aln, "genome_name1")
        genome_name2 = getattr(blastz_aln, "genome_name2")
        orient = getattr(blastz_aln, "orient")
        blocks = getattr(blastz_aln, "blocks")
        
        last_ungapped = blocks[-1]
        last_ungapped_start_top = getattr(last_ungapped, "start_top")
        last_ungapped_start_bot = getattr(last_ungapped, "start_bot")
        last_ungapped_end_top = getattr(last_ungapped, "end_top")
        last_ungapped_end_bot = getattr(last_ungapped, "end_bot")
        last_ungapped_ident = getattr(last_ungapped, "ident")

        self.assertEqual(score, 74457)
        self.assertEqual(start_top, 40)
        self.assertEqual(start_bot, 41)
        self.assertEqual(end_top, 1120)
        self.assertEqual(end_bot, 1120)
        self.assertEqual(genome_name1, 'testgenome1')
        self.assertEqual(genome_name2, 'testgenome2')
        self.assertEqual(orient, 1)
        self.assertEqual(len(blocks), 4)
        
        self.assertEqual(last_ungapped_start_top, 302)
        self.assertEqual(last_ungapped_start_bot, 302)
        self.assertEqual(last_ungapped_end_top, 1120)
        self.assertEqual(last_ungapped_end_bot, 1120)
        self.assertEqual(last_ungapped_ident, 84)


class Blastz_NLMSA_test(unittest.TestCase):

    def setUp(self):

        self.buf = open('output').read()
        self.buf = self.buf.replace("\r\n","\n")

        thisdir = os.path.abspath(os.path.dirname(__file__))
        self.db = seqdb.SequenceFileDB(os.path.join(thisdir,
                                                    'test_genomes.fna'))
        matches, genome_names = blastz_NLMSA.parse_blastz(self.buf)
        
        alignment = cnestedlist.NLMSA('test', mode='memory', seqDict=self.db,
                                      use_virtual_lpo=True)
        alignment += self.db[genome_names[0]]

        self.temp_nlmsa = blastz_NLMSA.create_NLMSA_blastz(self.buf,
                                                   self.db, alignment)
    
        
        

    def test_align_manual1(self):
        """
        In this test, manual alignments from the blastz alignment file are
        read and tested against the alignments read and built into the NLMSA
        """
        
        s1=self.db['testgenome1']
        temp_lst=[]
        
        for s2 in self.temp_nlmsa[s1[40:50]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst,['TGGTTGAAAA'])

        s2=self.db['testgenome2']
        temp_lst=[]
        
        for s2 in self.temp_nlmsa[s2[:10]]:
            temp_lst.append(str(s2))
        self.assertEqual(temp_lst,[])

        # can add additional manual tests
    
        

    
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Blastz_test))
    suite.addTest(unittest.makeSuite(Blastz_NLMSA_test))
    return suite


if __name__=="__main__":
    # unittest.main()
    unittest.TextTestRunner(verbosity=2).run(suite())
