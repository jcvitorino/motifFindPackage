import InOut
import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class TestInOut(unittest.TestCase):
        
    def test_next_vertex(self):
        s = Seq("AAAA",IUPAC.unambiguous_dna)
        i = 0
        L = 4
    
        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAA", str(s))
        self.assertEqual(1, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAA", str(s))
        self.assertEqual(2, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAA", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAC", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAG", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAAT", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AACT", str(s))
        self.assertEqual(2, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AACA", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AACC", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AACG", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AACT", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAGT", str(s))
        self.assertEqual(2, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAGA", str(s))
        self.assertEqual(3, i)

        (s,i) = InOut.next_vertex(s,i,L,'T')
        self.assertEqual("AAGC", str(s))
        self.assertEqual(3, i)

    def test_bypass(self):
        s = Seq("AAAA",IUPAC.unambiguous_dna)
        i = 0

        (s,i) = InOut.bypass(s,i,'T')
        self.assertEqual("CAAA", str(s))
        self.assertEqual(0, i)

        s = Seq("CTAA",IUPAC.unambiguous_dna)
        (s,i) = InOut.bypass(s,1,'T')
        self.assertEqual("GTAA", str(s))
        self.assertEqual(0, i)

        (s,i) = InOut.bypass(s,2,'T')
        self.assertEqual("GTCA", str(s))
        self.assertEqual(2, i)

    def test_total_distance(self):
        s = Seq("A",IUPAC.unambiguous_dna)
        i = 0
        seq = InOut.read('seqT.fasta')
        r = InOut.total_distance(s,seq,1)
        self.assertEqual(0, r)

        s = Seq("AA",IUPAC.unambiguous_dna)
        i = 0
        seq = InOut.read('seqT.fasta')
        r = InOut.total_distance(s,seq,2)
        self.assertEqual(1, r)


        s = Seq("TT",IUPAC.unambiguous_dna)
        i = 0
        seq = InOut.read('seqT.fasta')
        r = InOut.total_distance(s,seq,2)
        self.assertEqual(0, r)

        s = Seq("TAGTT",IUPAC.unambiguous_dna)
        i = 0
        seq = InOut.read('seqT.fasta')
        r = InOut.total_distance(s,seq,5)
        self.assertEqual(0, r)

        s = Seq("TC",IUPAC.unambiguous_dna)
        i = 0
        seq = InOut.read('seqT.fasta')
        r = InOut.total_distance(s,seq,2)
        self.assertEqual(1, r)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
