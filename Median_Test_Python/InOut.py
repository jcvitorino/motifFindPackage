from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
from time import clock, time

def read (path):
    try:
        records = list(SeqIO.parse(path, "fasta"))
    except IOError, e:
        print 'could not open file:', e
        exit()
    return records

def hamming_distance(s1, s2):
    if len(s1) == len(s2):
        assert len(s1) == len(s2)
        return sum([ch1 != ch2 for ch1, ch2 in zip(s1, s2)])
    else:
        print '++++++++ string not equal-----------------'
        return 200000

def total_distance(prefix, seq, sizeP):
    result = list()
    for i in seq:
        dist_1 = 2000000
        #print '((len(i.seq) - sizeP )+1)'
        #print len(i.seq)
        for j in xrange(0,(len(i.seq) - sizeP )+1):
        #    print (i.seq[0:18])
            dist_2 = hamming_distance(prefix,i.seq[j:(j+((sizeP+1) - 1))])
            if dist_1 > dist_2:
                dist_1 = dist_2
        result.append(dist_1)
    return sum(result)
    
def bypass(a,i,k):
    j = i
    while j != -1:
        if a[j] != k:
            a = a.tomutable()
            if a[j] == 'A':
                a[j] = 'C'
                a = a.toseq()
                return (a,j)
            if a[j] == 'C':
                a[j] = 'G'
                a = a.toseq()
                return (a,j)
            if a[j] == 'G':
                a[j] = 'T'
                a = a.toseq()
                return (a,j)
        j = j - 1
    return (a,-1)

def next_vertex(a,i,L,k):
    j = L - 1
    a = a.tomutable()
    if i < L - 1:
        a[i+1] = 'A'
        a = a.toseq()
        return (a,(i+1))
    else:
        while j != -1:
            if a[j] != k:
                if a[j] == 'A':
                    a[j] = 'C'
                    a = a.toseq()
                    return (a,j)
                if a[j] == 'C':
                    a[j] = 'G'
                    a = a.toseq()
                    return (a,j)
                if a[j] == 'G':
                    a[j] = 'T'
                    a = a.toseq()
                    return (a,j)
            j = j - 1
        a = a.toseq()
    return (a,-1)

def median_string(seq,L,mis):
    # number of samples
    t = len(seq)
    # start position
    s = Seq("",IUPAC.unambiguous_dna)
    for i in range(L):
        s = s + Seq("A",IUPAC.unambiguous_dna)
        
    i = 0
    best_distance = 100
    best_word = Seq("",IUPAC.unambiguous_dna)
    while i > -1:
        if i < L - 1:
            prefix = s[0:i]
            optimistic_distance = total_distance(prefix, seq, len(prefix))
        #    print optimistic_distance
            if optimistic_distance > best_distance:
                (s,i) = bypass(s,i,'T')
            else:
                (s,i) = next_vertex(s,i,L,'T')

        else:
            total_dist = total_distance(s, seq, len(s))
            if total_dist <= best_distance:
                best_distance = total_dist
                best_word = s
                if best_distance <= mis:
                    file_seq = open('result.txt','a')
                    file_seq.write(str(best_word))
                    file_seq.write('\n')
                    file_seq.close()
            (s,i) = next_vertex(s,i,L,'T')
    return best_word
                
def main():
    params = sys.argv[1:]
    try:
        # open and read .fasta file
        sequences = read(params[0])
        # length L of the motif
        L = params[1]
        # mismatches number
        mismatches = params[2]
    except IndexError:
        print 'wrong number of parameters'
        print 'Right way: \'program.py your_file.fasta motif_length mismatches_number\' '
        exit()
    ret = median_string(sequences,int(L),int(mismatches))
    print ret
    
if __name__ == '__main__':
    start = time()
    main()
    r = time() - start
    print r
