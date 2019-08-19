import sys
import os
import subprocess
from subprocess import PIPE
#Provides classes and data structures to parse EricScript output file (TSV format)
#perform blastN

class LocalBlastN():
    def __init__(self, blastnpath, blastdb):
        self.blastnpath_=blastnpath
        self.blastndb_=blastdb
    def executeBlastnOutfmt6Wrapper(self, sequence , numofalignments):
        blastcmd="echo -e \">"+sequence.get_name()+"\\n"+sequence.get_sequence().upper()+"\" | "+self.blastnpath_+" -db "+self.blastndb_+" -outfmt 6 -num_alignments "+str(numofalignments)
        popenobj=subprocess.Popen(blastcmd, shell=True, stdout=PIPE, stderr=PIPE)
        blastnformat6recs=[]
        for bl in popenobj.stdout:
            blastnformat6recs.append(BlastOutputFormat6Record(bl.decode("utf-8")))
        return blastnformat6recs
#Single record for blast output format 6 - a 12 column record
class BlastOutputFormat6Record():
    def __init__(self, rec):
        self.srec=rec.strip().split("\t")
        self.qseqid=self.srec[0]
        self.sseqid=self.srec[1]
        self.pident=float(self.srec[2])
        self.length=int(self.srec[3])
        self.mismatch=int(self.srec[4])
        self.gapopen=int(self.srec[5])
        self.qstart=int(self.srec[6])
        self.qend=int(self.srec[7])
        self.sstart=int(self.srec[8])
        self.send=int(self.srec[9])
        self.evalue=float(self.srec[10])
        self.bitscore=float(self.srec[11])
    def to_string(self):
        return "\t".join(self.srec)
    def is_identical(self, thresh=99.99):
        return self.pident >= thresh
    def is_significant(self, thresh=0.0001):
        return self.evalue < thresh
