import sys
import os
import subprocess
from subprocess import PIPE
#Provides classes and data structures to parse EricScript output file (TSV format)
#perform blastN
def executeBlastnOutfmt6Wrapper( ericrec, blastnpath, blastndb, numofalignments):
    blastcmd="echo -e \">"+ericrec.get_fusion_name()+"\\n"+ericrec.junctionsequence_.upper()+"\" | "+blastnpath+" -db "+blastndb+" -outfmt 6 -num_alignments "+str(numofalignments)
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

    #to check whether a fusion is a pseudogene we need to blast it against a pseudogenedb and analyse the results (outputformat6)
    def is_mapping(self, blastnpath, blastndb):
        ispseudo=False
        blasttop5rec=executeBlastnOutfmt6Wrapper(self, blastnpath, blastndb, 5)
        if len(blasttop5rec)>0:
            for brec in blasttop5rec:
                if brec.pident > 99.9999 and brec.evalue < 0.00001 and brec.length > 98:
                    ispseudo=True
                    break
        return ispseudo
    def get_mapping_database_id(self, blastnpath, blastndb):
        pseudoid=""
        blasttop5rec=executeBlastnOutfmt6Wrapper(self, blastnpath, blastndb, 5)
        if len(blasttop5rec)>0:
            for brec in blasttop5rec:
                if brec.pident > 99.9999 and brec.evalue < 0.00001 and brec.length > 98:
                    pseudoid=brec.sseqid
                    break
        return pseudoid
    def get_mapping_database_blast_records(self, blastnpath, blastndb):
        pseudorec=[]
        blasttop5rec=executeBlastnOutfmt6Wrapper(self, blastnpath, blastndb, 5)
        if len(blasttop5rec)>0:
            for brec in blasttop5rec:
                print(brec.to_string())
                if brec.pident > 98 and brec.evalue < 0.00001 and brec.mismatch <= 2 and brec.length >= 98 :
                    pseudorec.append(brec)
                    break
        return pseudorec

    #return a vector of pseudogene alignments - format 6
    #arguments: completepath of blastN, pseudogene database
    #start from here, look at the example saved in ..
    #  def blast_pseudogene_db(blastpath, psedogenedb):
    #      blastcmd="echo -e \">"+self.get_fusion_name()+"\n"+self.junctionsequence_+"\"| "+blastpath+" -db "+pseudogenedb+" -outfmt 6 -num_alignments 5"
    #      blastresult=subprocess.run(blastcmd, )

