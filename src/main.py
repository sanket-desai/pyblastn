from localblastn import *
from seq import *
import pysam
import sys

def main():
    #blastnpath=
    #blastdbpath=
    if len(sys.argv) < 4:
        print("Please check the arguments!")
        print("Usage: python main.py <samfile> <blastnpath> <blastndb> <outputfile>")
        sys.exit()
    else:
        blastnpath=sys.argv[2]
        blastdbpath=sys.argv[3]
        localblastn=LocalBlastN(blastnpath, blastdbpath)
        samfile=pysam.AlignmentFile(sys.argv[1],'r')
        ofi=open(sys.argv[4],'w')
        for r in samfile:
            sq = Sequence(r.query_name+"#"+r.reference_name, r.seq)
            blastoutrecs=localblastn.executeBlastnOutfmt6Wrapper(sq, 1)
            if len(blastoutrecs)>0:
                blast6out=BlastOutputFormat6Record(blastoutrecs[0])
                print(blast6out.to_string())
                if blast6out.is_identical(10) and blast6out.is_significant():
                    #print(blast6out.to_string())
                    ofi.write("%s\n" %(blast6out.to_string()) )
            else:
                print("No hits found for : %s" %(sq.get_name()))
        ofi.close()
if __name__=="__main__":
    main()
