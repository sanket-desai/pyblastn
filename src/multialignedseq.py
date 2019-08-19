import sys
import os.path
from alignedseq import AlignedSequence

class MultiAlignedSequence(object):
	def __init__(self,fn=None):
		self.aseqs=[]
		if fn is None:
			self.aseqs=[]
		elif fn != None and os.path.exists(fn):
			fi = open(fn,'r')
			rec=[]
			b=False
			for ln in fi:
				fl = ln.rstrip()
				if fl.startswith(">") and b:
					sn,ss="",""
					for i in rec:
						if i.startswith(">"):
							sn=i[1:]
						else:
							ss=ss+i.replace(" ","")
					asq = AlignedSequence(sn,ss)
					del rec[:]
					self.aseqs.append(asq)
					rec.append(fl)
				elif fl.startswith(">") and not b and not rec:
					b=True
					rec.append(fl)
				elif not fl.startswith(">"):
					rec.append(fl)
			sn,ss="",""
			for i in rec:
				if i.startswith(">"):
					sn=i[1:]
				else:
					ss=ss+i.replace(" ","")
			asq = AlignedSequence(sn,ss)
			del rec[:]
			self.aseqs.append(asq)
		else:
			print( "File",fn,"does not exist!!")
			sys.exit()
	def __len__(self):
		return len(self.aseqs[0])
	def __getitem__(self,i):
		return self.aseqs[i]
	def __setitem__(self,i,j):# i is the index, j is the AlignedSeq object
		self.aseqs[i]=j
	def column(self, i): #Given index return the column in the msa
		col = list()
		for as_ in self.aseqs:
			col.append(as_[i])
		return col
	def number_of_sequences(self):
		return len(self.aseqs)
