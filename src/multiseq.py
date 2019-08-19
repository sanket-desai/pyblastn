#!/usr/bin/python
import sys
import os.path
import re
from seq import Sequence

class MultiSequence(object):

	def __init__(self, fn=None):
		self.seqs = []
		if fn is None:
			self.seqs=[]
		elif fn != None and os.path.exists(fn):
			fi = open(fn,'r')
			rec=[]
			b=False
			for ln in fi:
				fl = ln.rstrip()
				if fl.startswith(">") and b:
					#sq = self.createSeq(rec)
					sn,ss="",""
					for i in rec:
						if i.startswith(">"):
							sn=i[1:]
						else:
							ss=ss+i.replace(" ","")
					sq = Sequence(sn,ss)						
					del rec[:]
					self.seqs.append(sq)
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
		return len(self.seqs)
	def __getitem__(self,i):
		return self.seqs[i]
	def __setitem__(self,i,j):
		self.seqs[i] = j
	def average_sequence_length(self):
		totlen=0
		for sq in self.seqs:
			totlen+=len(sq)
		return float(totlen/len(self.seqs))
	def append(self, sqn):
		if type(sqn) is Sequence:
			self.seqs.append(sqn)
		else:
			TypeError("Sequence object to be appended. Object type not 'Sequence'")
	def write(self, fo): # file handle with append permissions
		for s in self.seqs:
			s.write(fo)
		fo.close()
