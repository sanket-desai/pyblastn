#!/usr/bin/python
import sys
import os.path
import re

class ResidueError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class Sequence(object):
	def __init__(self, fn=None, ssq=None):
		self.seq_=""
		self.name_=""
		if fn == None and ssq == None:
			self.name_="Unknown"
		elif fn != None and ssq == None and os.path.exists(fn):
			fi = open(fn,'r')
			fl = fi.readline().rstrip()
			if fl.startswith(">"):
				self.name_ = fl[1:]
				for ln in fi:
					self.seq_ = self.seq_ + ln.rstrip()
				self.seq_ = self.seq_.replace(" ","")
			else:
				print (fn,"Not Fasta format file!")
				sys.exit()
		elif fn != None and ssq != None:
			self.name_=fn
			self.seq_=ssq
			self.seq_ = self.seq_.replace(" ","")
		elif ssq == None and fn != None and not os.path.exists(fn):
			print ("File",fn,"does not exist!!")
			sys.exit()
		if not self.is_dna():
			self.is_protein()
			
	def is_dna(self):
		return not re.search(r"[^ATGC]",self.seq_)
	def is_protein(self):
		if self.is_dna():
			return False
		else:
			for i in self.seq_:
				if i not in ['U','G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','K','R','H','D','E','X','Z','B']: #'X' a feature where the identity of the amino acid is unknown (an X is shown at this position in the sequence) and the only information concerning the modification is that the N-terminus is blocked: P80979 (Blocked amino end (Xaa))
#'Z' - Note: Pyro-Glu is often indicated in papers as ‘pGlu’ and sometimes, in one-letter code as “U”, although this is now used for selenocysteine. In figures of publications, it may be cited as Z, pQ or E
					raise ResidueError("Residue '%s' cannot be identified as either a nucleotide or amino acid for sequence %s."%(i, self.name_))
			return True		
	def write(self): #File handle, write permissions
		n = ">"+self.name_
		print(n)
		c=0
		n=self.seq_
		while c<len(n):
			print(n[c:c+60])
			c=c+60
	def write(self,fo): #File handle, write permissions
		n = ">"+self.name_
		fo.write(n)
		fo.write("\n")
		c=0
		n=self.seq_
		while c<len(n):
			fo.write(n[c:c+60])
			fo.write("\n")
			c=c+60
		#fo.close()
	def __len__(self):
		return len(self.seq_)
	def __str__(self):
		sq = self.name_+","+self.seq_
		return sq
	def __getitem__(self,i):
		return self.seq_[i]
	def __gt__(self, sq):
		return len(self.seq_)>len(sq)
	def __lt__(self, sq):
		return len(self.seq_)<len(sq)
	def __ge__(self, sq):
		return len(self.seq_)>=len(sq)
	def __le__(self,sq):
		return len(self.seq_)<=len(sq)
	def x_percent(self,ch):
		xc = self.seq_.count(ch)
		return (xc * 100 ) / self.length()
	def get_segment(self, st, en):
		ss = ""
		if st <= en:
			ss=self.seq_[st-1:en]
		else:
			raise ValueError("Start is greater than end.")
		return ss
