import sys

f = open("head100000.fastq","r")
fout = open("head50000.fastq","w")
for i in range(0, 50000):
	fout.write(f.readline())
fout.close()
  
