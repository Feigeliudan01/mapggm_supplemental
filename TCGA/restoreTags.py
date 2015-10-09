import sys

def main():
	outfile = open('cleanedFinal.txt','w')
	with open('premerged.txt','r') as infile:
		outfile.write(infile.readline() + "\n")
		for line in infile:
			splitLine = line.split()
			tag = splitLine[0]
			for n in xrange(8-len(tag)):
				tag = "0" + tag
			tag = "\"cg" + tag + "\""
			nextLine = tag
			for i in xrange(1,len(splitLine)):
				nextLine = nextLine + " " + splitLine[i]
			nextLine = nextLine + "\n"
			outfile.write(nextLine)
	outfile.close()
if __name__ == '__main__':
	main()
