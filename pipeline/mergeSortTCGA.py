import glob
import os
import numpy

def main():
    x = 0
    os.chdir("TCGA")
    numFiles = len(glob.glob("*.txt"))
    numVals = 0
    #Establish the book-keeping columns
    files = glob.glob("*.txt")[0]
    #Count lines, create matrix
    with open(files,"r") as temp:
        for line in temp:
            numVals = numVals + 1
        matrix = [["NA" for x in xrange(numFiles+4)] for x in xrange(numVals-1)]
        matrix[0][0] = "Chromosome"
        matrix[0][1] = "Location"
        matrix[0][2] = "Reference"
        matrix[0][3] = "Gene"
        #initialize starting locations
        emptyRow = 1
        emptyControl = 4
        emptyTreated = numFiles+3
    #Populate first columns
    with open(files,"r") as temp:
        cnt = 1
        for line in temp:
            if cnt > 2:
                parts = line.split()
                if len(parts) == 5:
                    #chrom
                    matrix[emptyRow][0] = parts[3]
                    #loc
                    matrix[emptyRow][1] = parts[4]
                    #ref
                    matrix[emptyRow][2] = parts[0]
                    #gene
                    matrix[emptyRow][3] = parts[2]
                else:
                    matrix[emptyRow][0] = parts[2]
                    matrix[emptyRow][1] = parts[3]
                    matrix[emptyRow][2] = parts[0]
                    matrix[emptyRow][3] = "NA"
                #advance row
                emptyRow = emptyRow + 1
            #advance line
            cnt = cnt + 1
    #Populate data columns
    complete = 0.0
    for files in glob.glob("*.txt"):
        print '%Complete:' + repr(100.0*complete/numFiles)
        with open(files,"r") as myfile:
            cnt = 1
            for line in myfile:
                if cnt == 1:
                    sample = line.split()[2]
                    code = sample.split("-")[3]
                    code = code[0:2]
                    if int(code) > 9:
                        index = emptyControl
                        emptyControl = emptyControl + 1
                    else:
                        index = emptyTreated
                        emptyTreated = emptyTreated - 1
                    matrix[0][index] = sample
                if cnt > 2:
                    matrix[cnt-2][index] = line.split()[1]
                cnt = cnt + 1
            complete = complete + 1.0
    print("Merge Complete. Outputting File.")
    with open("finalOutput.txt","w") as output:
        for i in xrange(len(matrix)):
            myLine = matrix[i][0]
            for j in xrange(1,len(matrix[0])):
                myLine = myLine + " " + matrix[i][j]
            output.write(myLine + "\n")
    print("File Created.")

if __name__ == '__main__':
    main()
