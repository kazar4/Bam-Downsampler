import os
import glob
import sys
import getopt
import argparse
import pandas


# required arg

print(sys.argv)

parser = argparse.ArgumentParser()
   
parser.add_argument('-b','--bam', required=True, help="Bam files to downsample, can use *")
parser.add_argument('-d','--depth', required=True, help="Number to downsample by")
parser.add_argument('-p','--primers', required=True, help="list of primer sequences by newlines")

args = parser.parse_args() 

print(f'Hello {args.bam}')

newPath = os.path.join(os.getcwd(), args.bam)
print("BRUH" + os.getcwd())
print(newPath)

bamsPaths = glob.glob(newPath)
depth = int(args.depth)
print(depth)
primersPath = os.path.join(os.getcwd(), args.primers)

print(glob.glob(newPath))

os.system("mkdir downsampled_bams")

aFiles = []
for files in bamsPaths:

    fileName = files.split("/")[-1]

    samFile = os.path.join("downsampled_bams/", fileName.replace('.bam', '.sam'))
    samFile = os.path.join(os.getcwd(), samFile)
    aFiles.append(samFile)
    print(samFile)

    test = os.system("samtools view -h -o {} {}".format(samFile, files))

print(aFiles)

def getPrimers(filePath):
    test = pandas.read_csv(filePath, sep="\t")
    #test.head()
    primers = {}
    print(len(list(test.items())))
    for row in test.iterrows():
        #print(row[1][0])
        matches = [row[1][1], row[1][2]]
        #print(row["forward"][0])
        #print(matches)
        #print(row["target"][0])
        primers[row[1][0]] = matches
    return primers

p = getPrimers(primersPath)  
print(p)      

def doThings(filePath, primers):
    downsampleCounts = {}
    #test = pandas.read_csv(filePath, sep="\t")
    file1 = open(filePath, 'r+')
    count = 0

    ok = file1.readlines()
    print("lines of file: " + str(len(ok)))

    
    keyCount = {}
    for prime in primers.keys():
        print(primers[prime])
        keyCount[prime] = 0

    index = 0
    rc = 0
    for row in ok:

        toDelete = True
        #print(row)

        if row[0] == "@":
            index = index + 1
            continue

        for key in primers.keys():
            pMatch = primers[key]
            if ((pMatch[0] in row) or (pMatch[1] in row)):
                if keyCount[key] < depth or True:
                    #toDelete = False
                    keyCount[key] = keyCount[key] + 1
                    print(key + " : " + str(index))
                    count = count + 1
                    #break
        if toDelete:
            #ok.remove(row)
            rc = rc + 1;
            ok.pop(index)
    
        index = index + 1

    print(count)
    print(keyCount)
    print("rc: " + str(rc))
    file1.close()
    file2 = open(filePath, 'w')
    file3 = open(filePath.replace(".sam", "f.sam"), 'w+')
    file3.writelines(ok)
    print("lines of new file: " + str(len(ok)))
    file2.writelines(ok)
    
    
for file in aFiles:
    doThings(file, p)

    print(file)

    samFile = os.path.join("downsampled_bams/", file)
    samFile = os.path.join(os.getcwd(), samFile)
    print(samFile)

    bamFile = os.path.join("downsampled_bams/", file.replace(".sam", ".bam"))
    bamFile = os.path.join(os.getcwd(), bamFile)
    print(bamFile)

    please = os.system("samtools view -bS {} > {}".format(samFile, bamFile))


"""
# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]

print(argumentList)

# Options
options = "bdp:"

# Long options
long_options = ["Bam", "Depth", "Primers"]

try:
    # Parsing argument
    arguments, values = getopt.getopt(argumentList, options, long_options)
    print(arguments)
     
    # checking each argument
    for currentArgument, currentValue in arguments:

        print(currentArgument, currentValue)
 
        if currentArgument in ("-b", "--Bam,"):
            print(currentValue)
             
        if currentArgument in ("-d", "--Depth"):
            print(currentValue)
             
        if currentArgument in ("-p", "--Primers"):
            print(currentValue)
             
except getopt.error as err:
    # output error, and return with an error code
    print(str(err))
    print("got here")
"""



import glob
path = "/home/mydir/*.txt"
for filename in glob.glob(path):
    with open(filename, 'r') as f:
        for line in f:
            print(line)


# Take -list of bam files or a folder
# Take a depth number
# Take a list of Primers

# create a folder in this location of the new files to save

# Create a dictionary of primers (where each key contains forward and reverse)
# go through each line of each file
# if it contains primer key, and is dictonary value < depth, keep it 
# otherwise remove from file
# save new file in folder with same name as orginal
# Go onto next file and repeat
