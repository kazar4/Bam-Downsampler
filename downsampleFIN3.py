import os
import argparse
import pandas as pd
import csv

# Just wanted to add cool colors :|
# https://stackoverflow.com/questions/287871/how-to-print-colored-text-to-the-terminal
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Set up argument parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#print(sys.argv)

parser = argparse.ArgumentParser()
   
parser.add_argument('-b','--bam', required=True, nargs='+', help="Bam files to downsample, can use *")
parser.add_argument('-d','--depth', required=True, help="Number to downsample by")
parser.add_argument('-p','--primers', required=True, help="list of primer sequences by newlines")
parser.add_argument('-o', '--output', help="if only one bam is given the output name can be chosen, include .bam in name")

args = parser.parse_args() 

# Save values from parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# combines full path with bam path
bamsPaths = list(map(lambda x: os.path.join(os.getcwd(), x), args.bam))
depth = int(args.depth)
primersPath = os.path.join(os.getcwd(), args.primers)
outputName = args.output

if len(bamsPaths) > 1 and outputName != None:
    outputName = None
    print("")
    print("Output name can only be set if only one bam file is given, output argument will be ignored")

# Create a folder to save results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

os.system("mkdir downsampled_bams")

# Create a sam file from the bam file and move it to the new folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aFiles = []
print("")
print("Files found: ")
for file in bamsPaths:

    fileName = file.split("/")[-1]

    samFile = os.path.join("downsampled_bams/", fileName.replace('.bam', '.sam'))
    samFile = os.path.join(os.getcwd(), samFile)
    aFiles.append(samFile)

    print(bcolors.BOLD + samFile + bcolors.ENDC)

    create_sam = os.system("samtools view -h -o {} {}".format(samFile, file))

print("")

#print(aFiles)

# Get the primers to filter by
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def getPrimers(filePath):
    test = pd.read_csv(filePath, sep="\t")
    #test.head()
    primers = {}
    #print(len(list(test.items())))
    
    for row in test.iterrows():
        matches = [row[1][1].split(","), row[1][2].split(",")]

        primers[row[1][0]] = matches
    return primers

p = getPrimers(primersPath)

# Filters the sam file by depth and primers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gets number of rows that are headers
def headerRowCount(data, delim):
    headerCount = 0
    for row in data.iterrows():
        if (row[1][0])[0] == delim:
            headerCount = headerCount + 1
        else:
            break
    return headerCount

# main filtering logic
def filterSam(filePath, primers):

    df = pd.read_csv(filePath, header=None, sep='\n')

    commentCount = headerRowCount(df, "@")
    if commentCount < 0:
        commentCount = 0 

    comments = df.head(commentCount)

    #df = df[0].str.split('\t', expand=True).iloc[commentCount:, :]
    df = df.iloc[commentCount:, :]

    # can improve with
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/text.html

    filtered = []
    p_filter = df

    # Creates boolean series to check if the dataframe row contains a primer sequence
    def cPrimer(df, key):
        statement = df[0].str.contains(p[key][0][0], na=False) | df[0].str.contains(p[key][1][0], na=False)
        for pNum in range(1, len(p[key][0])):
            statement = statement | df[0].str.contains(p[key][0][pNum], na=False)
        for pNum in range(1, len(p[key][1])):
            statement = statement | df[0].str.contains(p[key][1][pNum], na=False)

        return statement

    for key in primers.keys():
        p_filter = df[cPrimer(df, key)].head(depth)
  
        print("Found " + str(p_filter.shape[0]) + " matches of " + key)
        filtered.append(p_filter)
    
    downsampled = pd.concat(filtered)

    newFilePath = filePath.replace(".sam", "f.sam")

    downsampled.to_csv(newFilePath, index=False, sep="\n", encoding='utf-8', header=None, quoting=csv.QUOTE_NONE, escapechar = '\n')
    
    f = open(newFilePath, 'r')
    original = f.readlines()

    f = open(newFilePath, 'w+')
    print("There were {} header rows".format(str(comments.shape[0])))
    for row in comments.iterrows():
        #print(row[1][0])
        f.write(row[1][0] + '\n')
     
    f.writelines(original) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Main logic, runs previous method and then converts the new filtered sam files back into bam
fileCount = 0
for file in aFiles:
    fileCount = fileCount + 1

    filterSam(file, p)

    bamFileName = file

    filterStrSam = "f"
    filterStrBam = "f"
    if (outputName != None):
        filterStrBam = ""
        bamFileName = outputName.replace(".bam", ".sam")

    samFile = os.path.join("downsampled_bams/", file.replace(".sam", filterStrSam + ".sam"))
    samFile = os.path.join(os.getcwd(), samFile)
    print("Filtered Sam File: " + samFile)

    bamFile = os.path.join("downsampled_bams/", bamFileName.replace(".sam", filterStrBam + ".bam"))
    bamFile = os.path.join(os.getcwd(), bamFile)

    please = os.system("samtools view -bS {} > {}".format(samFile, bamFile))

    print("{}File {} Completed{}".format(bcolors.OKCYAN, str(fileCount), bcolors.ENDC))
    print("File Saved to: " + bcolors.OKGREEN + bamFile + bcolors.ENDC)
    print("")
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
