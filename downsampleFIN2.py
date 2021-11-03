import os
import glob
import sys
import getopt
import argparse
import pandas as pd
import time
import csv


#print(sys.argv)

# Set up argument parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser()
   
parser.add_argument('-b','--bam', required=True, nargs='+', help="Bam files to downsample, can use *")
parser.add_argument('-d','--depth', required=True, help="Number to downsample by")
parser.add_argument('-p','--primers', required=True, help="list of primer sequences by newlines")
parser.add_argument('-o', '--output', help="if only one bam is given the output name can be chosen")

args = parser.parse_args() 

# Save values from parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# combines full path with bam path
bamsPaths = list(map(lambda x: os.path.join(os.getcwd(), x), args.bam))
depth = int(args.depth)
primersPath = os.path.join(os.getcwd(), args.primers)
outputName = args.bam
print(outputName)
if len(bamsPaths) > 1 and outputName != "":
    outputName = ""
    print("Output name can only be set if only one bam file is given, output argument will be ignored")


# Create a folder to save results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

os.system("mkdir downsampled_bams")

# Create a sam file from the bam file and move it to the new folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aFiles = []
for files in bamsPaths:

    fileName = files.split("/")[-1]

    if (outputName != ""):
        fileName = outputName

    samFile = os.path.join("downsampled_bams/", fileName.replace('.bam', '.sam'))
    samFile = os.path.join(os.getcwd(), samFile)
    aFiles.append(samFile)
    print(samFile)

    test = os.system("samtools view -h -o {} {}".format(samFile, files))

#print(aFiles)

# Get the primers to filter by
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def getPrimers(filePath):
    test = pd.read_csv(filePath, sep="\t")
    #test.head()
    primers = {}
    print(len(list(test.items())))
    for row in test.iterrows():
        #print(row[1][0])
        matches = [row[1][1].split(","), row[1][2].split(",")]
        #print(row["forward"][0])
        #print(matches)
        #print(row["target"][0])
        primers[row[1][0]] = matches
    return primers

p = getPrimers(primersPath)  
print(p)     

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
def doThings(filePath, primers):

    df = pd.read_csv("HLA-LA2/200918_barcode08.sam", header=None, sep='\n')

    commentCount = headerRowCount(df, "@")
    if commentCount < 0:
        commentCount = 0 

    comments = df.head(commentCount)
    #print(comments)

    #df = df[0].str.split('\t', expand=True).iloc[commentCount:, :]
    df = df.iloc[commentCount:, :]
    #print(df.head())

    # print("Shape: " + str(df.shape))
    # print(df.head())

    # print(p1[1][1])

    # can improve with
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/text.html

    filtered = []
    p_filter = df

    def containsPrimer(series, key):
        newArray = []
        for x in range(0, len(series)):
            anyMatch = False;
            for pNum in range(0, len(p[key][0])):
                anyMatch = anyMatch or (p[key][0][pNum] in series[x])
                anyMatch = anyMatch or (p[key][1][pNum] in series[x])

                #anyMatch.append(x.str.contains(p[key][0][pNum], na=False)) 
                #anyMatch.append(x.str.contains(p[key][1][pNum], na=False))
            #print(type(anyMatch))
            newArray.append(anyMatch)
        return newArray

    def cPrimer(df, key):
        #df[df[9].str.contains(p[key][0], na=False) | df[9].str.contains(p[key][1], na=False)]
        statement = df[0].str.contains(p[key][0][0], na=False) | df[0].str.contains(p[key][1][0], na=False)
        for pNum in range(1, len(p[key][0])):
            statement = statement | df[0].str.contains(p[key][0][pNum], na=False)
        for pNum in range(1, len(p[key][1])):
            statement = statement | df[0].str.contains(p[key][1][pNum], na=False)

        return statement

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    for key in primers.keys():
        p_filter = df[cPrimer(df, key)]
        #df[df.apply(lambda x: containsPrimer(x, key))]
        #df[containsPrimer(df[0], key)]
        #print(prime_filter.head())
        p_filter = p_filter.head(depth)
        print(p_filter.shape)
        filtered.append(p_filter)
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    downsampled = pd.concat(filtered)
    #print(downsampled)
    #print(downsampled.shape)

    newFilePath = filePath.replace(".sam", "f.sam")

    # sep="\t"
    downsampled.to_csv(newFilePath, index=False, sep="\n", encoding='utf-8', header=None, quoting=csv.QUOTE_NONE, escapechar = '\n')
    print(downsampled.shape)
    #print(downsampled[downsampled[9].str.len() == downsampled[10].str.len()].shape)
    
    f = open(newFilePath, 'r')
    original = f.readlines()

    f = open(newFilePath, 'w+')
    print(comments.shape)
    for row in comments.iterrows():
        #print(row[1][0])
        f.write(row[1][0] + '\n')
     
    f.writelines(original) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Main logic, runs previous method and then converts the new filtered sam files back into bam
for file in aFiles:
    doThings(file, p)

    print(file)

    filterStr = "f"
    if (outputName != ""):
        filterStr = ""

    samFile = os.path.join("downsampled_bams/", file.replace(".sam", filterStr + ".sam"))
    samFile = os.path.join(os.getcwd(), samFile)
    print(samFile)

    bamFile = os.path.join("downsampled_bams/", file.replace(".sam", filterStr + ".bam"))
    bamFile = os.path.join(os.getcwd(), bamFile)
    print(bamFile)

    please = os.system("samtools view -bS {} > {}".format(samFile, bamFile))
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~