import os
import glob
import sys
import getopt
import argparse
import pandas as pd
import time


#print(sys.argv)

# Set up argument parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser()
   
parser.add_argument('-b','--bam', required=True, help="Bam files to downsample, can use *")
parser.add_argument('-d','--depth', required=True, help="Number to downsample by")
parser.add_argument('-p','--primers', required=True, help="list of primer sequences by newlines")

args = parser.parse_args() 

# Save values from parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

newPath = os.path.join(os.getcwd(), args.bam)

bamsPaths = glob.glob(newPath)
depth = int(args.depth)
primersPath = os.path.join(os.getcwd(), args.primers)

#print(glob.glob(newPath))

# Create a folder to save results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

os.system("mkdir downsampled_bams")

# Create a sam file from the bam file and move it to the new folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aFiles = []
for files in bamsPaths:

    fileName = files.split("/")[-1]

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
        matches = [row[1][1], row[1][2]]
        #print(row["forward"][0])
        #print(matches)
        #print(row["target"][0])
        primers[row[1][0]] = matches
    return primers

p = getPrimers(primersPath)  
#print(p)     

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

    df = df[0].str.split('\t', expand=True).iloc[commentCount:, :]
    #print(df.head())

    # print("Shape: " + str(df.shape))
    # print(df.head())

    # print(p1[1][1])

    # can improve with
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/text.html

    filtered = []
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    for key in primers.keys():
        prime_filter = df[df[9].str.contains(p[key][0], na=False) | df[9].str.contains(p[key][1], na=False)].head(depth)
        #print(prime_filter.head())
        filtered.append(prime_filter)
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    downsampled = pd.concat(filtered)
    #print(downsampled)
    #print(downsampled.shape)

    newFilePath = filePath.replace(".sam", "f.sam")

    downsampled.to_csv(newFilePath, index=False, sep="\t", encoding='utf-8', header=None)
    print(downsampled.shape)
    print(downsampled[downsampled[9].str.len() == downsampled[10].str.len()].shape)
    
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

    samFile = os.path.join("downsampled_bams/", file.replace(".sam", "f.sam"))
    samFile = os.path.join(os.getcwd(), samFile)
    print(samFile)

    bamFile = os.path.join("downsampled_bams/", file.replace("f.sam", "f.bam"))
    bamFile = os.path.join(os.getcwd(), bamFile)
    print(bamFile)

    please = os.system("samtools view -bS {} > {}".format(samFile, bamFile))
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~