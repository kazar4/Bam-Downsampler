import pandas as pd

#df = pd.read_csv("HLA-LA2/200918_barcode08.sam", skiprows=2, header=None, sep='\t',
#    names=list(range(21)))

df = pd.read_csv("HLA-LA2/200918_barcode08.sam", header=None, sep='\n')

# print("Shape: " + str(df.shape))
# print(df.head())

def headerRowCount(data, delim):
    headerCount = 0
    for row in data.iterrows():
        if (row[1][0])[0] == delim:
            headerCount = headerCount + 1
        else:
            break
    return headerCount

print(headerRowCount(df, "@"))

df = df[0].str.split('\t', expand=True)

p1 = ("HLA-A", ["ATCCTGGATACTCACGACGCGGAC", "CATCAACCTCTCATGGCAAGAATTT"])
p2 = ("HLA-B", ["AGGTGAATGGCTCTGAAAATTTGTCTC", "AGAGTTTAATTGTAATGCTGTTTTGACACA"])
p3 = ("HLA-C", ["CAGCACGAAGATCACTGGAA", "TGAGGAAAAGGAGCAGAGGA"])

# print(p1[1][1])

# can improve with
# https://pandas.pydata.org/pandas-docs/stable/user_guide/text.html
p1df = df[df[9].str.contains(p1[1][0], na=False) | df[9].str.contains(p1[1][1], na=False)].head(25)
#print(p1df.shape)

#print(p1df)
#print(p1df.iloc[:, 2:])