# coding: utf-8
# check python version
import sys
if sys.version_info[0] < 3:
    sys.stdout.write("Sorry, requires Python 3.x\n")
    sys.exit(1)

# In[21]:

import re
import sys
import os
import csv
import itertools
import numpy


# In[22]:
gtfin = sys.argv[1]
gtfout = sys.argv[2]

print("Processing " + gtfin)

fp = open(gtfin) # open file on read mode
lines = fp.read().split("\n") # create a list containing all lines
fp.close() # close file


# In[23]:

#def add gene_name and ref_gene_id to field 
def ExtractUpdateFeature(x):
  x = re.split('\t| ',x.rstrip())
  if x[2]=="transcript":
    if len(x)==12:
      x.extend(["gene_name",x[9],"ref_gene_id",x[9]])
  if x[2]=="exon":
    if len(x)==14:
      x.extend(["gene_name",x[9],"ref_gene_id",x[9]])
    elif len(x)==16:
      x.extend(["ref_gene_id",x[9]])
  return x

#pretty head for dictionary, lazy!
def glance(myDict):
  return dict(list(myDict.items())[0:5])

#count unique gene_id and unique ref_gene_ids in the feature list
def CountFeature(x):
  # input is list of lines from GTF parse, 16 fields for transcript and 18 fields for exons
  r = []
  g = []
  for i in range(0,len(x)-1):
    if isinstance(x[i], str):
      GTFfeature = ExtractUpdateFeature(x[i])#re.split('\t| ',x[i].rstrip())
    elif isinstance(x[i], list):
      GTFfeature = list(x[i])
    refGeneId = GTFfeature[len(GTFfeature)-1]
    geneId = GTFfeature[9]
    #if len > 1 then take refgeneid and geneid if not update refgeneid and geneid according to len
    r.append(refGeneId)
    g.append(geneId)
  return len(list(set(r))),len(list(set(g)))

def CountENS(x):
  return len([i for i in x if i[0:4]=='"ENS'])

def CountMST(x):
  return len([i for i in x if i[0:4]=='"MST'])


# In[24]:

## Count number of exon and transcript features in the original GTF
lines = list(filter(None, lines))
print("Total number of features:", len(list(filter(None, lines))))
print("Total number of transcript features:",
      len([x for x in lines if ExtractUpdateFeature(x)[2]=="transcript"]))
print("Total number of exon features:",
      len([x for x in lines if ExtractUpdateFeature(x)[2]=="exon"]))


# In[25]:

## Create a dictionary mapping old and newly created ref_gene_id (last item) to gene_id (10th field)
featureMapNovel = {}
for i in range(0, len(lines)-1):
    GTFfeature = list(ExtractUpdateFeature(lines[i]))
    if len(GTFfeature) == 16:
        refGeneId = GTFfeature[len(GTFfeature)-1]
        geneId = GTFfeature[9]
        
        featureMapNovel.setdefault(geneId,[]).append(refGeneId)
        featureMapNovel[geneId] = list(set(featureMapNovel.get(geneId)));

glance(featureMapNovel)
# In[27]:

newGTF = []
ENS = list(numpy.zeros(4))
MST = list(numpy.zeros(6))
c = 0.0
f = 0
for i in range(0,len(lines)):
    GTFfeature = list(ExtractUpdateFeature(lines[i]))
    refGeneId = GTFfeature[len(GTFfeature)-1]
    geneId = GTFfeature[9]
    newGTF.append(list(GTFfeature))
    #total number of multi gene_id <-> re_gene_id relation
    if newGTF[i][2]=="transcript":
        c += 1
    
    if (geneId[0:4]=='"ENS' and refGeneId[0:4]=='"ENS' and
        len(featureMapNovel[geneId])==1):
        ENS[0] += 1
#        if ENS[0]==1:
#            print("ENS-ENS_unique: ",newGTF[i])
    elif (geneId[0:4]=='"ENS' and refGeneId[0:4]=='"ENS' and
          len(featureMapNovel[geneId])>1):
        ENS[1] += 1
#        if ENS[1]==1:
#            print(newGTF[i])
    elif (geneId[0:4]=='"ENS' and refGeneId[0:4]=='"MST' and
          len(featureMapNovel[geneId])==1):
        ENS[2] += 1
#        if ENS[2]==1:
#            print(newGTF[i])
    elif (geneId[0:4]=='"ENS' and refGeneId[0:4]=='"MST' and
          len(featureMapNovel[geneId])>1):
        ENS[3] += 1
#        if ENS[3]==1:
#            print(newGTF[i])
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"MST' and
          len(featureMapNovel[geneId])==1):
        MST[0] += 1
#        if MST[0]==1:
#            print("MST-MST_unique: ",newGTF[i])
        #Novel gene and novel transcript. No update needed.
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"MST' and
          (geneId in featureMapNovel[geneId]) and len(featureMapNovel[geneId])>1):
        MST[1] += 1
#        if CountENS(featureMapNovel[geneId])>1 and f==0:
#            print("MST to MST, mutliple (including MST and ENS): ",newGTF[i])
#            f=1
        #if there is only one ENS, the new novel transcript is for the annotated gene. Since this might be an
        #ambiguity, and we are not interested in novel isoforms, the MST feature should be removed.
        if (CountENS(featureMapNovel[geneId])==1 and CountMST(featureMapNovel[geneId])==1):
            newGTF[i] = list()
            
        #if there are two ENS, remove MST feature as well. The reasonf for multiple condition is future expansions
        elif (CountENS(featureMapNovel[geneId])>1 and CountMST(featureMapNovel[geneId])==1):
            newGTF[i] = list()
        
        #if there are more than one MST, throw a warning and break the loop. this calls for an investigation
        elif CountMST(featureMapNovel[geneId])>1:
            print(newGTF[i])
            print("More than one novel features, investigation needed. Breaking the loop")
            break
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"MST' and
          not(geneId in featureMapNovel[geneId]) and len(featureMapNovel[geneId])>1):
        MST[2] += 1
#        if MST[2]==2000:
#            print("MST to MST, multiple (only other MSTs): ", newGTF[i])
        #if this case is found, remove it. This should _NOT_ happen.
        #increment the count, remove it, and continue
        newGTf[i] = list()
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"ENS' and
          len(featureMapNovel[geneId])==1):
        MST[3] += 1
#        if MST[3]==1:
#            print("MST to ENS, one to one: ",newGTF[i])
        #Well, this is a bit weird. There are no novel genes or isoforms, but ref_gene_id "ENS" now has a new
        #gene_id. probably due to a novel isoform. Updating gene_id to ref_gene_id should fix the issue.
        newGTF[i][9] = refGeneId
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"ENS' and
          (geneId in featureMapNovel[geneId]) and len(featureMapNovel[geneId])>1):
        MST[4] += 1
#        if MST[4]==2000:
#            print("MST to MST, mutliple (MST to itself and other ENS but not other MSTs): ", newGTF[i])
        #This is a rather special case. Multiple conditions have to be check.
        #1. If there is another MST in the featureMapNovel this is red flag. this calls for a closer inspection.
        #   Breakt the loop and print feature.
        if CountMST(featureMapNovel[geneId])>1:
            print(newGTF[i])
            print("More than one novel features mapped to eachother, investigation needed. Breaking the loop")
            break
            
        #2. If there is only one ENS and of course obviously one MST itself, just update gene_id to ref_gene_id.
        #   This is similar to one of the cases above, where MST<->ENS exists. MST<->MST for this case is handelded
        #   above and removed promptly. Here we just update it back to ref_gene_id.
        #3. If there are more than one ENS, then just update gene_id to ref_gene_id
        #Since each feature is inspected independently, we can just update them independently as well.
        #So here if there is a MST<->anotherMST case, it will break the loops, otherwise, we update it
        #to ref_gene_id.
        newGTF[i][9] = refGeneId
        
    elif (geneId[0:4]=='"MST' and refGeneId[0:4]=='"ENS' and
          not(geneId in featureMapNovel[geneId]) and len(featureMapNovel[geneId])>1):
        MST[5] += 1
#        if MST[5]==20000:
#            print("MST to MST, mutliple (MST to other ENSs only, and not itself nor other MSTs: ", newGTF[i])
        #This case is just a usual case, transcript assembler somehow merges nearby annotated genes.
        #This exists due to failure in transcriptome assembler. Just update gene_id to ref_gene_id in each case.
        newGTF[i][9] = refGeneId
        
    else:
        print(newGTF[i])
        print("None of the conditions were met")
        print(featureMapNovel[geneId])
        sys.exit()
        break
            

# It's possible that one gene_name has two ref_gene_id.
# This is not a dealbreaker, as it only updates the gene_id based on novel IDs (i.e. MSTRGXXXX)
# so assign ref_gene_id anyway and resolve it later on.
# Note: Cases like these, where there is a single gene name for two ref_gene_id, are not unique.
#       and they are not considered as anything new.
print("number of features in new GTF: ", len(newGTF))
print("number of features in original GTF: ", len(lines))
newGTF = list(filter(None, newGTF))
print("number of features in new GTF after clean up: ", len(newGTF))


# In[28]:

print("ENS to ENS, unique one to one: ", ENS[0])
print("ENS to ENS, mutliple: ", ENS[1])
print("ENS to MST, unique one to one: ", ENS[2])
print("ENS to MST, mutliple: ", ENS[3])

print("MST to MST, unique one to one: ", MST[0])
print("MST to MST, mutliple (including MST and ENS): ", MST[1])
print("MST to MST, multiple (only other MSTs): ", MST[2])
print("MST to ENS, one to one: ", MST[3])
print("MST to MST, mutliple (MST to itself and other ENS but not other MSTs): ", MST[4])
print("MST to MST, mutliple (MST to other ENSs only, and not itself nor other MSTs: ", MST[5])

print("ENS features's array:")
print(ENS)

print("MST features's array:")
print(MST)

print("Total features processed:")
print(sum(ENS)+sum(MST))


# In[29]:

print("ref_gene_id and gene_id after processing: ", CountFeature(newGTF))
print("ref_gene_id and gene_id before processing: ", CountFeature(lines))


# In[30]:

print("Final number of features (exon+transcript) \n after removing redundant ones: ", len(newGTF))


# In[31]:

with open(gtfout, mode='tw') as myfile:
    for L in newGTF:
        myfile.write("\t".join(L[0:8])+"\t"+" ".join(L[8:len(L)])+"\n")
# In [32]:

print("Output is written to: " + gtfout)
