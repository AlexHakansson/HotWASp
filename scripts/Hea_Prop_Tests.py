import numpy as np
import os
from scipy.stats.stats import pearsonr

def scoreOverlap (genes, heat,overlap , top =1000, outFile = None):
    
    hgTup= []
    with open(genes) as data:
        geneList = [x.strip() for x in data.readlines()]
    
    with open(heat) as data:
        heatList = [x.strip() for x in data.readlines()]
        
    with open(overlap) as data:
        overlapList = [x.strip() for x in data.readlines()]
    
    for i in range(len(heatList)):
        hgTup.append((heatList[i], geneList[i]))
        
    hgTup.sort(reverse=True)
    
    if top == "all":
        top = len(hgTup)
    
    count = 0
    if outFile != None:
        out = open(outFile,"w")
        
    for  i in hgTup[:top]:
        if i[1] in overlapList:
            count +=1
        if outFile != None:
            out.write(i[1] + "\n")
    
    return count
    
def rankCorrelation (genes, heat, compareRank, top = 100):

    hgTup= []
     
    with open(genes) as data:
        geneList = [x.strip() for x in data.readlines()]
    
    with open(heat) as data:
        heatList = [x.strip() for x in data.readlines()]
        
    with open(compareRank) as data:
        compareRank = [x.strip() for x in data.readlines()]
        
    for i in range(len(heatList)):
        hgTup.append((heatList[i], geneList[i]))
        
    hgTup.sort(reverse=True)
    
    geneToRank = {}
    
    for i in range(len(hgTup)):
        geneToRank[hgTup[i][1]] = i
        
    r1List = []
    r2List = []
    for i in range(top):
    
        g = compareRank[i]
        if g not in geneToRank:
            continue
        else:
            r1List.append(i)
            r2List.append(geneToRank[g])
    
    r1List = [x for x in r1List]
    r2List = [x for x in r2List]
    return pearsonr(r1List,r2List)[0]
    
    
    
def runTestsOnFolder (folder,outFile):
    files = [folder + x for x in os.listdir(folder)]
    
    gwasRankFile = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\Original_GWAS_Ordered.txt"
    netWASFile = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\NET_WAS_Positive_GWAS.txt"
    secondGwas = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\GWAs.txt"
    geneList = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\Original_GWAS_Ordered.txt"
    
    with open(outFile,"w") as out:
        out.write("Time\tOriginal_GWAS_Correlation\tNetWas_Overlap\tOther_Gwas_Overlap\n")
        for f in files:
            name = os.path.basename(f).split("_")[-1].split(".tx")[0]
            cor = rankCorrelation(geneList,f,gwasRankFile)
            net_over = scoreOverlap(geneList,f,netWASFile)
            gwas_over = scoreOverlap(geneList,f,secondGwas)
            out.write("%s\t%f\t%i\t%i\n"%(name,cor,net_over,gwas_over))
          
          
def multiTimeTest (folder, overlap_out,correl_out):

    files = [folder + x for x in os.listdir(folder)]
    
    gwasRankFile = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\Original_GWAS_Ordered.txt"
    netWASFile = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\NET_WAS_Positive_GWAS.txt"
    secondGwas = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\GWAs.txt"
    geneList = "E:\\Alex_MM_tmp\\GWAS_Analysis\\Test_Files\\Original_GWAS_Ordered.txt"
    
    number = [1000,750,500,250,100,50,10]
    
    GWAS_File = "hi"
    with open(overlap_out,"w") as out:
        out.write("Time\t" + "\t".join([str(x) for x in number]) + "\n")
        for f in files:
            name = os.path.basename(f).split("_")[-1].split(".tx")[0]
            out.write(name+"\t")
            for n in number:
                
                #cor = rankCorrelation(geneList,f,gwasRankFile, )
                #net_over = scoreOverlap(geneList,f,netWASFile)
                gwas_over = scoreOverlap(geneList,f,secondGwas, top = n)
                out.write(str(gwas_over)+"\t")
            out.write("\n")
            
    with open(correl_out,"w") as out:
        out.write("Time\t" + "\t".join([str(x) for x in number]) + "\n")
        for f in files:
            name = os.path.basename(f).split("_")[-1].split(".tx")[0]
            out.write(name+"\t")
            for n in number:
                
                cor = rankCorrelation(geneList,f,gwasRankFile, top = n )
                #net_over = scoreOverlap(geneList,f,netWASFile)
                #gwas_over = scoreOverlap(geneList,f,secondGwas)
                out.write(str(cor)+"\t")
            out.write("\n")

def makeMax(folder):
    
    files = [folder + x for x in os.listdir(folder)]
    maxList = []
    outFile = folder + os.path.basename(folder) +"_max.txt"
    with open(outFile,"w") as out:
        out.write("Time\tOriginal_GWAS_Correlation\tNetWas_Overlap\tOther_Gwas_Overlap\n")
        for f in files:
            curr = [x.strip() for x in open(f).readlines()]
            if len(maxList) ==0:
                for i in curr:
                    maxList.append([i])
            else:
                for i in range(len(curr )):
                    maxList[i].append(curr[i])
                    
        for i in maxList:
            out.write(str(max(i))+"\n")
