def netListToAdcancyMatrix2 (netList, outFile):
    #netData = open(netList).readlines()
    
    source = [x[0] for x in netList]
    target = [x[1] for x in netList]
    
    geneList = list(set(source) | set(target))
    geneList.sort()
    
    del source
    del target
    
    length = len(netList)
    count = 0
    netDict = {}
    
    print("Dict Making")
    while len(netList) > 0:
        i = netList.pop()
        if count% 5000 == 0:
            print ("...%i/%i"%(count,length))
        count +=1
        if i[0] not in netDict:
            netDict[i[0]] = set()
            
        if i[1] not in netDict:
            netDict[i[1]] = set()
        
        netDict[i[0]] = i[1]
        netDict[i[1]] = i[0]
    
    print("writing")
    count = 0
    length = len(geneList)
    
    with open(outFile,"w") as out:
        out.write("gene\t"+"\t".join(geneList) +"\n")
        for i in geneList:
            if count% 5000 == 0:
                print ("...%i/%i"%(count,length))
            count +=1
            out.write(i + "\t")
            for k in geneList:
                if k in netDict[i]:
                    out.write("1\t")
                else:
                    out.write("0\t")
            out.write("\n")
       
    
def subnetConNeigh (network, genes,minConnections = None, outFile = None):

    netDict,netGenes  = rn(network)

def heatProp2 (adjMatrix, heat, prop = 3, time = 1):

    degree = np.sum(adjMatrix,1)
    heatProp = np.identity(adjMatrix.shape[0])
    laplace = heatProp
    
    print("Calculating Laplace")
    for i in range(len(degree)):
        laplace[i,i] = degree[i]
    
    laplace = adjMatrix - laplace
    laplace = laplace*time
    
    #lpTest = np.sum(laplace,1)
    
    lpNeg = []
    #for i in lpTest:
      #  if i < 0:
        #    lpNeg.append(i)
   # print (lpNeg)
            
    
    laplacePower = laplace
    
    p = 1
    while prop > 0:
        print("Calculating Propagation %i"%p)
        heatProp += laplacePower/ math.factorial(p)
        p += 1
        prop -= 1
        
        if prop != 0:
            laplacePower = np.matmul(laplacePower,laplace)
    
    print("Calculating new heat")
    newHeat = np.matmul(heat,heatProp)
    return newHeat
