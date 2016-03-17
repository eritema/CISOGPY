import optparse, sys
import networkx as nx
import matplotlib.pyplot as plt
#from  math import *
import math
import sys
import re
import fileinput
import scipy.special
from decimal import *


def exportGraph(extFile,G):
	for i in range(0,len(lista)):
		campi=campi+","+lista[i]
	extFile.write(str(source)+","+str(target)+","+campi+"\n")
	return


def getValue(valori,intensity):
	numCampi=len(intensity.split(','))
	campiNumerici=intensity.split(',')	
	if valori==[]:
		for i in range(3,numCampi):
			valori.append(0.0)
	for i in range(3,numCampi):
		if campiNumerici[i]!='':
			valori[i-3]=valori[i-3]+float(campiNumerici[i])
	return valori


def getExpectedValue100K(x,order,cons):
	tableMean=[[0,    5,20  ,200   , 600  ,2000],
                   [0,    1.5,7   ,70    ,200   ,600],
                   [0,0.3  ,1.5 ,15    ,55    ,160],
                   [0,0.1  ,0.3 ,3     ,11    ,30],
                   [0,0.04 ,0.06 ,0.6  ,2     ,6],
                   [0,0.01 ,0.012 ,0.12 ,0.7   ,2],    #7
                   [0,0.004,0.004 ,0.04 ,0.2   ,0.8],
                   [0,0    ,0    ,0.01 ,0.08  ,0.2],
                   [0,0    ,0    ,0.004,0.01  ,0.06],
                   [0,0    ,0    ,0.001,0.003 ,0.02]]
	tableMax=[[1,  33,  120, 600,1500,10000],
		[1,   3,   9,90,390,7500],
		[1, 1, 1.7, 15,100,3000],
		[1,1, 1,   10,45,1200],
		[1,1,1, 5,15,600],
		[1,   1,1,3,10,400],
		[1,   1,   1,1,5,150],
		[1,   1,   1,   1,3,50],	
		[1,   1,   1,   1,1,20],	
		[1,   1,   1,   1,1,10]]	
		
	if cons:
		table=tableMax
	else:
		table=tableMean
	if order<2:
		order=2
	if order>11:
		order=11

	if x<500:
		i=0
		x0=0.0
		x1=500.0
	elif x>=500 and x<1000:
		i=1
		x0=500.0
		x1=1000.0
	elif x>=1000 and x<2500:
		i=2
		x0=1000.0
		x1=2500.0
	elif x>=2500 and x<5000:
		i=3
		x0=2500.0
		x1=5000.0
	elif x>=5000 and x<20000:
                i=4
                x0=5000.0
                x1=20000.0
	expect=(x-x0)/(x1-x0)*(table[order-2][i+1]-table[order-2][i])+table[order-2][i]
	#print (table[order-2][i+1]-table[order-2][i])
	#return (table[order-2][i+1]-table[order-2][i])
	return expect
 
def getExpectedValue50K(x,order,cons):
		#   0,500  ,1000,2500 ,5000 ,20000
	tableMean=[[0,    5,20  ,200   , 600  ,2000],
                   [0,    1.5,7   ,70    ,200   ,600],
                   [0,0.3  ,1.5 ,15    ,55    ,160],
                   [0,0.1  ,0.3 ,3     ,11    ,30],
                   [0,0.04 ,0.06 ,0.6  ,2     ,6],
                   [0,0.01 ,0.012 ,0.12 ,0.7   ,2],    #7
                   [0,0.004,0.004 ,0.04 ,0.2   ,0.8],
                   [0,0    ,0    ,0.01 ,0.08  ,0.2],
                   [0,0    ,0    ,0.004,0.01  ,0.06],
                   [0,0    ,0    ,0.001,0.003 ,0.02]]
        tableMax=[[0,6,50,300,1000,2500],
                [0,2,5,40,150,600],
                [0,1,1,2,50,150],
                [0,1,1,1,3,15],
                [0,1,1,1,1,2],
                [0,1,1,1,1,1],
                [0,1,1,1,1,1],
                [0,1,1,1,1,1],
                [0,1,1,1,1,1],
                [0,1,1,1,1,1]]	
	if cons:
		table=tableMax
	else:
		table=tableMean
	if order<2:
		order=2
	if order>11:
		order=11

	if x<500:
		i=0
		x0=0.0
		x1=500.0
	elif x>=500 and x<1000:
		i=1
		x0=500.0
		x1=1000.0
	elif x>=1000 and x<2500:
		i=2
		x0=1000.0
		x1=2500.0
	elif x>=2500 and x<5000:
		i=3
		x0=2500.0
		x1=5000.0
	elif x>=5000 and x<20000:
                i=4
                x0=5000.0
                x1=20000.0
	elif x>=20000:
		i=4
                x0=5000.0
                x1=20000.0
	expect=(x-x0)/(x1-x0)*(table[order-2][i+1]-table[order-2][i])+table[order-2][i]
	#print (table[order-2][i+1]-table[order-2][i])
	#return (table[order-2][i+1]-table[order-2][i])
	return expect


def getErf(density,mu,sigma,priori):
	return (0.5+0.5*scipy.special.erf((density-mu)/(math.sqrt(2)*sigma)))*priori
	#return (0.5+0.5*scipy.special.erf((density-mu)/(math.sqrt(2)*sigma)))

def getCisProb(diameter,numberOfNodes,cisNumberOrder,numIS,cons,thr):
	mu50=[17600,15700,15700,14800,14500,9300,9100,8700,7800,7340]
	sigma50=[9700,7600,6500,6400,4800,4800]
	
	order=numberOfNodes
	if numIS>000:
          if numberOfNodes==3:
		iSigma=0
		iMu=0
	  if numberOfNodes==4:
		iSigma=1
		iMu=1
	  if numberOfNodes==5:
		iSigma=2
		iMu=2
	  if numberOfNodes==6:
		iSigma=3
		iMu=3
	  if numberOfNodes==7:
		iSigma=4
		iMu=4
	  if numberOfNodes>=8:
		iSigma=5
		iMu=4
		numberOfNodes=8
	else:
          if numberOfNodes==3:
		iSigma=0
		iMu=5
	  if numberOfNodes==4:
		iSigma=1
		iMu=6
	  if numberOfNodes==5:
		iSigma=2
		iMu=7
	  if numberOfNodes==6:
		iSigma=3
		iMu=8
	  if numberOfNodes==7:
		iSigma=4
		iMu=9
	  if numberOfNodes>=8:
		iSigma=5
		iMu=9
		numberOfNodes=8
	mu=mu50
	sigma=sigma50
	expectVal=float(getExpectedValue50K(numIS,order,cons))
	
	#priori=expectVal/(expectVal+cisNumberOrder[numberOfNodes])
	if numberOfNodes==1: 
		priori=expectVal
	elif expectVal>cisNumberOrder[numberOfNodes]:
		priori=1
	else:
		priori=expectVal/(cisNumberOrder[numberOfNodes])
		#print str(numberOfNodes)+" "+str(expectVal)+" "+str(priori)+" "+str(cisNumberOrder[numberOfNodes])
	
	if numberOfNodes>2:
		return getErf(diameter/(numberOfNodes-1),mu[iMu],sigma[iSigma],priori)		
		#return getExpectedValue(numIS,order) 
	else:
	    if int(thr)<=75000:
		return (float(diameter)/50000)*priori
	    return (float(diameter)/100000)*priori
		#return getExpectedValue(numIS,order)+cisNumberOrder[numberOfNodes]	

def getCisType(diameter,numberOfNodes,cisNumberOrder,numIS,cons,thr,verbose,minStatNode):
	mu50=[17600,15700,15700,14800,14500,9300,9100,8700,7800,7340]
	sigma50=[9700,7600,6500,6400,4800,4800]
	muExp=[23000,24000,25000,26000,26000,23000,24000,25000,26000,26000]
	#mu100=[48000,47000,46000,44500,44000]
	#sigma100=[20000,16600,15100,12200,10700,9600]
	#mu50=[23500,23000,23000,22500,22500]
	#sigma50=[20000,16600,15100,12200,10700,9600]
	#muExp=[17000,16000,14500,11500,11500]
	order=numberOfNodes
	if numIS<8000:
          if numberOfNodes==3:
		iSigma=0
		iMu=0
	  if numberOfNodes==4:
		iSigma=1
		iMu=1
	  if numberOfNodes==5:
		iSigma=2
		iMu=2
	  if numberOfNodes==6:
		iSigma=3
		iMu=3
	  if numberOfNodes==7:
		iSigma=4
		iMu=4
	  if numberOfNodes>=8:
		iSigma=5
		iMu=4
		numberOfNodes=8
	else:
          if numberOfNodes==3:
		iSigma=0
		iMu=5
	  if numberOfNodes==4:
		iSigma=1
		iMu=6
	  if numberOfNodes==5:
		iSigma=2
		iMu=7
	  if numberOfNodes==6:
		iSigma=3
		iMu=8
	  if numberOfNodes==7:
		iSigma=4
		iMu=9
	  if numberOfNodes>=8:
		iSigma=5
		iMu=9
		numberOfNodes=8
	#if int(thr)<=75000:
	mu=mu50
	sigma=sigma50
	
	expectVal=float(getExpectedValue50K(numIS,order,cons))
	
	#else:
	#	mu=mu100
        #        sigma=sigma100
	#	expectVal=float(getExpectedValue100K(numIS,order,cons))
	if verbose:
		print "Expected for a "+str(order)+" CIS: "+str(expectVal)
		print "Found CIS: "+str(cisNumberOrder[numberOfNodes])
	#priori0=expectVal/(expectVal+cisNumberOrder[numberOfNodes])
	if numberOfNodes==1: 
		priori0=0
		priori1=1
	elif expectVal>cisNumberOrder[numberOfNodes]:
		priori0=1
		priori1=0
	else:
		priori0=expectVal/(cisNumberOrder[numberOfNodes])
		priori1=float(cisNumberOrder[numberOfNodes]-expectVal)/(cisNumberOrder[numberOfNodes])
		#priori1=float(cisNumberOrder[numberOfNodes])/(expectVal+cisNumberOrder[numberOfNodes])
	
	if numberOfNodes<=minStatNode:
                CISvalue="LowOrder"
		#like0=0.01*priori0
		#like1=(-0.00000013*float(diameter)+0.013)*priori1
        else:
		like0=getErf(diameter/(numberOfNodes-1),mu[iMu],sigma[iSigma],priori0)
		like1=getErf(diameter/(numberOfNodes-1),muExp[iMu],21000,priori1)
        	if (like0/like1)>1:
			#CISvalue=str(like0/like1)
			#CISvalue=str(getExpectedValue(numIS,order,cons))+" "+str(cisNumberOrder[numberOfNodes])
                	CISvalue="Bad"
        	else:
			#CISvalue=str(getExpectedValue(numIS,order,cons))+" "+str(cisNumberOrder[numberOfNodes])
			#CISvalue=str(like0/like1)
                	CISvalue="Good"
		#print str(numberOfNodes)+" "+str(diameter/(numberOfNodes-1))+" "+str(like0)+" "+str(like1)+" "+str(priori0)+" "+str(priori1)
	return CISvalue,expectVal

def getEntropy(specie,specieNumber):
	tot=0
	p=0
	S=0
	for s in specie:
		tot=tot+specie[s]
	for s in specie:
		p=float(specie[s])/tot
		S=S+p*math.log(p)
	if len(specie)==1:
		return 0
	#S=-S/math.log(len(specie))
	S=-S/math.log(specieNumber)
	return S

def saveAttributes(comp,numberOfNodes,chromo,G,clusterPosition,diameter,outputFileColorTable,outputFileISTable,noLowOrder,stat,cisNumberOrder,numIS,cons,prob,thr,entropy,verbose,specieNumber,minStatNode,optGeneString):
	valori=[]
	valoriStr=""
	prec=[]
	nome=""
	primo=1
	specie={}
	S=0
	numberOfNodesApp=numberOfNodes
	for n in nx.nodes_iter(G):
		chromosome=G.node[n]['chromosome']
		position=G.node[n]['position']
		intensity=G.node[n]['intensity']
		gene=G.node[n]['gene']
		annot1=G.node[n]['ANNOTATION1']
		annot2=G.node[n]['ANNOTATION2']
		annot3=G.node[n]['ANNOTATION3']
		if optGeneString!="":
			if annot1==optGeneString:
				numberOfNodesApp=numberOfNodesApp-1	
		if (gene not in prec):
		   prec.append(gene)
		if (entropy):
			if(annot3 not in specie):
				specie[annot3]=1
			else:
				specie[annot3]=specie[annot3]+1
		valori=getValue(valori,intensity)
		if not stat:
			outputFileISTable.write(str(comp)+","+str(numberOfNodesApp)+","+str(chromosome)+","+str(position)+","+str(clusterPosition)+","+gene+","+annot1+','+annot2+','+annot3+','+intensity+"\n")
	prec.sort()
	primo=1
	if (entropy):
		S=getEntropy(specie,specieNumber)
	for name in prec:
		if primo:
                        nome=name
                        primo=0
                else :
                        nome=nome+"/"+name
	primo=1
	for valore in valori:
	    if primo:
		valoriStr=str(valore)
		primo=0
	    else:
		valoriStr=valoriStr+","+str(valore)
		#print ("\tcluster->%i,%s" % (comp,intensity))
	[CISvalue,expectVal]=getCisType(diameter,numberOfNodes,cisNumberOrder,numIS,cons,thr,verbose,minStatNode)
	CISProb=getCisProb(diameter,numberOfNodes,cisNumberOrder,numIS,cons,thr)
	if (noLowOrder) and numberOfNodes<6:
		CISvalue="Low_Order"
	getcontext().prec = 10
	CISProb=Decimal(str(CISProb))
	if prob!="" and float(prob)>float(CISProb):
		#print float(prob),CISProb
		outputFileColorTable.write(str(comp)+","+str(numberOfNodesApp)+","+str(chromo)+","+str(clusterPosition)+","+str(diameter)+","+nome+","+str(S)+","+CISvalue+","+str(CISProb)+","+valoriStr+"\n")
	elif prob=="":
		outputFileColorTable.write(str(comp)+","+str(numberOfNodesApp)+","+str(chromo)+","+str(clusterPosition)+","+str(diameter)+","+nome+","+str(S)+","+CISvalue+","+str(CISProb)+","+valoriStr+"\n")
		
		
	return CISProb,CISvalue,expectVal,S,numberOfNodesApp

def printAttributes(comp,chromo,G):
	for n in nx.nodes_iter(G):
		intensity=G.node[n]['intensity']
		gene=G.node[n]['gene']
		print ("\tcluster->%i,%s" % (comp,intensity))
	return

def scanLine(lista,start):
	patList=""
	primo=1
	for pa in range(start,len(lista)):
		if primo:
			patList=lista[pa]
			primo=0
		else:
			patList=patList+","+lista[pa]
	return patList

def getRegion(lista):
	for i in range(1,25):
		if long(lista[1])-1000000000<0:
			return i
		lista[1]=lista[1]-1000000000
	print "Maybe your genome is not human "+str(lista[1])
	return 

def isISinPrioriRange(ISAnnot,noIS):
	chromosome=str(ISAnnot[0])
	pos=ISAnnot[1]
	if not noIS.has_key(chromosome):
		return False
	for p in noIS[chromosome].keys():
		if (int(p)-noIS[chromosome][p]<=int(pos))and(int(p)+noIS[chromosome][p]>=int(pos)):
			return True
	return False


def filterComponent(G,maxCluster,noLowOrder,verbose):
	#ypos=nx.graphviz_layout(g)
	H=nx.Graph()
	lista=G.nodes()
	lista.sort()
	minmax=[0,10000000000000000]
	numNodes=0
	if verbose:
        	print "---------------------------------------------"
		print "comp has "+str(nx.number_of_nodes(G))+"nodes"
	#if nx.number_of_nodes(G)==1:
		
	#	return H
	for i in range(0,len(lista)-3):
	    for j in range(i+5,len(lista)):
		#if (lista[j]-lista[i]<=minmax[1]-minmax[0]) and (lista[j]-lista[i] < maxCluster) and (j-i>numNodes): 
		if (lista[j]-lista[i] < maxCluster) and (j-i>numNodes): 
			minmax=[lista[i],lista[j]]
			numNodes=j-i
	  		if verbose:
				print str(lista[i])+"->"+str(lista[j])
		elif (lista[j]-lista[i] < maxCluster) and (j-i==numNodes):
			if lista[j]-lista[i]<minmax[1]-minmax[0]:
				minmax=[lista[i],lista[j]] 
				if verbose:
					print str(lista[i])+"->"+str(lista[j])
		elif (lista[j]-lista[i] > maxCluster):
			break
	lista=G.nodes()
        lista.sort()	
	
	#if nx.number_of_nodes(G)==1:
	#	return H
	if not noLowOrder:
	    return G
	else: 
		for node in lista:
                	if node<minmax[0] or node>minmax[1]:
				if verbose:
                       	 		print "remove node ->"+str(node)
                        	G.remove_node(node)
		if verbose:
            		print "comp has "+str(nx.number_of_nodes(G))+"nodes"
            	return G

def getCisNumber(G,cisNumberOrder):
	for i in range(2,9):
		for gg in G:
        		if len(gg.nodes())==i:
				cisNumberOrder[i]=cisNumberOrder[i]+1
	for gg in G:
                        if len(gg.nodes())>8:
                                cisNumberOrder[8]=cisNumberOrder[8]+1
	return cisNumberOrder




parser = optparse.OptionParser()
parser.add_option( '-i', '--in', dest='input',help='The set to be translated' )
parser.add_option( '-o', '--out', dest='output', default="filtered.csv",help='The output table [filtered.csv]' )
parser.add_option( '-t', '--thr', dest='thr', default="50000", help='Threshold for distance [50]Kbp' )
parser.add_option( '-l', '--len', dest='maxCluster', default="1000000", help='Threshold for maximal CIS diameter [1]Mbp' )
parser.add_option( '-L', '--noLow', action='store_true', dest='noLowOrder', default=False, help='Filter Low Order CIS with thresholds [False]' )
parser.add_option( '-m', '--minNode', dest='minStatNode', default="3", help='Minimum value for compute CIS p-value and MAP classifier [3]' )
parser.add_option( '-s', '--stat', action='store_true', dest='statistical', default=False, help='enter in statistical modality [False]' )
parser.add_option( '-r', '--order', dest='order', default="0", help='min order on graph [0]' )
parser.add_option( '-R', '--orderCis', dest='orderCis', default="0", help='min order of returned CIS [0]' )
parser.add_option( '-b', '--bigger', action='store_true', dest='bigger', default=False, help='also bigger order on graph [False]' )
parser.add_option( '-c', '--cons', action='store_true', dest='cons', default=False, help='Use max values in p-value computation. More conservative [False]' )
parser.add_option( '-g', '--graph', dest='graphOutputFile', default="", help='The output graph files prefix [None]' )
parser.add_option( '-e', '--export', dest='extFile', default="", help='The external graph file in csv format (gml format)' )
parser.add_option( '-p', '--prob', dest='prob', default="1.1", help='Return only CIS with p < prob [1]' )
parser.add_option( '-a', '--annot', dest='annot', default="", help='Do not connect nodes where first nodes has annotation1=annot (ie IS and genes). IS at the head of the file otherwise no connection [None]' )
parser.add_option( '-A', '--annot2', action='store_true',dest='annot2', default=False, help='Do not connect nodes where first node has annotation2=annotation4 (ie elements on different dna strands) [False]' )
parser.add_option( '-H', '--entropy', action='store_true',dest='entropy', default=False, help='Compute normalized entropy on node composition. Species are determined by names in annotation field 3. The value in this field must start with an alphabetic character ([A-Za-z]) [False]' )
parser.add_option( '-v', '--verbose', action='store_true',dest='verbose', default=False, help='Return lot of info about IS and CIS [False]' )
parser.add_option( '-P', '--priori', dest='priori', default="", help='Introduce apriori file on genomic positions that were not considered for CIS formation [""]' )
parser.add_option( '-T', '--thresholdPriori', dest='thresholdPriori', default="50000", help='radius of influence on IS of the forbidden apriori regions [50]Kbp' )
parser.add_option( '-d', '--prioriDiameter', action='store_true',dest='prioriDiameter', default=False, help='The radius of influence on IS is in the 3th column into the apriori file [False]' )
parser.add_option( '-n', '--chrNumber', dest='chrNumber', default="24", help='Number of chromosomes into the specie (Homo s. [24])' )
parser.add_option( '-B', '--binWidth', dest='binWidth', default="1000000", help='Dimension of the window in which the theoretical CIS expression is binned [1]Mbp' )
parser.add_option( '-S', '--separator', dest='separator', default=",", help='Field separator for the input file. Tabs are automatically converted into separator value so you don\'t need to set it [,]' )
parser.add_option( '-O', '--separatorPriori', dest='separatorPriori', default=",", help='Field separator for the priori file. Tabs are automatically converted into separator value so you don\'t need to set it [,]' )


(options, args) = parser.parse_args()
try:
        inputFile = open (options.input, "rU")
except:
                raise 'Errore in lettura','File'
inputFile.close()
if not options.statistical:
	outputFileColor = open("Color_"+options.output,"w")
	outputFile = open (options.output, "w")
	outputFileExp= open("Exp"+options.output,"w")
else:
	outputFileColor = open(options.output,"w")
	outputFileExp=0;
   	outputFile=0	
	externalFile=0
thr=long(options.thr)
#IS=[]
IS={}
if options.extFile!="":
	E=nx.Graph()
print "extFile="+options.extFile
expDict={}
noIS={}
maxVal={}
diameter=00
font=5
region=550
binWidth=int(options.binWidth)
chrNumber=int(options.chrNumber)
graphs=[]
separator=options.separator
separatorPriori=options.separatorPriori
annotation1=""
annotation2=""
annotation3=""
annotation4=""
annotation5=""
annotation6=""
start=3
species={}
specieNumber=0
removedAprioriIS=0
regular=re.compile('[A-Za-z_]+')
numIS=0
# if known hotspot regions are known and you want to avoid CIS in that regions 
# you have to construct a map that mask some genome regions
# noIS is a double dictionary that has the structure noIS[chr][pos]
if (options.priori!=""):
	for linePriori in fileinput.input(options.priori):
		noVirgolettePriori=linePriori.replace('\"','')
       		noVirgolettePriori=noVirgolettePriori.replace('\t',separatorPriori)
		value=noVirgolettePriori.strip('\n')
		indices=value.split(separatorPriori)
		if not noIS.has_key(indices[0]):
			noIS[indices[0]]={}
		if options.prioriDiameter:
			noIS[indices[0]][indices[1]]=int(indices[2])
			
		else:
			noIS[indices[0]][indices[1]]=int(options.thresholdPriori)
		
		#print indices[0]+" "+indices[1]+" "+str(noIS[indices[0]][indices[1]])
	
for i in range(0,chrNumber+1):
	IS[str(i+1)]=[]
	IS[str(i+1)].append(str(i+1)+separator+str(0)+separator+"")
#for chrNum in range(0,chrNumber):
#	print str(chrNum)+" "+str(chrNum+1)+" "+str(len(IS[str(chrNum+1)]))
#	for i in range(len(IS[str(chrNum+1)])-1):
#		print str(i)

for lineFirst in fileinput.input(options.input):
	#IS.append(lineFirst.split()[0])
	noVirgolette=lineFirst.replace('\"','')	
	noVirgolette=noVirgolette.replace('\t',separator)	
	
#Is the apriori constrain set on?
	# No, append the IS to the list
	
	noVirgolette=noVirgolette.strip('\n')
        ISAnnot=noVirgolette.split(separator)
	ISAnnot[0]=ISAnnot[0].replace('chr','')
	if (ISAnnot[0]=="X") or (ISAnnot[0]=="x"):
        	ISAnnot[0]=str(chrNumber-1)
        elif (ISAnnot[0]=="Y") or (ISAnnot[0]=="y"):
        	ISAnnot[0]=str(chrNumber)
        elif (ISAnnot[0]=="M") or (ISAnnot[0]=="m"):
        	ISAnnot[0]=str(chrNumber+1)
	if (options.priori==""):
		if not IS.has_key(ISAnnot[0]):
			IS[ISAnnot[0]]=[]
		IS[ISAnnot[0]].append(noVirgolette.strip('\n'))
		numIS=numIS+1
		#print str(ISAnnot[0])+" "+str(IS[ISAnnot[0]])
		if(options.verbose):
			print noVirgolette.strip('\n')
        # Yes, check if the IS falls into a forbidden region
	else:
		if isISinPrioriRange(ISAnnot,noIS):
			removedAprioriIS=removedAprioriIS+1
			if options.verbose:
				print "Removed IS: "+noVirgolette
			continue
		else:
			if not IS.has_key(ISAnnot[0]):
				IS[ISAnnot[0]]=[]
			#print str(ISAnnot[0])+" "+str(IS[ISAnnot[0]])
			IS[ISAnnot[0]].append(noVirgolette.strip('\n'))
			numIS=numIS+1
			if(options.verbose):
				print noVirgolette

if (options.priori!=""):
	print "\nNumber of apriori filtered IS: "+str(removedAprioriIS)

# construc a void graph for each chromosome
for chrNum in range(0,chrNumber+1):
    graphs.append(nx.Graph())
# Populate the graphs
    for i in range(len(IS[str(chrNum+1)])-1):
		
		ISAnnot=IS[str(chrNum+1)][i].split(separator)
		ISAnnot[0]=ISAnnot[0].replace('chr','')
		
		if (ISAnnot[0]=="X") or (ISAnnot[0]=="x"):
			ISAnnot[0]=str(chrNumber-1)
		elif (ISAnnot[0]=="Y") or (ISAnnot[0]=="y"):
                        ISAnnot[0]=str(chrNumber)
        	elif (ISAnnot[0]=="M") or (ISAnnot[0]=="m"):
        		ISAnnot[0]=str(chrNumber+1)
		chromo1=int(ISAnnot[0])-1
		IS1pos=long(ISAnnot[1])
		#IS1=long(ISAnnot[1])+1000000000*(chromo1)
		IS1=IS1pos+1000000000*(chromo1)
		#print IS[i].split()[2]
		#codice sporco x selezionare due eventuali campi di annotazione
		if len(ISAnnot)>3:
			if(regular.match(ISAnnot[3])!=None):
				annotation1=ISAnnot[3]
				start=4
		if len(ISAnnot)>4:
			if(regular.match(ISAnnot[4])!=None) or (ISAnnot[4]=="-") or (ISAnnot[4]=="+"):
				annotation2=str(ISAnnot[4])
				start=5
		if len(ISAnnot)>5:
			if(regular.match(ISAnnot[5])!=None):
				annotation3=ISAnnot[5]
				if not species.has_key(annotation3):	
					species[annotation3]=1
					specieNumber=specieNumber+1
				start=6
		nameIS1=ISAnnot[2]
		if (options.annot!="") and (options.annot==annotation1):
			continue
		patientList1=scanLine(ISAnnot,start)
		for j in range(i+1,len(IS[str(chrNum+1)])):
			ISAnnot=IS[str(chrNum+1)][j].split(separator)
			ISAnnot[0]=ISAnnot[0].replace('chr','')
			if (ISAnnot[0]=="X") or (ISAnnot[0]=="x"):
				ISAnnot[0]=str(chrNumber-1)
			elif (ISAnnot[0]=="Y") or (ISAnnot[0]=="y"):
                        	ISAnnot[0]=str(chrNumber)
        		elif (ISAnnot[0]=="M") or (ISAnnot[0]=="m"):
        			ISAnnot[0]=str(chrNumber+1)
			chromo2=int(ISAnnot[0])-1
			IS2pos=long(ISAnnot[1])
			IS2=IS2pos+1000000000*(chromo2)
			#IS2=long(ISAnnot[1])+1000000000*(chromo2)
			nameIS2=ISAnnot[2]
			if len(ISAnnot)>3:
				if(regular.match(ISAnnot[3])!=None):
					start=4
					annotation4=ISAnnot[3]
			if len(ISAnnot)>4:
				if(regular.match(ISAnnot[4])!=None) or (ISAnnot[4]=="-") or (ISAnnot[4]=="+"):
					start=5
					annotation5=str(ISAnnot[4])
			if (options.annot2) and (annotation5!=annotation2):
				#print annotation5+"!="+annotation2
				continue
			#else: print annotation5+"="+annotation2
			if len(ISAnnot)>5:
				if(regular.match(ISAnnot[5])!=None):
					annotation6=ISAnnot[5]
					if not species.has_key(annotation6):	
						species[annotation6]=1
						specieNumber=specieNumber+1
					start=6
			patientList2=scanLine(ISAnnot,start)
			dist=abs(IS1-IS2);
			if (dist<thr):
				sys.stdout.write("\rProcessing chromosome -> %i  " % int(chromo1+1))
				graphs[chromo1].add_node(IS1,gene=nameIS1,ANNOTATION1=annotation1,ANNOTATION2=annotation2,ANNOTATION3=annotation3,chromosome=str(chromo1+1),position=IS1pos,intensity=patientList1)
				graphs[chromo2].add_node(IS2,gene=nameIS2,ANNOTATION1=annotation4,ANNOTATION2=annotation5,ANNOTATION3=annotation6,chromosome=str(chromo2+1),position=IS2pos,intensity=patientList2)
				#G.add_edge(long(IS[i]),long(IS[j]),weight=thr-dist)
				graphs[chromo1].add_edge(IS1,IS2,weight=dist)
				#print IS[i]+" "+IS[j]+"->"+str(dist)

cisNumberOrder=[0,0,0,0,0,0,0,0,0,0]
for l in range(0,chrNumber):
    C=nx.connected_component_subgraphs(graphs[l])
    cisNumberOrder=getCisNumber(C,cisNumberOrder)

for l in range(0,chrNumber):
    C=nx.connected_component_subgraphs(graphs[l])
    expDict[l]={}
    comp=0
    maxVal[l]=0

    for gg in C:
	if len(gg.nodes())<=int(options.orderCis):
		continue
	g=filterComponent(gg,int(options.maxCluster),options.noLowOrder,options.verbose)
	if len(g.nodes())==0: 
		continue
	# each CIS contibute to the expression in the bin with his order
	# The expression is computed in a BIN that has a width equal to binWidth 
	# Take into account the paramenter t (threshold) for graph construction
	listNodes=g.nodes()
        listNodes.sort()
	# mavVal is a dictionary that contains the position of the last CIS found in the chromosome l
        if maxVal[l]<int((listNodes[len(listNodes)-1]-1000000000*l)/binWidth)+1:
		maxVal[l]=int((listNodes[len(listNodes)-1]-1000000000*l)/binWidth)+1
	compPos=int((listNodes[len(listNodes)-1]+listNodes[0])/2)
        position=int((compPos-1000000000*l)/binWidth)+1
	if not expDict[l].has_key(position):
		#print position
		expDict[l][position]=0
	expDict[l][position]=expDict[l][position]+nx.number_of_nodes(g)
	clusterPosition=int(compPos-1000000000*l)+1
	diameterCluster=int(listNodes[len(listNodes)-1]-listNodes[0])
	order=nx.number_of_nodes(g)+g.number_of_selfloops()
	if (order==1): order=2
	[CISProb,CISValue,expectVal,S,newOrder]=saveAttributes(comp+1,order,l+1,g,clusterPosition,diameterCluster,outputFileColor,outputFile,options.noLowOrder,options.statistical,cisNumberOrder,numIS,options.cons,options.prob,options.thr,options.entropy,options.verbose,specieNumber,int(options.minStatNode),options.annot)
	if CISProb>1: CISProb=1
	if(options.extFile!="")and(not(options.statistical)):
	    if (options.prob!="")and(float(options.prob)>=float(CISProb)): #and float(order>options.order):
		[E.add_node(str(u),edata,entropy=str(S),CISProb=str(CISProb),MAPflag=CISValue,DIAMETER=diameterCluster,CISPosition=clusterPosition,NodeNumber=newOrder,CISexpect=expectVal) for u,edata in g.nodes(data=True) if 'gene' in edata ]
		[E.add_edge(str(u),str(v),edata,prob=str(CISProb),MAP=CISValue) for u,v,edata in g.edges(data=True) if 'weight' in edata ]
	comp=comp+1	
if not options.statistical:
    for i in range(0,chrNumber):
	for j in range(1,250):
		if not expDict[i].has_key(j):
			outputFileExp.write(str(i+1)+"-"+str(j)+",0,"+"\n")
		else:
			outputFileExp.write(str(i+1)+"-"+str(j)+","+str(expDict[i][j])+","+"\n")
    outputFile.close()
    outputFileExp.close()
    if (options.extFile!=""):
	nx.write_gml(E,options.extFile)
outputFileColor.close()
print "\n"
