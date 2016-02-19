
"""
TRAP - Time-series RNA-seq Analysis Package

Created by Kyuri Jo on 2014-02-05.
Copyright (c) 2014 Kyuri Jo. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import math
import random
import copy
import numpy as np
import networkx as nx
import scipy.stats as stats
import matplotlib.pyplot as plt

import TRAP

def new_hypergeom_sf(k, *args, **kwds):
    (M, n, N) = args[0:3]
    try:
        return stats.hypergeom.sf(k, *args, **kwds)
    except Exception as inst:
        if k >= n and type(inst) == IndexError:
            return 0 ## or conversely 1 - hypergeom.cdf(k, *args, **kwds)
        else:
            raise inst
def calPF(g, gene, redic, PFdic, t, a, recur) :
        if (g in redic) :
                PFsum_pre = 0
                PFsum_curr = 0
                if (t>0) :
                        for asc in redic[g] :
                                PFsum_pre = PFsum_pre + asc[2]*(PFdic[t-1][asc[0]]/asc[1])
                for asc in redic[g] :
                        if (asc[0] not in PFdic[t]) :
				if (asc[0] in recur) :
					PFdic[t][asc[0]]=gene[asc[0]]
				else :
					recur.add(g)
	                                calPF(asc[0], gene, redic, PFdic, t, a, recur)
                        PFsum_curr = PFsum_curr + asc[2]*(PFdic[t][asc[0]]/asc[1])
                PFdic[t][g] = a*PFsum_pre + (1-a)*PFsum_curr + gene[g]
        else :
                PFdic[t][g] = gene[g]


def pathwayAnalysis(outPath, wgene, wredic, DEG, idDic, pnameDic, timeLag, timeLen, ind, fcList) : 
	fileN = len(ind)
        tA = []
        status = []
        pORA = []
        pOFDR = []
        pPERT = []
        pG = []
        pFDR = []
	pMIX = []
	totWgene = []
	for t in range(timeLen) :
		totWgene.append([])
		for g,exp in fcList.iteritems():
			totWgene[t].append(exp[t])
        for i in range(0, fileN) :
                tA.append(0)
                status.append([])
                pORA.append(0)
                pOFDR.append(0)
                pPERT.append(0)
                pG.append(0)
                pFDR.append(0)
		pMIX.append(0)
	
		if wredic[i]=={} :
			continue
	
		# pPERT
		# Calculation of PF
		tempPF = []
		currtA = 0
		recur = set()
		for t in range(0, timeLen) :
			tempPF.append({})
			for gene in wgene[t][i] :
				calPF(gene, wgene[t][i], wredic[i], tempPF, t, timeLag, recur)
			currtA = currtA + sum(tempPF[t].values())-sum(wgene[t][i].values())
			status[i].append(sum(tempPF[t].values()))
		tA[i] = currtA

		# Calculation of tA (Null dist)
		nulltA = []
		repeat = 2000
		for j in range(0, repeat) :
			nullTemp = 0
			randPF = []
			tempFC = []
			recur = set()
			for t in range(0, timeLen) :
				tempFCt = copy.copy(wgene[t][i])
#				sh = tempFCt.values()
#				random.shuffle(sh)
				for key, value in tempFCt.iteritems() :
#					tempFCt[key]=sh[random.randint(0, len(tempFCt)-1)]
					tempFCt[key]=totWgene[t][random.randint(0, len(totWgene[t])-1)]
				tempFC.append(tempFCt)
			for t in range(0, timeLen) :
				randPF.append({})
				for g in tempFCt :
					calPF(g, tempFC[t], wredic[i], randPF, t, timeLag, recur)
				nullTemp = nullTemp + sum(randPF[t].values())-sum(tempFC[t].values())
			nulltA.append(nullTemp)

		def above(x):
			return round(x, 5)>=round(currtA, 5)
		def below(x):
			return round(x, 5)<=round(currtA, 5)

		avgtA = np.median(nulltA)
		if (currtA >=avgtA) :
			pPERT[i]=float(len(filter(above, nulltA)))/float(repeat)
		else :
			pPERT[i]=float(len(filter(below, nulltA)))/float(repeat)
		for t in range(timeLen) :
			if status[i][t] >=0 :
				status[i][t]="Activated"
			else :
				status[i][t]="Inhibited"

	# pORA
	genesum = {}
	DEGset = []
	DEGsum = 0
	for i in range(0, fileN) :
		genesum.update(wgene[0][i])
	for t in range(0, timeLen) :
		DEGset.append(set())
		for i in range(0, fileN) :
			DEGset[t] = DEGset[t].union(DEG[t][i])
		DEGsum = DEGsum + len(DEGset[t])
	totG = len(genesum)*timeLen
	totD = DEGsum

	geneNum = []
	DEGnum = []
	for i in range(0, fileN) :
		geneNum.append(0)
		DEGnum.append(0)
		geneNum[i]=len(wgene[0][i])*timeLen
		for t in range(0, timeLen) :
			DEGnum[i] = DEGnum[i] + len(DEG[t][i])
	for i in range(0, fileN):
		pORA[i]=new_hypergeom_sf(DEGnum[i], totG, totD, geneNum[i], loc=0)

	# pG
	for i in range(0, fileN) :
		c = pORA[i]*pPERT[i]
		if (c<=0) :
			pG[i]==0
		else :
			pG[i] = c-c*math.log(c)

        pFDR = TRAP.cal_FDR(pG)
        pOFDR = TRAP.cal_FDR(pORA)

        for i in range(0, fileN) :
                if (wredic[i]=={}) :
                        pMIX[i]=pOFDR[i]
                else :
                        pMIX[i]=pFDR[i]

	# Text result
        outDEG = open(outPath+"_DEG.txt", "w")
	tempDEG = set()
	outDEG.write("GeneID\t")
	for t in range(0, timeLen) :
		outDEG.write(str(t)+"\t")
		tempDEG = tempDEG.union(DEGset[t])
	outDEG.write("\n")
	for gene in tempDEG :
       	        if (gene in idDic) :
               	        outDEG.write(idDic[gene][0]+"\t")
	        else :
       	                outDEG.write(gene+"\t")
		for t in range(0, timeLen):
			if (gene in DEGset[t]) :
				outDEG.write("O\t")
			else :
				outDEG.write("X\t")
		outDEG.write("\n")
        outDEG.close()

	outPathway = open(outPath+"_pathway.txt", "w")
	sthead = []
	for t in range(timeLen) :
		sthead.append('Status'+str(t+1))
	outPathway.write("PathwayID\tPathwayName        \tGeneNum\tDEGNum\tpORA\tpORAfdr\ttA\tpPERT\tpG\tpG_fdr\t"+'\t'.join(sthead)+"\n")
	sortedkey = sorted(ind, key = lambda x : pMIX[ind[x]])
	for sk in sortedkey :
		i = ind[sk]
		ststr = []
		for t in range(timeLen) :
			ststr.append('.')
		pathwayName = ""
		if (sk in pnameDic) :
			pathwayName = pnameDic[sk]

                nameLen = len(pathwayName)
                if (nameLen<15) :
                        pathwayName = pathwayName+TRAP.addStr(18-nameLen)
                else :
                        pathwayName = pathwayName[0:15]+"..."
		if (wredic[i]=={}) :
			outPathway.write(sk+"\t"+pathwayName+"\t"+str(geneNum[i])+"\t"+str(DEGnum[i])+"\t"+str(round(pORA[i],3))+"\t"+str(round(pOFDR[i], 3))+"\t.\t.\t.\t.\t"+'\t'.join(ststr)+"\n")
		else :
			outPathway.write(sk+"\t"+pathwayName+"\t"+str(geneNum[i])+"\t"+str(DEGnum[i])+"\t"+str(round(pORA[i],3))+"\t"+str(round(pOFDR[i],3))+"\t"+str(round(tA[i],3))+"\t"+str(round(pPERT[i],3))+"\t"+str(round(pG[i],3))+"\t"+str(round(pFDR[i],3))+"\t"+'\t'.join(status[i])+"\n")
	outPathway.close()


	# Graph result
	G = nx.Graph()
	for f,i in ind.iteritems() :
		if (wredic[i]=={}) :
			pval=pOFDR[i]
			if (pval<=0.01) :
				color='#B2FFB2'
			elif (pval<=0.05) :
				color='#4CFF4C'
			else :
				color='#FFFFFF'
		else :
			pval=pFDR[i]
			if (status[i]=="Activated") :
				if (pval<=0.01) :
					color='#FFB2B2'
				elif (pval<=0.05) :
					color='#FF4C4C'
				else :
					color='#FFFFFF'
			else :
				if (pval<=0.01) :
					color='#B2B2FF'
				elif (pval<=0.05) :
					color='#4C4CFF'
				else :
					color='#FFFFFF'
		if (len(wgene[0][i])>=300) :
			size = 4500
		elif (len(wgene[0][i])<=50) :
			size = 750
		else :
			size = len(wgene[0][i])*15
		if pval>=0.1 :
			continue
		G.add_node(f[:8], color=color, size=size, pval=pval)

	edgelist = []
	plist = nx.get_node_attributes(G, 'pval')
	for f1,i1 in ind.iteritems() :
		for f2,i2 in ind.iteritems() :
			if (f1!=f2 and f1[:8] in plist and f2[:8] in plist and (set([f1,f2]) not in edgelist)):
				edgelist.append(set([f1,f2]))
				inter = dict.fromkeys(x for x in wgene[0][i1] if x in wgene[0][i2])
				intsize = len(inter)
				if (intsize!=0) :
					G.add_edge(f1, f2, label=intsize)

	nodeN = G.number_of_nodes()
	if (nodeN>0) :
		graph_pos=nx.fruchterman_reingold_layout(G, dim=2, pos=None, fixed=None, iterations=50, weight='label', scale=1)
		nodes,ncolors = zip(*nx.get_node_attributes(G,'color').items())
		nodes,sizes = zip(*nx.get_node_attributes(G, 'size').items())
		colorDict = dict(zip(nodes, ncolors))
		sizeDict = dict(zip(nodes, sizes))
		if (G.number_of_edges()>0) :
			edges,labels = zip(*nx.get_edge_attributes(G, 'label').items())
			labeldic = dict(zip(edges, labels))
			plt.figure(figsize=(100,100))
			nx.draw(G, graph_pos, nodelist=nodes, edgelist=edges, node_color=ncolors, node_size=sizes)
			nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=labeldic)
			nx.draw_networkx_labels(G, graph_pos)
		else :
			plt.figure(figsize=(20,20))
			nx.draw(G, graph_pos, nodelist=nodes, node_color=ncolors, node_size=sizes)
			nx.draw_networkx_labels(G, graph_pos)
		plt.savefig(outPath+"_pathway.png")

		GfileN = open(outPath+"_node.txt", 'w')
		GfileE = open(outPath+"_edge.txt", 'w')
		GfileN.write('Node\tColor\tSize\n')
		for n in G.nodes() :
			GfileN.write(n+'\t'+colorDict[n]+'\t'+str(sizeDict[n])+'\n')
		GfileE.write('Node1\tNode2\tCommonGenes\n')
		for (a, b) in G.edges() :
			GfileE.write(a+'\t'+b+'\t'+str(labeldic[(a,b)])+'\n')
		GfileN.close()
		GfileE.close()
