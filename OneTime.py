
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
import scipy.stats as stats

import TRAP

colorCode = ["#FFAAAA", "#FF5555", "#FF0000", "#AAAAFF", "#5555FF", "#0000FF"]

def new_hypergeom_sf(k, *args, **kwds):
    (M, n, N) = args[0:3]
    try:
        return stats.hypergeom.sf(k, *args, **kwds)
    except Exception as inst:
        if k >= n and type(inst) == IndexError:
            return 0 ## or conversely 1 - hypergeom.cdf(k, *args, **kwds)
        else:
            raise inst

def calPF_one(g, wgene, redic, PFdic, recur) :
        if (g in redic) :
                PFsum = 0
                for alist in redic[g] :
                        if (alist[0] not in PFdic) :
				if (alist[0] in recur) :
					PFdic[alist[0]]=wgene[alist[0]]
				else :
					recur.add(g)
        	                        calPF_one(alist[0], wgene, redic, PFdic, recur)
                        PFsum = PFsum + alist[2]*(PFdic[alist[0]]/alist[1])
                PFdic[g]=PFsum+wgene[g]
        else :
                PFdic[g]=wgene[g]

def pickColor(fc, cut) :
	if (fc==0) :
		return "#FFFFFF"
	elif (fc>0) :
		index = int(fc/(cut/2))
		if (index >2) : 
			index =2
		return colorCode[index]
	else :	
		index = int(abs(fc)/(cut/2))+3
		if (index >5) :
			index =5
		return colorCode[index]

def pathwayAnalysis(outPath, fileN, wgene, wredic, DEG, DEGCut, idDic, pnameDic, ind) :
        tA = []
	status = []
	pORA = []
	pOFDR = []
	pPERT = []
	pG = []
	pFDR = []
	pMIX = []
        for i in range(0, fileN) :
	        tA.append(0)
        	status.append("")
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
		tempPF = {}
		recur = set()
                PFsum = 0
                for gene in wgene[i] :
                        calPF_one(gene, wgene[i], wredic[i], tempPF, recur)

                status[i] = sum(tempPF.values())
                currtA = sum(tempPF.values())-sum(wgene[i].values())
		tA[i] = currtA
                # Calculation of tA from H0
                nulltA = []
                repeat = 2000
                tempFC = copy.copy(wgene[i])
                sh = tempFC.values()
		recur = set()
                for j in range(0, repeat) :
                        randPF = {}
                        random.shuffle(sh)
                        for key, value in tempFC.iteritems() :
                                tempFC[key]=sh[random.randint(0, len(tempFC)-1)]
                        for g in tempFC :
                                calPF_one(g, tempFC, wredic[i], randPF, recur)
                        nulltA.append(sum(randPF.values())-sum(tempFC.values()))

                def above(x):
                        return round(x, 5)>=round(currtA, 5)
                def below(x):
                        return round(x, 5)<=round(currtA, 5)

                avgtA = np.median(nulltA)
                if (currtA >=avgtA) :
                        pPERT[i]=float(len(filter(above, nulltA)))/float(repeat)
                else :
                        pPERT[i]=float(len(filter(below, nulltA)))/float(repeat)
		if status[i]>=0 :
			status[i]="Activated"
		else :
			status[i]="Inhibited"

        # pORA
        genesum = {}
        DEGsum = set()
        for i in range(0, fileN) :
                genesum.update(wgene[i])
                DEGsum = DEGsum.union(DEG[i])
        totG = len(genesum)
        totD = len(DEGsum)

        for i in range(0, fileN) :
                pORA[i]=new_hypergeom_sf(len(DEG[i]), totG, totD, len(wgene[i]), loc=0)

        # pG
        for i in range(0, fileN) :
                c = pORA[i]*pPERT[i]
                if (c==0) :
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
	for gene in DEGsum :
		if (gene in idDic) :
			outDEG.write(idDic[gene][0]+"\n")
		else :
			outDEG.write(gene+"\n")
	outDEG.close()

	outColor = open(outPath+"_color.txt", "w")
	for g,fc in genesum.iteritems() :
		outColor.write(g+"\t"+pickColor(fc, DEGCut)+"\n")
	outColor.close()

        outPathway = open(outPath+"_pathway.txt", "w")
        outPathway.write("PathwayID\tPathwayName        \tGeneNum\tDEGNum\tpORA\tpORAfdr\ttA\tpPERT\tpG\tpG_FDR\tStatus\n")
        sortedkey = sorted(ind, key = lambda x : pMIX[ind[x]])
        for sk in sortedkey :
                i = ind[sk]
		pathwayName = ""
		if (sk in pnameDic) : 
			pathwayName = pnameDic[sk]
		nameLen = len(pathwayName)
		if (nameLen<15) :
			pathwayName = pathwayName+TRAP.addStr(18-nameLen)
		else : 
			pathwayName = pathwayName[0:15]+"..."
		if (wredic[i]=={}) :
			outPathway.write(sk+"\t"+pathwayName+"\t"+str(len(wgene[i]))+"\t"+str(len(DEG[i]))+"\t"+str(round(pORA[i],3))+"\t"+str(round(pOFDR[i], 3))+"\t.\t.\t.\t.\t.\n")
		else : 
	                outPathway.write(sk+"\t"+pathwayName+"\t"+str(len(wgene[i]))+"\t"+str(len(DEG[i]))+"\t"+str(round(pORA[i],3))+"\t"+str(round(pOFDR[i], 3))+"\t"+str(round(tA[i],3))+"\t"+str(round(pPERT[i],3))+"\t"+str(round(pG[i],3))+"\t"+str(round(pFDR[i],3))+"\t"+status[i]+"\n")
	outPathway.close()
