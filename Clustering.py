
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

import re
import math
import itertools
import numpy as np
import scipy.stats as stats
import TRAP 


def ucdstring(fcvalue, pval, cut, pcut, diff) :
	if (diff==0 and fcvalue >= cut) or (diff==1 and pval < pcut and fcvalue >= cut) :
		return "U"
	elif (diff==0 and fcvalue <= cut*(-1)) or (diff==1 and pval < pcut and fcvalue <= cut*(-1)) :
		return "D"
	else :
		return "C"

def ucdlist(fclist, plist, cut, pcut, cuffdiff) :
	tempstr=""
	for i in range(len(fclist)) :
		if cuffdiff=='yes' :
			tempstr = tempstr+ucdstring(fclist[i], plist[i], cut, pcut, 1)
		else :
			tempstr = tempstr+ucdstring(fclist[i], 1, cut, pcut, 0)
	return tempstr

def clusteringAnalysis(outPath, wgene, fcList, plist, idDic, pnameDic, cut, pcut, timeLen, ind, cuffdiff):
	# Statistical test for each cluster
	ucddic = {} 	# ucddic = { string : index }
	oppodic = {}
	strings = ["U", "D", "C"]
	for i in range(pow(3, timeLen)) :
		num = i
		clstr = ""
		for t in range(timeLen-1, -1, -1) :
			clstr=strings[num/pow(3, t)]+clstr
			num = num%pow(3,t)	
		ucddic[clstr]=i
		oppodic[i]=clstr

        permList = list(itertools.permutations(np.arange(timeLen)))
	permN = len(permList)
	proN = int(math.pow(3, timeLen))

	pval = []	# pval[profN]
	FDR = []
	profiles = []	# profiles[profN][permN]
	original = []	
	clDic = {}      # clDic = { cluster : [ genes ] }

	for i in range(0, proN) :
		pval.append(0)
		FDR.append(0)
		profiles.append([])
		original.append(0)
		clDic[oppodic[i]]=[]
		for j in range(0, permN) :
			profiles[i].append(0)

	for gene,fc in fcList.iteritems() :
		if cuffdiff=='yes' :	
			proStr = ucdlist(fc, plist[gene], cut, pcut, cuffdiff)
		else :
			proStr = ucdlist(fc, [], cut, pcut, cuffdiff)
		clDic[proStr].append(gene)
		original[ucddic[proStr]] = original[ucddic[proStr]]+1
		for i in range(0, permN) :
			fcTemp = []
			pTemp = []
			for p in permList[i] :
				fcTemp.append(fc[p])
				if cuffdiff=='yes' :
					pTemp.append(plist[gene][p])
				else:
					pTemp.append(1)
			pInd = ucddic[ucdlist(fcTemp, pTemp, cut, pcut, cuffdiff)]
			profiles[pInd][i]=profiles[pInd][i]+1

	n = len(fcList)
	for i in range(0, proN):
		x = original[i]
		p = float(sum(profiles[i]))/(float(permN)*float(n))
		pval[i]=stats.binom.sf(x,n,p,loc=0)

	FDR=TRAP.cal_FDR(pval)

	output = open(outPath+"_cluster.txt", "w")
	output.write("Cluster\tgeneNum\tP-value\tFDR P-value\n")
	sortedList = [i[0] for i in sorted(enumerate(FDR), key=lambda x:x[1])]
	for j in sortedList :
		output.write(oppodic[j]+"\t"+str(original[j])+"\t"+str(round(pval[j],3))+"\t"+str(round(FDR[j],3))+"\n")
	output.close()

	# Pathway analysis for each cluster
	fileN = len(ind)
	ORAlist = []
	FDRlist = []
	kyolist = []
	output1 = open(outPath+"_gene.txt", "w")
	output2 = open(outPath+"_pathway.txt", "w")
	output2.write("Cluster\tSize\tPathwayID\tPathwayName\tGeneNum\tIntersection\tpORA\tpORA_FDR\n")

	def special_match(strg, search=re.compile(r'[^C]').search):
		return not bool(search (strg))

	g_path = []
	for j in range(0, fileN) :
		g_path.append([])
		for gene in wgene[0][j] : 
			if (gene in idDic) :
				for name in idDic[gene] :
					if name in fcList : 
						g_path[j].append(name)
						break
			else :
				g_path[j].append(gene)
	
	k=0
	for string,glist in clDic.iteritems() :
		for g in glist :
			output1.write(string+'\t'+g)
			for t in range(timeLen) :
				output1.write('\t'+str(fcList[g][t]))
			output1.write('\n')

		ORAlist.append([])	# ORAlist[proN][pathN]
		FDRlist.append([])
		kyolist.append([])
		for j in range(0, fileN) :
			kyolist[k].append(len([val for val in g_path[j] if val in glist]))
			ORAlist[k].append(stats.hypergeom.sf(kyolist[k][j], n, len(glist), len(wgene[0][j]), loc=0))
		FDRlist[k]=TRAP.cal_FDR(ORAlist[k])

		sortedInd = sorted(ind, key = lambda x : FDRlist[k][ind[x]])
		for si in sortedInd :
			s = ind[si]
			pathwayName = ""
			if (si in pnameDic):
				pathwayName = pnameDic[si]
	                nameLen = len(pathwayName)
        	        if (nameLen<15) :
                	        pathwayName = pathwayName+TRAP.addStr(15-nameLen)+"\t"
	                else :
        	                pathwayName = pathwayName[0:15]+"..."

			if (special_match(string)==False and FDRlist[k][s] and kyolist[k][s]!=0 and FDRlist[k][s]<=0.05) :
				output2.write(string+"\t"+str(len(glist))+"\t"+si+"\t"+pathwayName+"\t"+str(len(wgene[0][s]))+"\t"+str(kyolist[k][s])+"\t"+str(round(ORAlist[k][s],3))+"\t"+str(round(FDRlist[k][s],3))+"\n")
		k = k+1
	output1.close()
	output2.close()

