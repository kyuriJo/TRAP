
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

import sys
import os
import math
import numpy as np
from xml.dom.minidom import parseString

import OneTime as OT
import TimeSeries as TS
import Clustering as CL

def cal_FDR(arr) :
        plen = len(arr)
        p_arr = np.copy(arr)
        q_arr = [1]*plen
        sort_ind = np.argsort(p_arr)
        p_arr.sort()

        for i in range(plen-1, -1, -1) :
                if (i==plen-1) :
                        q_arr[i]=p_arr[i]
                else :
                        rank = i+1
                        q_arr[i]=min(p_arr[i]*(float(plen)/float(rank)), q_arr[i+1])
        q = [0]*plen
        for i in range(0, plen) :
                q[sort_ind[i]] = q_arr[i]
        return q

def cal_FWER(arr) :
        plen = len(arr)
        bon_arr = np.copy(arr)
        for i in range(0, plen) :
                bon_arr[i]=min(bon_arr[i]*float(plen), 1.0)
        return bon_arr

def addStr(num) :
	temp = ""
	while num>0 :
		temp = temp+" "
		num = num-1
	return temp 

def oriGene(gene) :
	tp = gene.split(".")
	return tp[0]

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def median(l) :
	if len(l)==2 :
		return sum(l)/2
	else :
		temp = sorted(l)
		return temp[len(l)/2]

def main() :
	# Start of TRAP
	sys.stdout.write("\n###### Pathway and clustering analysis ######\n\n")
	sys.stdout.flush()

	cuffPath = "cufflinks_result"
	diffPath = "cuffdiff_result"
	resultPath = "TRAP_result"

	controlList = []
	caseList = []
	diffList = []
	timeLen = 0
	geneIDPath = ""
	pnamePath = ""
	kgmlPath = ""
	xmlPath = ""
	cuffdiff= ""
	pCut = 0.05
	DEGCut = 2.0
	clusterCut = 2.0
	timeLag = 1.0

	fcList = {}	# fcList[geneID]=[fc_0, ... , fc_t]
	pVal = {} 	# pVal[geneID]=[p_0, ..., p_t]
	idDic = {}	# idDic[keggID]=geneID
	pnameDic = {}	# pnameDic[pID]=pathwayName
	 

	# Reading configuration file
	sys.stdout.write("Reading configuration file\t......\t")
	sys.stdout.flush()

	try :
		config = open("config.txt", "r")
		while True :
			cl = config.readline()
			if cl=="" :
				break
			tp = cl.split("=")
			if (len(tp)<2) :
				continue
			key=tp[0]
			val=tp[1]
			if (key[:7]=="control") :
				controlList.append(val.strip().split(','))	
			elif (key[:9]=="treatment") :
				caseList.append(val.strip().split(','))
			elif (key=="numTP") :
				timeLen=int(val.strip())	
			elif (key=="convfilePath") :
				geneIDPath = val.strip()
			elif (key=="pnamePath") :
				pnamePath = val.strip()
			elif (key=="kgmlPath") :
				kgmlPath = val.strip()
			elif (key=="cuffdiff") :
				cuffdiff = val.strip()
			elif (key=="pVal") :
				pCut = float(val.strip())
			elif (key[:4]=="diff") :
				diffList.append(val.strip())
			elif (key=="DEGCut") :
				DEGCut = float(val.strip())
			elif (key=="clusterCut") :
				clusterCut = float(val.strip())
			elif (key=="timeLag") :
				timeLag = float(val.strip())
			else :
				continue
		idFile = open(geneIDPath, "r")
		pnameFile = open(pnamePath, "r")
		xmlPath = os.walk(kgmlPath)

		if (cuffdiff=="no" and len(controlList)!=len(caseList)) :
			raise

	except IOError:
		print "Check if the configuration file exists"
	except :
		print "Configuration file error"
		raise

	# Reading ID-conversion / pathway name file
	for ids in idFile.readlines() :
		tp = ids.split("\t")
		tp2 = tp[1].split(";")
		tp3 = tp2[0].split(", ")
		if tp[0] in idDic :
			for name in tp3 :
				idDic[tp[0]].append(name.strip())
		else : 
			idDic[tp[0]]=[]
			for name in tp3 :
				idDic[tp[0]].append(name.strip())
	idFile.close()

	for path in pnameFile.readlines() :
		tp = path.split("\t")
		tp2 = tp[0].split(":")
		tp3 = tp[1].split(" - ")
		pnameDic[tp2[1]]=tp3[0]
	pnameFile.close()

	sys.stdout.write("Done\n")
	sys.stdout.flush()

	# Reading fpkm file
	sys.stdout.write("Reading expression files\t......\t")
	sys.stdout.flush()

	geneSum = set()
	if cuffdiff=="yes" :
                for j in range(timeLen) :
                        pfile = open(os.path.join(diffPath, diffList[j], "gene_exp.diff"), "r")
                        for l in pfile.readlines() :
                                tp = l.split()
				if not is_number(tp[9]) :
					continue
                                geneSum.add(tp[2])
                        pfile.close()
		for gene in geneSum :
			fcList[gene]=[]
                        pVal[gene]=[]
                for j in range(timeLen) :
                        pfile = open(os.path.join(diffPath, diffList[j], "gene_exp.diff"), "r")
                        temp = {}
                        temp2 = {}
                        for l in pfile.readlines() :
                                tp = l.split()
				if not is_number(tp[9]) :
					continue
				if ( tp[9]=='inf' or tp[9] == '-inf') :
					temp[tp[2]]=0
				else: 
	                                temp[tp[2]]=float(tp[9])
                                temp2[tp[2]]=float(tp[12])

                        for gene in geneSum :
                                if gene in temp :
                                        fcList[gene].append(temp[gene])
                                        pVal[gene].append(temp2[gene])
                                else :
                                        fcList[gene].append(0)
                                        pVal[gene].append(1)
                        pfile.close()

	else : 
                for j in range(timeLen) :
			for con in controlList[j] : 
	                        pfile = open(os.path.join(cuffPath, con, "genes.fpkm_tracking"), "r")
	                        for l in pfile.readlines() :
        	                        tp = l.split()
                	                if not is_number(tp[9]) :
                        	                continue
                                	geneSum.add(tp[4])
	                        pfile.close()
                        for case in caseList[j] :
                                pfile = open(os.path.join(cuffPath, case, "genes.fpkm_tracking"), "r")
                                for l in pfile.readlines() :
                                        tp = l.split()
                                        if not is_number(tp[9]) :
                                                continue
                                        geneSum.add(tp[4])
                                pfile.close()
                for gene in geneSum :
                        fcList[gene]=[]
		for j in range(timeLen) :
			temp1 = {}
			temp2 = {}
			for con in controlList[j] :
                                pfile = open(os.path.join(cuffPath, con, "genes.fpkm_tracking"), "r")
                                for l in pfile.readlines() :
					tp = l.split()
					if (tp[9]=="FPKM") :
						continue
					if tp[4] in temp1 : 
						temp1[tp[4]].append(float(tp[9]))
					else : 
						temp1[tp[4]]=[float(tp[9])]
				pfile.close()
			for case in caseList[j] :
                                pfile = open(os.path.join(cuffPath, case, "genes.fpkm_tracking"), "r")
				for l in pfile.readlines() :
					tp = l.split()
					if (tp[9]=="FPKM") :
						continue
					if tp[4] in temp2 :
                                                temp2[tp[4]].append(float(tp[9]))
                                        else :
                                                temp2[tp[4]]=[float(tp[9])]
				pfile.close()
			
			for gene in geneSum :
				med1 = 0
				med2 = 0
				if gene in temp1 and gene in temp2 :
					med1 = median(temp1[gene])
                                        med2 = median(temp2[gene])
				elif gene in temp1 :
                                        med1 = median(temp1[gene])
				elif gene in temp2 :
                                        med2 = median(temp2[gene])
				else :
					med1 = med1
					med2 = med2
				
				if (abs(med2-med1)<1.0) : 
					fcList[gene].append(0)
				else :
					fcList[gene].append(math.log((med2+0.01)/(med1+0.01),2))

	sys.stdout.write("Done\n")
	sys.stdout.flush()

	# Parsing xml file to get gene and relation information
	sys.stdout.write("Reading xml files\t\t......\t")
	sys.stdout.flush()

	i=0
	ind = {}
	DEG = []
	wgene = []
	wredic = []
	empty = []
	empty2 = []

	for t in range(0, timeLen) :
		wgene.append([])	#wgene[t][i]={keggID:fc}
		DEG.append([])		#DEG[t][i]=set(keggID)
		empty.append(0)
		empty2.append(1)

	for root,dirs,files in xmlPath:
	   for file in files:
		filetp = file.split(".")
		ind[filetp[0]]=i
		for j in range(0, timeLen) :
			wgene[j].append({})
			DEG[j].append(set())
		wredic.append({})	#wredic[i]={keggID:(list of [asc, length, j])}

		xmlfile = open(os.path.join(kgmlPath, file), "r")
		xmldata = xmlfile.read()
		dom = parseString(xmldata)
		xmlfile.close()

		geneSet = set()
		entrydic = {}
		entries = dom.getElementsByTagName("entry")
		for e in entries :
			if (e.attributes.getNamedItem("type").nodeValue == 'gene') :
				id = e.attributes.getNamedItem("id").nodeValue
				id = str(id)
				genes = e.attributes.getNamedItem("name").nodeValue
				genes = str(genes)
				genelist = genes.split()
				entrydic[id]=[]
				for g in genelist : 
					entrydic[id].append(g)
					geneSet.add(g)
			elif (e.attributes.getNamedItem("type").nodeValue == 'group') :
				id = e.attributes.getNamedItem("id").nodeValue
				id = str(id)
				comps = e.getElementsByTagName("component")
				entrydic[id]=[]
				for c in comps :
					geneId =c.attributes.getNamedItem("id").nodeValue
					for g in entrydic[geneId] :
						entrydic[id].append(g)
						geneSet.add(g)
		for g in geneSet :
			if (g in idDic) :
				nameExist = 0
				tpName = ""
				for name in idDic[g] :
					if name in fcList.keys() :
						nameExist = 1
						tpName = name
						break
				if nameExist==1 :
					for t in range(0, timeLen) :
						foldchange = fcList[tpName][t]
						wgene[t][i][g]=foldchange
						if (cuffdiff=="yes" and pVal[tpName][t]<pCut and abs(foldchange)>=DEGCut) :
							DEG[t][i].add(g)
						elif (cuffdiff=="no" and abs(foldchange)>=DEGCut) :
							DEG[t][i].add(g)
				else :
					for t in range(0, timeLen) :
						wgene[t][i][g]=0
						fcList[idDic[g][0]]=empty
						if (cuffdiff=="yes") :	
							pVal[idDic[g][0]]=empty2
			else :
				for t in range(0, timeLen) :
					wgene[t][i][g]=0
					fcList[g]=empty
					if (cuffdiff=="yes") :
						pVal[g]=empty2

		redic = wredic[i]
		relations = dom.getElementsByTagName("relation")
		for r in relations :
			subs = r.getElementsByTagName("subtype")
			ent1 = r.attributes.getNamedItem("entry1").nodeValue
			ent2 = r.attributes.getNamedItem("entry2").nodeValue
			if (not (subs==[])) :
				for s in subs :
					type = s.attributes.getNamedItem("name").nodeValue
					if (type=="activation" or type=="expression") :
						j=1
					elif (type=="inhibition" or type=="repression") :
						j=-1
					else :
						j=0
					if (j!=0 and (ent1!=ent2) and (ent1 in entrydic) and (ent2 in entrydic)) :
						for desc in entrydic[ent2] :
							length = len(entrydic[ent2])
							for asc in entrydic[ent1] :
								if (desc in redic) :
									redic[desc].append([asc, length, j])
								else :
									redic[desc]=[[asc, length, j]]
		i=i+1

	fileN = i
	sys.stdout.write("Done\n")
	sys.stdout.flush()


	# 1. One time point SPIA analysis
	sys.stdout.write("One time point SPIA analysis\n")
	sys.stdout.flush()

	for t in range(0, timeLen) :
		sys.stdout.write("\t"+str(t+1)+"th time point\t......\t")
		sys.stdout.flush()

		OT.pathwayAnalysis(os.path.join(resultPath, "OneTime_"+str(t+1)), fileN, wgene[t], wredic, DEG[t], DEGCut,  idDic, pnameDic, ind)

		sys.stdout.write("Done\n")
		sys.stdout.flush()


	# 2. Time-series SPIA analysis
	sys.stdout.write("Time-series SPIA analysis\t......\t")
	sys.stdout.flush()

	TS.pathwayAnalysis(os.path.join(resultPath, "TimeSeries"), wgene, wredic, DEG, idDic, pnameDic, timeLag, timeLen, ind, fcList)

	sys.stdout.write("Done\n")
	sys.stdout.flush()

	# 3. Clustering analysis
	sys.stdout.write("Clustering Analysis\t\t......\t")
	sys.stdout.flush()

	CL.clusteringAnalysis(os.path.join(resultPath, "Clustering"), wgene, fcList, pVal, idDic, pnameDic, clusterCut, pCut, timeLen, ind, cuffdiff)

	sys.stdout.write("Done\n")
	sys.stdout.flush()

if __name__ == "__main__" :
	main()

