
"""
TimeTP

Created by Kyuri Jo on 2016-01-24.
Copyright (c) 2016 Kyuri Jo. All rights reserved.

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

# Microarray
def makeRcode(dataPath, resPath, gFile, numTP, single) :
	groupF= open(os.path.join(dataPath, gFile), 'r')
	Rcode = open(os.path.join(dataPath, 'limma.R'), 'w')
	g = {}
	for l in groupF.readlines():
		tp = l.split()
		n = int(tp[0])
		if n not in g :
			g[n]=set()
		g[n].add(tp[1].rstrip())

	Rcode.write('library(limma)\nlibrary(affy)\n')
	Rcode.write('setwd("'+dataPath+'")\n')
	i=0
	while (i<numTP) :
		if single==0 :
			n1 = len(g[i])
			n2 = len(g[i+numTP])
		else :
			n1 = len(g[0])
			n2 = len(g[i+1])
		s1 = []
		s2 = []
		for j in range(n1):
			s1.append('C')
		for j in range(n2):
			s2.append('T')	
		if single==0 :	
			Rcode.write('Data <- ReadAffy("'+'","'.join(list(g[i]))+'","'+'","'.join(list(g[i+numTP]))+'")\n')
		else :
			Rcode.write('Data <- ReadAffy("'+'","'.join(list(g[0]))+'","'+'","'.join(list(g[i+1]))+'")\n')
		Rcode.write('eset <- rma(Data)\nwrite.exprs(eset, file="Expr_T'+str(i+1)+'.txt")\npData(eset)\ndesign <- model.matrix(~factor(c("'+'","'.join(s1)+'","'+'","'.join(s2)+'")))\ncolnames(design) <- c("C","T")\n')
		Rcode.write('fit <- lmFit(eset, design)\nfit <- eBayes(fit)\noptions(digits=2)\ntable <- topTable(fit, coef=2, adjust="fdr", n=Inf)\n')
		
		Rcode.write('write.table(table, file="'+resPath+'/DEG_result.T'+str(i+1)+'.txt", sep="\\t", col.names=TRUE, row.names=TRUE)\n')
		i+=1

# RNA-seq
def makeRcode2(dataPath, resPath, gFile, countF, numTP, single, read) :
        groupF= open(os.path.join(dataPath, gFile), 'r')
        Rcode = open(os.path.join(dataPath, 'deseq2.R'), 'w')
        g = {}
        for l in groupF.readlines():
                tp = l.split()
                n = int(tp[0])
                if n not in g :
                        g[n]=[]
                g[n].append(tp[1].rstrip())

	Rcode.write('library("DESeq2")\n')
	Rcode.write('setwd("'+dataPath+'")\n')
	i=0
	while (i<numTP) :
		if single==0 :
                        c1 = g[i]
                        c2 = g[i+numTP]
                else :
                        c1 = g[0]
                        c2 = g[i+1]
		rows = []
		tp = []
		for j in range(len(c2)) :
                        rows.append('"'+c2[j]+'"')
                        tp.append('"treated"')
                        tp.append('"'+read+'"')

		for j in range(len(c1)) :
			rows.append('"'+c1[j]+'"')
			tp.append('"untreated"')
			tp.append('"'+read+'"')
		Rcode.write('colData <- matrix(c('+','.join(tp)+'), nrow='+str(len(c1)+len(c2))+', ncol=2, byrow=TRUE)\n')
		Rcode.write('rownames(colData) <- c('+','.join(rows)+') \n')
		Rcode.write('colnames(colData) <- c("condition", "type")\n')
		Rcode.write('colData <- data.frame(colData)\n')
		Rcode.write('countData <- read.csv("'+countF+'", header=TRUE, row.names=1)[, c('+','.join(rows)+')]\n')
		Rcode.write('dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)\n')
		Rcode.write('dds$condition <- relevel(dds$condition, "untreated")\n')
		Rcode.write('dds <- DESeq(dds)\nres <- results(dds)\n')
		Rcode.write('write.table(res, file="'+resPath+'/DEG_result.T'+str(i+1)+'.txt", sep="\\t", col.names=TRUE, row.names=TRUE)\n')
		i+=1
	Rcode.close()

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
	
def makeProfile(dataPath, resPath, geneF, countF, numTP, remove0, fdr, single, thre, profile) :
	gpDict = {}
	if geneF=='NA' :
		gFile = open(dataPath+'/'+countF, 'r')
		for l in gFile.readlines() :
			tp = l.split(',')
			gene = tp[0].replace('\"', '')
			gpDict[gene]=gene
	else :
		gFile = open(dataPath+'/'+geneF, 'r')
		gFile.readline()
		for l in gFile.readlines() :
			tp = l.rstrip().split()
			gene = tp[1]
			probe = tp[0].replace('\"','')
			if gene not in gpDict :
				gpDict[gene] = probe

	pgDict = {y:x for x,y in gpDict.iteritems() }

	outfile = open(resPath+'/'+profile+'_profile.txt','w')
	count = []
	for i in range(numTP) :
		count.append(0)

	if type=='RNA-seq' :	#DESeq2
		foldInd = 2
		pInd = 5
	else :			#limma
		foldInd = 1
		pInd = 4	

	T0 = set()
	prof = {}
	for i in range(numTP) :
		tfile = open(resPath+'/DEG_result.T'+str(i+1)+'.txt','r')
		tfile.readline()
		for l in tfile.readlines() :
			tp = l.split()
			if len(tp)>7 :
				tp = tp[1:]
			probe =  tp[0].replace("\"", "")
			if probe not in pgDict :
				continue
			gene = pgDict[probe]

			if gene not in prof :
				prof[gene]=[]
				for j in range(numTP) :
					prof[gene].append(0)
				
			if (remove0==1 and i>0 and gene in T0) :
				continue

			# fold profile : no need to check p-value
			if profile=='fold':
				if isNum(tp[foldInd]) :
					prof[gene][i]=float(tp[foldInd])
			# DEG profile (-1, 0, 1)
			elif profile=='DEG' :
				if (isNum(tp[pInd+fdr]) and float(tp[pInd+fdr])<thre) :
					if remove0==1 and i==0 :
						T0.add(gene)
						continue
					if (float(tp[foldInd])<0) :
						prof[gene][i]=-1	# mean(g2)<mean(g1)
					else :
						prof[gene][i]=1
					count[i]+=1	
		tfile.close()
	for gene,ps in prof.iteritems() :
		outfile.write(gene+'\t'+'\t'.join(map(str, ps))+'\n')

if __name__ == "__main__" :
	conf = open(sys.argv[1], 'r')
	for l in conf.readlines() :
        	tp = l.rstrip().split()
		if tp[0]=='dataDir' :
			dataPath = tp[1]
		elif tp[0]=='outDir' :
			resPath = tp[1]
		elif tp[0]=='geneConvFile' :
			geneF = tp[1]
		elif tp[0]=='countFile' :
			countF = tp[1]
		elif tp[0]=='numTP' :
			numTP = int(tp[1])
		elif tp[0]=='removeT0' :
			remove0 = int(tp[1])
		elif tp[0]=='FDR' :
			fdr = int(tp[1])
		elif tp[0]=='single' :
			single = int(tp[1])
		elif tp[0]=='type' :
			type = tp[1]
		elif tp[0]=='groupFile' :
			gFile = tp[1]
		elif tp[0]=='threshold' :
			thre = float(tp[1])
		elif tp[0]=='readType' :
			read = tp[1]
	if not os.path.exists(resPath) :
                os.mkdir(resPath)
	if single==1 :
		numTP-=1
	
	if sys.argv[2]=='0' :
		if type=='Microarray' :
			makeRcode(dataPath, resPath, gFile, numTP, single)
		else :
			makeRcode2(dataPath, resPath, gFile, countF, numTP, single, read)
		# Run R
	else :
		makeProfile(dataPath, resPath, geneF, countF, numTP, remove0, fdr, single, thre, 'DEG')

