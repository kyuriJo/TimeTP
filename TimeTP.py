
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

import copy
import re
import os
import shutil
import scipy.stats
import scipy.spatial
from scipy.signal import correlate 
from scipy.optimize import curve_fit
from scipy.misc import factorial
import networkx as nx
from xml.dom.minidom import parseString
import numpy as np
import sys
import math
import random

keggDict = {}
kgDict = {}
gkDict = {}

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def genDict(file1, KEGGgeneFilePath, numTP) :
	pfile = open(file1, 'r')

	# Gene ID - Profile
	# Gene ID include gene symbol or entrez gene ID
	pDict = {}
	for l in pfile.readlines() :
		tp = l.split()
		pDict[tp[0]]=map(float, tp[1:])

	# keggDict 	: Kegg ID - Profile
	# kgDict 	: Kegg ID - Gene ID
	temp = {}
	keggFile = open(KEGGgeneFilePath, 'r')
	
	nullP = []
	for i in range(numTP) :
		nullP.append(0)
	for l in keggFile.readlines() :
		tp = l.split("\t")
                tp2 = tp[1].split(";")
                tp3 = tp2[0].split(", ")
		species = tp[0].split(':')[0]
		ent = tp[0].split(':')[1]
		keggID = tp[0]
		firstName = ''

		tp3.append(ent)	
		# If keggID have gene expression
		for name in tp3 :
			if name==tp3[0] :
				firstName = name
			if name not in pDict :
				continue
			keggDict[keggID]=pDict[name]
			if name==ent :
				kgDict[keggID]=firstName
			else :
				kgDict[keggID]=name
			break

		# store gene ID is ith of keggID
		for i in range(len(tp3)) :
			if tp3[i] not in temp :
				temp[tp3[i]]=[]
			temp[tp3[i]].append([i, keggID])
	
		# If keggID doesn't have gene expression
		if keggID not in keggDict :
			keggDict[keggID]=nullP
			if firstName=='' :
				kgDict[keggID]='None'
			else :
				kgDict[keggID]=firstName

	# if a geneID is match to several keggID,
	# choose the keggID with minimum i
	for g,list in temp.iteritems() :
		min = 10000 
		for l in list :
			if l[0]<min :
				gkDict[g]=l[1]
				min = l[0]
	
	return species

def poisson_new(k, lamb):
        return (lamb**k/factorial(k)) * np.exp(-lamb)

# checkProfile fixed to use cross-correlation
def checkProfile(p1, p2, type, maxDelay) :
	if np.count_nonzero(p1)==0 or np.count_nonzero(p2)==0 or (type!=1 and type!=-1) :
                return False

	delay = -1
	p1 = np.array(p1)
	p2 = np.array(p2)
	re = correlate(p1, p2)
	if (type==1) :
		corr = max(re)
	else :
		corr = min(re)
	if corr==0 :
		return False
	tot = len(p1)-1
	temp = []
	for i in range(len(re)) :
		if (re[i]==corr) :
			temp.append(tot-i)
	if min(temp)>=0 :	# delay is positive
		delay = min([x for x in temp if x>=0])
	else :			# delay is negative
		delay =	max([x for x in temp if x<0]) 
	return (delay>0 and delay<=maxDelay)

# checkPath fixed to return network
def checkPathway(network, entry, keggDict, maxDelay, skip) :

	netlist_e = []	# set of entry for each path [ set(0,1,2), set(0,1,2), set(2,3,4) ... ]
	netlist_eg = []	# connected components for each path [ (0!g1->1!g2->2!g3), (0!g4->1!g5->2!g6), ... ]
	netlist_cc = []
	net_e = nx.DiGraph()
	net_eg = nx.DiGraph()

	for (a,b) in network.edges() :
		type = network[a][b]['type']
		e1 = entry[a]
		e2 = entry[b]
		for g1 in e1 :
			for g2 in e2 :
				if g1 not in keggDict or g2 not in keggDict :
					continue 
				if checkProfile(keggDict[g1], keggDict[g2], type, maxDelay) :
					net_e.add_edge(a, b)
					net_eg.add_edge(a+'!'+g1, b+'!'+g2)

	comps = nx.connected_components(net_eg.to_undirected())
	temp = []
	for c in comps :
		real_c = map(lambda x: x.split('!')[1], c)
		real_ce = set(map(lambda x: x.split('!')[0], c))
		# Check if the same gene set w/ different entry name exists
		if real_c not in temp and ((not skip) or len(real_ce)>2):
			temp.append(real_c)
			netlist_e.append(set(map(lambda x: x.split('!')[0], c)))
			netlist_eg.append(net_eg.subgraph(c))

	for eg in netlist_eg :
		cc_e = {}
		tempcc = 0
		for (a, b) in eg.edges() :
			e1 = a.split('!')[0]
			e2 = b.split('!')[0]
			g1 = a.split('!')[1]
			g2 = b.split('!')[1]
			t1 = np.array(keggDict[g1])
			t2 = np.array(keggDict[g2])
			if network[e1][e2]['type']==1 :
				tempcc = abs(max(correlate(t1, t2)))
			else :
				tempcc = abs(min(correlate(t1, t2)))
			if (e1, e2) not in cc_e :
				cc_e[(e1, e2)]=[]
			cc_e[(e1, e2)].append(tempcc)
		cc = 0
		for e, cclist in cc_e.iteritems() :
			cc+= sum(cclist)/float(len(cclist))
		netlist_cc.append(cc)			
 
	return netlist_e, netlist_eg, netlist_cc
 
# pathwayAnalysis fixed to use checkPathway and fit poisson
def pathwayAnalysis(numTP, networks, entries, DEGs, DEGPval, maxDelay, profile, threshold) :
	
        labeled = {}
	labeled_path = {}
        pathDict = {}   # [pathway] - pathlist, genepathlist, firstIDs
        pathPval = {}

        for pathway,net in networks.iteritems() :
		print ".. analyzing "+pathway
		pathPval[pathway] = []
		ent = entries[pathway]
		
		net_e, net_eg, net_cc = checkPathway(net, ent, keggDict, maxDelay, True)
		pathDict[pathway] = (net_e, net_eg)	

		if len(net_e)==0 :
			continue	
		# maximum number of node = number of entries in the pathway
		bin = np.array(range(len(ent)+1))
		bin2 = np.array(range(len(ent)*numTP))
		prob = np.zeros(len(ent)+1)
		prob2 = np.zeros(len(ent)*numTP)

		repeat = 1000
		#randDict = { key: keggDict[key] for key in DEGs[pathway][1] }
		randDict = copy.copy(keggDict)
		temp = []
		temp2 = []
		temp3 = []
                for j in range(repeat) :
                	values = randDict.values()
                        random.shuffle(values)
                        it = iter(values)
	                for key in randDict :
        	        	randDict[key]=next(it)

			temp_e, temp_eg, temp_cc = checkPathway(net, ent, randDict, maxDelay, False)

			for c in temp_e :
				prob[len(c)]+=1
				temp.append(len(c))
			for c in temp_cc :
				if profile =='DEG' :
					prob2[int(round(c))]+=1	# for poisson
				temp2.append(c)		# for hist

		# Fitting the Poisson dist (node)
		#prob = prob / np.sum(prob)
		#param, cov = curve_fit(poisson_new, bin, prob)

		# (cc)
		if profile=='DEG' :
			if (np.sum(prob2)==0) : # no subpath found in permutation
        	                for i in range(len(net_e)) :
	                                pathPval[pathway].append(0.0)
                	        continue
			prob2 = prob2 / np.sum(prob2)
			param2, cov2 = curve_fit(poisson_new, bin2, prob2)

		for i in range(len(net_e)) :
			if profile == 'fold' :
				pval = len([x for x in temp2 if x>net_cc[i]])/float(len(temp2))        # cc and histogram
			else :
				pval = scipy.stats.poisson.sf(round(net_cc[i]), param2[0])	# cc and poisson
			pathPval[pathway].append(pval)
			print 'P-value:', pval
			# Fixed
			#if pval < threshold and DEGPval[pathway] < threshold :
			if pval < threshold : 
				source = set([x for x in net_eg[i].nodes() if net_eg[i].in_degree(x)==0])
				targets = set(map(lambda x: x.split('!')[1], net_eg[i].nodes()))
				print 'source '+str(len(source))
				for gene in source :
					realgene = gene.split('!')[1]
					if realgene not in labeled :
						labeled[realgene] = set()
					if pathway not in labeled_path :
						labeled_path[pathway] = set()
					labeled[realgene]|=targets
					labeled_path[pathway].add(kgDict[realgene])

	return pathPval, pathDict, labeled, labeled_path	

def parseXml(xmlPath, keggDict) :
	networks = {}
	entries = {}
	DEGs = {}	# pathName -> [DEGSet, geneSet]
        for root,dirs,files in os.walk(xmlPath):
           for file in files:
                pathName = file.split(".")[0]
                networks[pathName]=nx.DiGraph()
		
                xmlfile = open(os.path.join(xmlPath, file), "r")
                xmldata = xmlfile.read()
                dom = parseString(xmldata)
                xmlfile.close()

                geneSet = set()
                entries[pathName] = {}
		entrydic = entries[pathName]
		DEGs[pathName] = [set(),set()]

                domentry = dom.getElementsByTagName("entry")
                for e in domentry :
                        if (e.attributes.getNamedItem("type").nodeValue == 'gene') :
                                id = e.attributes.getNamedItem("id").nodeValue
                                id = str(id)
                                genes = e.attributes.getNamedItem("name").nodeValue
                                genes = str(genes)
                                genelist = genes.split()
                                entrydic[id]=set()
                                for g in genelist :
					DEGs[pathName][1].add(g)
					if g in keggDict and np.count_nonzero(keggDict[g])!=0 :
						DEGs[pathName][0].add(g)
                                        entrydic[id].add(g)
					
                        elif (e.attributes.getNamedItem("type").nodeValue == 'group') :
                                id = e.attributes.getNamedItem("id").nodeValue
                                id = str(id)
                                comps = e.getElementsByTagName("component")
                                entrydic[id]=set()
                                for c in comps :
                                        geneId =c.attributes.getNamedItem("id").nodeValue
                                        for g in entrydic[geneId] :
                                        	entrydic[id].add(g)

                relations = dom.getElementsByTagName("relation")
                for r in relations :
                        subs = r.getElementsByTagName("subtype")
                        ent1 = r.attributes.getNamedItem("entry1").nodeValue
                        ent2 = r.attributes.getNamedItem("entry2").nodeValue
			if (ent1==ent2) or (ent1 not in entrydic) or (ent2 not in entrydic) :
				continue
                        if (not (subs==[])) :
                                for s in subs :
                                        type = s.attributes.getNamedItem("name").nodeValue
					j = 0
                                        if (type=="activation" or type=="expression" or type=='indirect effect') and len(entrydic[ent1])!=0 and len(entrydic[ent2])!=0:
                                                j=1
                                        elif (type=="inhibition" or type=="repression") and len(entrydic[ent1])!=0 and len(entrydic[ent2])!=0 :
                                                j=-1
					elif (type=="compound" or type=="binding/association") and (len(entrydic[ent1])!=0 or len(entrydic[ent2])!=0) :
						j=2	# !!!
					if j==1 or j==-1 :
						networks[pathName].add_edge(ent1, ent2, type=j)
					elif j==2 :
						networks[pathName].add_edge(ent1, ent2, type=j)
						networks[pathName].add_edge(ent2, ent1, type=j)
	return (networks, entries, DEGs)

def new_hypergeom_sf(k, *args, **kwds):
    (M, n, N) = args[0:3]
    try:
        return scipy.stats.hypergeom.sf(k, *args, **kwds)
    except Exception as inst:
        if k >= n and type(inst) == IndexError:
            return 0 ## or conversely 1 - hypergeom.cdf(k, *args, **kwds)
        else:
            raise inst

def DEGAnalysis(DEGs) :
	totDEG = set()
	totGene = set()
	for path, l in DEGs.iteritems() :
		totDEG = totDEG | l[0]
		totGene = totGene | l[1]
	totDEGN = len(totDEG)
	totGeneN = len(totGene)

	DEGPval = {}
	for path, l in DEGs.iteritems() :
		A = len(l[0])
		B = totDEGN-A
		C = len(l[1])-A
		D = totGeneN-totDEGN-C
		odds, pval = scipy.stats.fisher_exact([[A,B],[C,D]])
		DEGPval[path]=pval
	return DEGPval


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def readFiles(file1, file2, species, network, labeled, maxDelay) :
	GRNfile = open(file1, 'r')
	PPIfile = open(file2, 'r')

	negProfit = 0

	TFset = set()
	for l in GRNfile.readlines() :
		tp = l.rstrip().split()
                TFset.add(tp[0])
		TFset.add(tp[0].replace("-", "").upper())
		if species+':'+tp[0] in kgDict :
			TFset.add(kgDict[species+':'+tp[0]])
	GRNfile.close()
	GRNfile = open(file1, 'r') 

	for l in GRNfile.readlines() :
		tp = l.rstrip().split()
		TF = tp[0]
		gene = tp[1]

		# add genes (found in the gkDict) to the network
		if TF not in gkDict and (TF.replace("-", "").upper() in gkDict) :
			TF = TF.replace("-", "").upper()
		# TF & gene should be in kegg
		if TF not in gkDict or gene not in gkDict :
			continue
		match = checkProfile(keggDict[gkDict[TF]], keggDict[gkDict[gene]], 1, maxDelay)
		match2 = checkProfile(keggDict[gkDict[TF]], keggDict[gkDict[gene]], -1, maxDelay)

		if match :
			type = 1
		elif match2 :
			type = -1
		else :
			continue

		label = len(labeled[gkDict[TF]]) if gkDict[TF] in labeled else negProfit
		TF = kgDict[gkDict[TF]]
		network.add_node(TF, label=label, isTF=1) 

		label = len(labeled[gkDict[gene]]) if gkDict[gene] in labeled else negProfit
		gene = kgDict[gkDict[gene]]
		isTF = 1 if gene in TFset else 0

		network.add_node(gene, label=label, isTF=isTF)
		network.add_edge(TF, gene, type=type, weight=1000)
	
	for l in PPIfile.readlines() :
		tp = l.rstrip().split()
		g1 = tp[0]
		g2 = tp[1]
		if g1 not in gkDict or g2 not in gkDict :
			continue
		
		# Change weight only for PPI
		rType = -1 if tp[3]=='inhibition' else 1
		match = checkProfile(keggDict[gkDict[g1]], keggDict[gkDict[g2]], rType, maxDelay)

		weight = int(tp[5]) if len(tp)>=6 and is_number(tp[5]) else int(tp[4])
                if not match :
			continue
                        #weight = weight*0.5
		tp1 = g1
		tp2 = g2
		# Change gene name for consistency with gene labels from expression data
		label = len(labeled[gkDict[g1]]) if gkDict[g1] in labeled else negProfit
		g1 = kgDict[gkDict[g1]]
		isTF = 1 if g1 in TFset else 0
		network.add_node(g1, label=label, isTF=isTF)

		label = len(labeled[gkDict[g2]]) if gkDict[g2] in labeled else negProfit
                g2 = kgDict[gkDict[g2]]
		isTF = 1 if g2 in TFset else 0
                network.add_node(g2, label=label, isTF=isTF)

		if tp[3]=='activation' :
			network.add_edge(g1, g2, type=1, weight=weight)
		elif tp[3]=='inhibition' :
			network.add_edge(g1, g2, type=-1, weight=weight)
		else : 
	                network.add_edge(g1, g2, type=1, weight=weight)
			network.add_edge(g2, g1, type=1, weight=weight)
		# edge type : -1/1 inhibition/activation, 0 TF-target(TF/gene), 2 else PPI
	
	return network

def influenceByNode(n, g, visited) :
	sList = g.successors(n)
	infSet = set([n])
	profit = g.node[n]['label']
	if len(sList)==0 :
		return infSet, profit
	
	for s in sList :
		if s in visited :
			continue
		visited.add(s)
		(a, b) = influenceByNode(s, g, visited)
		infSet = infSet | a
		profit += b
	return infSet, profit

def startInd(n) :
	p = keggDict[gkDict[n]]
	pSet = set()
	nSet = set()
	for i in range(len(p)) :
		if p[i]==1 :
			pSet.add(i)
		elif p[i]==-1 :
			nSet.add(i)
	return (pSet, nSet)

def nextInd(n, pair, type, maxDelay) :
	p = keggDict[gkDict[n]]
	if type==1 :
		prevP = pair[0]
		prevN = pair[1]
	elif type==-1 :
		prevP = pair[1]
		prevN = pair[0]
	pSet = set()
	nSet = set()
	for ind in prevP :
		for i in range(ind, min(ind+maxDelay+1,len(p))):
			if p[i]==1 :
				pSet.add(i)
	for ind in prevN :
                for i in range(ind, min(ind+maxDelay+1,len(p))):
                        if p[i]==-1 :
                                nSet.add(i)
	
	return (pSet, nSet)

def infByNode(n, g, maxDelay) :
	visited = set()
	infNet = nx.DiGraph() 
	pnDict = {}
	profit = 0
	stack = []
	stack.append(n)
	pnDict[n] = startInd(n)
	while len(stack)>0:
        	v = stack.pop()
		if v not in visited :
			visited.add(v)
			prevP, prevN = pnDict[v] 
			for s in g.successors(v) :
				p, n = nextInd(s, (prevP, prevN), g[v][s]['type'], maxDelay)
				if len(p)+len(n) > 0 :
					infNet.add_edge(v, s)
					profit += g.node[s]['label']
					stack.append(s)
					pnDict[s] = (p, n)
	return infNet, profit

def IM(network, maxDelay, k) :
	
	infNo = {}
	TFpath = {}

	TFSet = set([x for x in network.nodes() if network.node[x]['isTF']==1])
	repeat = 100
	seedSet = set()
	elseSet = TFSet
	newSeed = ''

	for j in range(k) :
		print 'TF '+str(j)
		for n in elseSet :
			infNo[n]=0
		for i in range(repeat) :
			# Produce G'
			g = network.copy()
			for (a,b) in network.edges() :
				if random.randint(1,1000)>g[a][b]['weight'] :
					g.remove_edge(a, b)

			# Influence set by seed
			seedInfSet = set()
			for TF in seedSet :
				(infNet, profit) = infByNode(TF, g, maxDelay)
				infSet = infNet.nodes()
				seedInfSet = seedInfSet | set(infSet)
			# If a TF is not influenced by current seeds then calculate profit
			for n in elseSet-seedInfSet :
				(infNet, profit) = infByNode(n, g, maxDelay)
				infSet = infNet.nodes()
				if len(infSet)>0 :
					infNo[n] += float(profit)/len(infSet)
		for n in elseSet :
			infNo[n]=float(infNo[n])/repeat
	
		newSeed = max(infNo.iterkeys(), key=(lambda x : infNo[x]))
		if infNo[newSeed]==0.0 :
			continue
		infNet, profit = infByNode(newSeed, network, maxDelay)
		TFpath[newSeed]=infNet
		print newSeed+'\t'+str(infNo[newSeed])
		seedSet.add(newSeed)
		elseSet.remove(newSeed)	
		del infNo[newSeed]
	return seedSet, TFpath

def keggDiagram(pathway, net_eg) :
	url = "http://www.kegg.jp/kegg-bin/show_pathway?map="+pathway+"&multi_query="
	for gene in net_eg.nodes() :
		url += gene.split('!')[1]+"+red%0D%0A"
	return url

def profilejs(profile) :
	temp = []
	for i in range(len(profile)) :
		color = '#ffffff'
		if profile[i]>=1.0 :
			color = '#ff0000'
		elif profile[i]>=0.8 and profile[i]<1.0:
			color = '#ff2a2a'
		elif profile[i]>=0.6 and profile[i]<0.8 :
			color = '#ff5555'
		elif profile[i]>=0.4 and profile[i]<0.6:
                        color = '#ff7f7f'
                elif profile[i]>=0.2 and profile[i]<0.4 :
                        color = '#ffaaaa'
		elif profile[i]>0 and profile[i]<0.2:
                        color = '#ffd4d4'
		elif profile[i]==0 :
			color = '#ffffff'
		elif profile[i]>-0.2 and profile[i]<0:
                        color = '#d4d4ff'
                elif profile[i]>-0.4 and profile[i]<=-0.2 :
                        color = '#aaaaff'
                elif profile[i]>0.6 and profile[i]<=-0.4:
                        color = '#7f7fff'
                elif profile[i]>0.8 and profile[i]<=-0.6 :
                        color = '#5555ff'
                elif profile[i]>-1.0 and profile[i]<=-0.8:
                        color = '#2a2aff'
		elif profile[i]<=-1.0 :
			color = '#0000ff'
		temp.append('p'+str(i+1)+':"'+color+'"')
	return ','.join(temp)

def cytoscapejs(webPath, pathway, pathNo, path, net_eg, knetwork, network, seedSet, numTP, TFpath) :
	file1 = open(webPath+'/'+pathway+'_'+str(pathNo+1)+'.html', 'w')
	file1.write('<head><link href="style.css" rel="stylesheet" />\n')
	file1.write('<script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js"></script>\n')
	file1.write('<script src="https://cdn.rawgit.com/cytoscape/cytoscape.js/d48294ac/dist/cytoscape.min.js"></script>\n')
	file1.write('<script src="https://unpkg.com/dagre@0.7.4/dist/dagre.js"></script>\n<script src="cytoscape-dagre.js"></script>\n')
	file1.write('<script src="'+pathway+'_'+str(pathNo+1)+'.js"></script></head>\n')
	file1.write('<body><div id="cy"></div></body></html>\n')
	file1.close()

	file2 = open(webPath+'/'+pathway+'_'+str(pathNo+1)+'.js', 'w')
	file2.write('window.addEventListener(\'DOMContentLoaded\', function(){\nvar cy = window.cy = cytoscape({\n  container: document.getElementById(\'cy\'),\n  boxSelectionEnabled: false, autounselectify: true, layout: {  name: \'dagre\'  },\n')
	file2.write("style: [ {  selector: \'node[type=\"single\"]\',\nstyle: {\"width\":\"30px\",\"height\":\"30px\",\"content\":\"data(id)\",\"border-color\":\"data(bc)\",\"border-width\": 3,\"pie-size\":\"100%\",\n")

	for i in range(numTP) :
		file2.write('"pie-'+str(i+1)+'-background-color": "data(p'+str(i+1)+')",\n')
		file2.write('"pie-'+str(i+1)+'-background-size": '+str(100/(numTP)))
		if i<numTP-1 :
			file2.write(',')
		else:
			file2.write('}\n},\n')
	file2.write('{ selector: \'edge\',\nstyle:{\n"width":"mapData(weight, 0, 1000, 1, 4)","target-arrow-shape":"data(shape)","opacity":1.0,"target-arrow-color":"black","curve-style":"bezier"}\n}],\n')

	nodes = set()
	edges = set()

	TFs = set()
	totTFs = set([x for x in network.nodes() if network.node[x]['isTF']==1])
	# TF to the first gene
	sources = set([x for x in net_eg.nodes() if net_eg.in_degree(x)==0])
	for gene in sources : 
		first = kgDict[gene.split('!')[1]]
		minLen = float('Inf')
		minSeed = ''
		for seed in seedSet :
			if first not in TFpath[seed].nodes() :
				continue
			try:
				length = nx.shortest_path_length(TFpath[seed], seed, first)
			except:
				length = float('Inf')
				continue
			if length < minLen :
				minLen = length
				minSeed = seed
		if minSeed =='' :
			continue
		TFs.add(minSeed)
		spath = nx.shortest_path(TFpath[minSeed], minSeed, first)
		for i in range(len(spath)-1) :
			if i==0 :	# TF
				nodes.add('{ data: { id: "'+spath[i]+'",'+profilejs(keggDict[gkDict[spath[i]]])+', bc: "blue", type: "single" }}')
			else :		# other
				if network.node[spath[i]]['isTF']==1 :
					color = 'blue'
				else :
					color = 'gray'
				nodes.add('{ data: { id: "'+spath[i]+'",'+profilejs(keggDict[gkDict[spath[i]]])+', bc: "'+color+'", type: "single" }}')
			if i<len(spath)-1 :
				shape = 'tee' if network[spath[i]][spath[i+1]]['type']==-1 else 'triangle'
				edges.add('{ data: { source: "'+spath[i]+'", target: "'+spath[i+1]+'", weight:'+str(network[spath[i]][spath[i+1]]['weight'])+', shape:"'+shape+'"}}')
	
	# perturbed path
	egDict = {}
	for (g1, g2) in net_eg.edges() :
		eg1 = g1.split('!')
		eg2 = g2.split('!')
		e1 = eg1[0]
		e2 = eg2[0]
		g1 = eg1[1]
		g2 = eg2[1]
		if e1 not in egDict :
			egDict[e1] = set()
		egDict[e1].add(g1)
		if e2 not in egDict :
                        egDict[e2] = set()
                egDict[e2].add(g2)
		
		g1 = kgDict[eg1[1]]
		g2 = kgDict[eg2[1]]

		shape = 'tee' if knetwork[e1][e2]['type']==-1 else 'triangle'
		edges.add('{ data: { source: "'+g1+'", target: "'+g2+'", weight: 1000, shape:"'+shape+'"}}')

	
	for e, gSet in egDict.iteritems() :
		# if it has to be a compound node
                if len(gSet)>1 :
                        nodes.add('{data:{id:"'+e+'", type: "compound"}}')
                        for g in gSet :
                                nodes.add('{data:{id:"'+kgDict[g]+'",'+profilejs(keggDict[g])+', parent:"'+e+'", bc: "green", type: "single" }}')
                # single node
                else :
                        g = gSet.pop()
                        nodes.add('{data:{id:"'+kgDict[g]+'",'+profilejs(keggDict[g])+', bc: "green", type: "single" }}')	

	# Print out
	file2.write('elements:{ nodes: ['+','.join(list(nodes))+'],\n edges: ['+','.join(list(edges))+']}\n});\n});')
	file2.close()
	
	return pathway+'_'+str(pathNo+1)+'.html', TFs

def visualize(resPath, KEGGPathFilePath, networks, network, pathDict, pathPval, DEGPval, seedSet, numTP, TFpath, threshold) :
	pnameFile = open(KEGGPathFilePath,'r')
	pnameDic = {}
	for path in pnameFile.readlines() :
                tp = path.split("\t")
                tp2 = tp[0].split(":")
                tp3 = tp[1].split(" - ")
                pnameDic[tp2[1]]=tp3[0]
        pnameFile.close()

	file1 = open(os.path.join(resPath, 'trap_pvalue.html'), 'w')
	file1.write('<!DOCTYPE html>\n<html>\n<head><link href="style.css" rel="stylesheet" /></head>\n<body>\n<table>\n')
	file1.write('<tr><th>Pathway</th><th class="name">Pathway name</th><th>DEG p-value</th><th>Path</th><th>Path p-value</th><th>Combined p-value</th><th>Link1</th><th>Link2</th><th>TFs</th></tr>\n')
	strs = []
	for pathway,tri in pathDict.iteritems() :
		pathlist = tri[0]	# the list of entry path
		genelist = tri[1]	# the network (ent!gene)
		if pathway not in pnameDic :
			continue
		p = pnameDic[pathway]
		if len(p)>30 :
			p = p[:30]+'...'
		s = ''
		row = str(1) if len(pathlist)<2 else str(len(pathlist))
		s+='<tr><td rowspan="'+row+'">'+pathway+'</td><td class="name" rowspan="'+row+'">'+p+'</td><td rowspan="'+row+'">'+ "{:5.3f}".format(DEGPval[pathway])+'</td>\n'
		
		# If there is perturbed path(s)
		for i in range(len(pathlist)) :
			st, comPval = scipy.stats.combine_pvalues([DEGPval[pathway], pathPval[pathway][i]])
			if i!=0 :
				s+='<tr>'
			pstr = "{:5.3f}".format(pathPval[pathway][i])
			cpstr = "{:5.3f}".format(comPval)
			if pathPval[pathway][i]< threshold :
				pstr = '<b>'+pstr+'</b>'
			if comPval< threshold :
				cpstr = '<b>'+cpstr+'</b>'
			s+='<td>Path'+str(i+1)+'</td><td>'+pstr+'</td><td>'+cpstr+'</td>\n'
			link1 = keggDiagram(pathway, genelist[i])
			link2, TFs = cytoscapejs(resPath, pathway, i, pathlist[i], genelist[i], networks[pathway], network, seedSet, numTP, TFpath) 
			s+='<td><a href="'+link1+'" target="diagram">KEGG</td>\n<td><a href="'+link2+'" target="cyto">CYTO</a></td><td>'+','.join(list(TFs))+'</td></tr>\n'

		# If there is no perturbed path
		if len(pathlist)==0 :
			for i in range(6) :
				if i==2 :
					s+='<td>'+"{:5.3f}".format(DEGPval[pathway])+'</td>'
				else :
					s+='<td></td>'
			s+='</tr>\n'
		strs.append(s)
	strs = sorted(strs, key = lambda x: float(re.findall(r"[0-1]\.[0-9]+", x)[0]))
	file1.write('\n'.join(strs))
	file1.write('</table></body></html>\n')
	file1.close()
	return		

def visualize2(webPath, networks, network, TFpath, labeled, labeled_path, pathDict, pathPval, DEGPval, threshold) :
        file1 = open(webPath+'/network_overview.html', 'w')
	file1.write('<head><link href="style.css" rel="stylesheet" />\n')
	file1.write('<script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js"></script>\n')
	file1.write('<script src="https://cdn.rawgit.com/cytoscape/cytoscape.js/d48294ac/dist/cytoscape.min.js"></script>\n')
	file1.write('<script src="https://unpkg.com/dagre@0.7.4/dist/dagre.js"></script>\n<script src="cytoscape-dagre.js"></script>\n')
	file1.write('<script src="network_overview.js"></script></head>\n')
	file1.write('<body><div id="cy"></div></body></html>\n')
	file1.close()

	file2 = open(webPath+'/network_overview.js', 'w')
	file2.write('window.addEventListener(\'DOMContentLoaded\', function(){\nvar cy = window.cy = cytoscape({\n  container: document.getElementById(\'cy\'),\n  boxSelectionEnabled: false, autounselectify: true, layout: {  name: \'dagre\'  },\n')
	file2.write("style: [ {  selector: \'node[type = \"single\"]',\nstyle: {\"width\":\"30px\",\"height\":\"30px\",\"content\":\"data(id)\",\"border-color\":\"data(bc)\",\"border-width\": 3,\"pie-size\":\"100%\",\n")

	for i in range(numTP) :
		file2.write('"pie-'+str(i+1)+'-background-color": "data(p'+str(i+1)+')",\n')
		file2.write('"pie-'+str(i+1)+'-background-size": '+str(100/(numTP)))
		if i<numTP-1 :
			file2.write(',')
		else:
			file2.write('}\n},\n')
	file2.write('{ selector: \'edge\',\nstyle:{\n"width":"mapData(weight, 0, 1000, 1, 4)","target-arrow-shape":"data(shape)","opacity":1.0,"target-arrow-color":"black","curve-style":"bezier"}\n}],\n')

	tempNet = nx.DiGraph()
        nodes = set()
        edges = set()

	# TF-Pathway
	for seed, infNet in TFpath.iteritems() :
		targets = set([x for x in infNet.nodes() if network.node[x]['label']>0])
		for t in targets :
			spath = nx.shortest_path(infNet, seed, t)
			for i in range(len(spath)) :
				color = 'gray'
				if network.node[spath[i]]['isTF']==1 :
					color = 'blue'
				if network.node[spath[i]]['label']>0 :
					color = 'green'
				nodes.add('{ data: { id: "'+spath[i]+'",'+profilejs(keggDict[gkDict[spath[i]]])+', bc: "'+color+'", type: "single" }}')
				if i<len(spath)-1 :
					shape = 'tee' if network[spath[i]][spath[i+1]]['type']==-1 else 'triangle'
					edges.add('{ data: { source: "'+spath[i]+'", target: "'+spath[i+1]+'", weight:'+str(network[spath[i]][spath[i+1]]['weight'])+', shape:"'+shape+'"}}')
	
	# Perturbed path
	parent = 1
	for pathway, tri in pathDict.iteritems() :
		if not DEGPval[pathway] < threshold :
			continue
		knetwork = networks[pathway]
		pathlist = tri[0]       # the list of entry path
                genelist = tri[1]       # the network (ent!gene)
		for i in range(len(pathlist)) :
			if not pathPval[pathway][i] < threshold :
				continue
			net_eg = genelist[i] 		 
			egDict = {}
			for (g1, g2) in net_eg.edges() :
				eg1 = g1.split('!')
				eg2 = g2.split('!')
				e1 = eg1[0]
				e2 = eg2[0]
				g1 = eg1[1]
				g2 = eg2[1]
				if e1 not in egDict :
					egDict[e1] = set()
				egDict[e1].add(g1)
				if e2 not in egDict :
					egDict[e2] = set()
				egDict[e2].add(g2)

				g1 = kgDict[eg1[1]]
				g2 = kgDict[eg2[1]]

				shape = 'tee' if knetwork[e1][e2]['type']==-1 else 'triangle'
				edges.add('{ data: { source: "'+g1+'", target: "'+g2+'", weight: 1000, shape:"'+shape+'"}}')

			for e, gSet in egDict.iteritems() :
				# if it has to be a compound node
				if len(gSet)>1 :
					nodes.add('{data:{id:"'+str(parent)+'", type: "compound"}}')
					for g in gSet :
						nodes.add('{data:{id:"'+kgDict[g]+'",'+profilejs(keggDict[g])+', parent:"'+str(parent)+'", bc: "green", type: "single" }}')
					parent +=1
				# single node
				else :
					g = gSet.pop()
					nodes.add('{data:{id:"'+kgDict[g]+'",'+profilejs(keggDict[g])+', bc: "green", type: "single" }}')

        # Print out
        file2.write('elements:{ nodes: ['+','.join(list(nodes))+'],\n edges: ['+','.join(list(edges))+']}\n});\n});')
        file2.close()

def makeFiles(webPath) :
	empty = open(webPath+'/empty.html','w')
	index = open(webPath+'/index.html','w')
	css = open(webPath+'/style.css','w')

	index.write('<!DOCTYPE html>\n<html><head><title>TRAP results</title><link href="style.css" rel="stylesheet" /></head>\n<frameset cols="850px,*">\n<frame src="trap_pvalue.html" name="left" id="left" frameborder="0" />\n<frameset rows="50%, 50%">\n<frame src="empty.html" name="diagram" id="diagram" frameborder="0" />\n <frame src="empty.html" name="cyto" id="cyto" frameborder="0" />\n</frameset>\n</frameset>\n<body></body>\n</html>\n')
	css.write('body { font: 14px helvetica neue, helvetica, arial, sans-serif; }\nth, td { border: 1px black; width: 100px; text-align:center; }\n .name { width: 300px;}\n td.name { text-align: left; }\n#cy { height: 100%; width: 100%; position: absolute; left: 0; top: 0; }\n')
	empty.close()
	index.close()
	css.close()

if __name__ == "__main__" :
	conf = open(sys.argv[1],'r')
	for l in conf.readlines() :
		tp = l.rstrip().split()
		if tp[0]=='outDir' :
			resPath = tp[1]
		elif tp[0]=='numTP' :
			numTP = int(tp[1])
		elif tp[0]=='threshold' :
			threshold = float(tp[1])
		elif tp[0]=='maxDelay' :
			maxDelay = int(tp[1])
		elif tp[0]=='single' :
			single = int(tp[1])
		elif tp[0]=='k' :
			k = int(tp[1])
		elif tp[0]=='species' :
			species = tp[1]
		elif tp[0]=='downloadFiles' :
			download = tp[1]
		elif tp[0]=='KEGGGeneFilePath' :
			KEGGgeneFilePath = tp[1]
		elif tp[0]=='KEGGPathFilePath' :
			KEGGPathFilePath = tp[1]
		elif tp[0]=='GRNFilePath' :
			GRNfile = tp[1]
		elif tp[0]=='PINFilePath' :
			PPIfile = tp[1]
		elif tp[0]=='xmlFilePath' :
			xmlPath = tp[1]
	if single==1 :
		numTP-=1

	webPath = os.path.join(resPath, 'webpage')
	if not os.path.exists(resPath) :
		os.mkdir(resPath) 
	if not os.path.exists(webPath) :
		os.mkdir(webPath)
	shutil.copy("legend.PNG", webPath)	
	shutil.copy("cytoscape-dagre.js", webPath)
	conf.close()
	
	pfile = 'DEG_profile.txt'
	profile = 'DEG'
	if download.upper()=="YES" :
		KEGGgeneFilePath = os.path.join(species, 'KEGG_GeneSymbols.txt')
		KEGGPathFilePath = os.path.join(species, 'KEGG_PathwayList.txt')
		xmlPath = os.path.join(species, 'xmlFiles')		
 
	print ". genDict"
	species = genDict(os.path.join(resPath, pfile), KEGGgeneFilePath, numTP)

	print ". parseXml"
	(networks, entries, DEGs) = parseXml(xmlPath, keggDict)

        print ". DEGAnalysis"
        DEGPval = DEGAnalysis(DEGs)
	print ". pathwayAnalysis"
	pathPval, pathDict, labeled, labeled_path = pathwayAnalysis(numTP, networks, entries, DEGs, DEGPval, maxDelay, profile, threshold)

	network = nx.DiGraph()
	print ". read PIN and GRN"
	network = readFiles(GRNfile, PPIfile, species, network, labeled, maxDelay)

	# Reweight the matched pair edges
	for (g1, g2) in network.edges() :
		if not checkProfile(keggDict[gkDict[g1]], keggDict[gkDict[g2]], network[g1][g2]['type'], maxDelay) :
			network[g1][g2]['weight']/=2
			network[g1][g2]['match']='k'
			network.remove_edge(g1, g2)
		else :
			network[g1][g2]['match']='m'
	
	print '.. the number of TFs and labeled genes in time bounded network'
	print '# of TFs :'+str(len([x for x in network.nodes() if network.node[x]['isTF']==1]))
	print '# of labeled genes :'+str(len([x for x in network.nodes() if network.node[x]['label']>0]))

	print ". Influence Maximization"
	seedSet, TFpath = IM(network, maxDelay, k)
	for s in seedSet :
		print gkDict[s]+'\t'+s

	print ". Visualization"
	visualize(webPath, KEGGPathFilePath, networks, network, pathDict, pathPval, DEGPval, seedSet, numTP, TFpath, threshold)
	visualize2(webPath, networks, network, TFpath, labeled, labeled_path, pathDict, pathPval, DEGPval, threshold)
	makeFiles(webPath)
