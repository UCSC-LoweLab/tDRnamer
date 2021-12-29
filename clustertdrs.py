#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import subprocess
from tdrdbutils import *
from distutils.spawn import find_executable
import time


from collections import defaultdict


'''
A_DSLE_CSLE_TSLE_Z_CCA

5027c,5027b,5026a,5026b,5026c
>5027c
-->5027b
---->5026a

>5026a
GTTTCCGTAGTGTAGTGG
>5026b
GTTTCCGTAGTGTAGTGGTCATC
>5026c
GTTTCCGTAGTGTAGTGGTCATCACGTTCGC
>5027b
GTTTCCGTAGTGTAGTGGTTATC
>5027c
GTTTCCGTAGTGTAGTGGTTATCACGTTCGC

defaultdict(<type 'set'>, {'5026b': set(['5026c']), '5027b': set(['5027c']), '5017a': set(['5017b', '5017c']), '5017b': set(['5017c']), '5026a': set(['5027c', '5027b', '5026b', '5026c'])})
5027c,5027b,5026a,5026b,5026c



'''

def getduplicates(group,duplicates):
    currgroup = set(group)
    dupedtrnas = set()
    for currgene in currgroup:
        if currgene not in dupedtrnas and currgene in duplicates:
            dupset = set(duplicates[currgene]) | set([currgene])
            yield  dupset
            dupedtrnas = dupedtrnas | dupset

def printchildren(group,parents, duplicates):
    recursechildren(group,parents, duplicates,0)
class tree:
    def __init__(self, name, children):
        self.name = name
        self.children = children
def recursechildren(group,parents,duplicates, depth):
    currgroup = set(group)
    dupedtrnas = set()
    for currgene in currgroup:
        if currgene in  dupedtrnas:
            continue
        currparents = set(curr for curr in parents[currgene] if curr in currgroup)
        #print >>sys.stderr, currparents

    
        if len(currparents) == 0:
            if currgene in duplicates:
                dupfrags = set([currgene]) | duplicates[currgene]
            else:
                dupfrags = set([currgene])
            print ("--"*depth)+'>'+",".join(dupfrags)
            currchildren = set(curr for curr in currgroup if currgene in parents[curr])
            #print >>sys.stderr, currparents
            #currtree = tree(dupfrags, list(recursechildren(currchildren,parents,duplicates, depth  = depth + 1)))
            #yield dupfrags

            recursechildren(currchildren,parents,duplicates, depth  = depth + 1)
            dupedtrnas = dupedtrnas | dupfrags

class sprinzelrange:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        
    def coverage(self, other):
        start = max([self.start,other.start])
        end = min([self.end,other.end])
        
        if end - start < 0:
            return 0
        else:
            #print >>sys.stderr, end - start
            return end - start
    def overlaps(self, other):
        return self.start <= other.start and self.end >= other.end
    def equals(self, other):
        return self.start == other.start and self.end == other.end
            
class allgroup:
   def __init__(self):
       self.allset = set()
       self.groupset = list()
   def addpair(self, fir, sec):
       if fir in self.allset and sec in self.allset:
           #print >>sys.stderr, "***"+fir+":"+sec
           firset = None
           secset = None
           for i, currset in enumerate(self.groupset):
               if fir in currset and sec in currset:
                   secset = i
                   firset = i
               elif fir in currset:
                   firset = i
               elif sec in currset:
                   secset = i
           if firset != secset:
               self.groupset[firset] = self.groupset[firset] | self.groupset[secset]
               del self.groupset[secset]
       elif fir in self.allset:
           for currset in self.groupset:
               if fir in currset:
                   currset.add(sec)
       elif sec in self.allset:
           for currset in self.groupset:
               if sec in currset:
                   currset.add(fir)
       else:
           self.groupset.append(set([fir, sec]))
       self.allset.add(fir)
       self.allset.add(sec)
   def addunlinkpair(self, fir, sec):
       if fir in self.allset:
           for currset in self.groupset:
               if fir in currset:
                   currset.add(sec)
       if sec in self.allset:
           for currset in self.groupset:
               if sec in currset:
                   currset.add(fir)
       else:
           self.groupset.append(set([fir, sec]))
       self.allset.add(fir)
       self.allset.add(sec)
   def addsingle(self, fir):
       if fir in self.allset:
           pass
       else:
           self.groupset.append(set([fir]))
           self.allset.add(fir)
   def allsets(self):
       return self.groupset
   def allgenes(self):
       return self.allset

def createdomains(domains):
    domainlist  = list(domains.split("_"))
    #print domainlist
    domainset = list()
    for currdomain in domainlist:
        if currdomain[0] == "D":
            for currsub in currdomain[1:]:
                domainset.append(currdomain[0]+currsub)
        elif currdomain[0] == "C": 
            for currsub in currdomain[1:]:
                domainset.append(currdomain[0]+currsub)
        elif currdomain[0] == "T":
            for currsub in currdomain[1:]:
                domainset.append(currdomain[0]+currsub)
        else:
            domainset.append(currdomain)
    return domain(domainset)
        
class domain:
    def __init__(self, domains = set()):

        self.domains = set(domains)
    def listdomains(self):
        return self.domains
    def comparedomains(self, other):
        pass
    def adddomain(self, other):
        self.domains.add(other)
    def __len__(self):
        return len(self.domains)
    def domainstring(self): #A_DSLE_CSL
        outstring = list()
        if "A" in self.domains:
            outstring.append("A")
        subdstr = ""
        for currsub in "SLE":
            if "D"+currsub in self.domains:
                subdstr = subdstr + currsub
        if len(subdstr) > 0:
            outstring.append("D"+subdstr)
        subdstr = ""
        for currsub in "SLE":
            if "C"+currsub in self.domains:
                subdstr = subdstr + currsub
        if len(subdstr) > 0:
            outstring.append("C"+subdstr)
        if "V" in self.domains:
            outstring.append("V")
        subdstr = ""
        for currsub in "SLE":
            if "T"+currsub in self.domains:
                subdstr = subdstr + currsub
        if len(subdstr) > 0:
            outstring.append("T"+subdstr)
        if "Z" in self.domains:
            outstring.append("Z")
        if "CCA" in self.domains:
            outstring.append("CCA")
        return "_".join(outstring)

        
    
overlapgroups = defaultdict(lambda: defaultdict(set))
domainconfigs = defaultdict(int)
sprinzelseqs = dict()
tdrfile = open(sys.argv[1])

stkfile = sys.argv[2]
fragstk = list(readrnastk(open(stkfile, "r")))[0]
overlaponly = False
classicfrags = False
allfrags = list()
fraglengths = dict()
aminos = dict()
fragtypes = dict()
badfrags = set()
mismatches = dict()
for currline in tdrfile:
    if currline.startswith("frag_name"):
        continue
    fields = currline.rstrip().split("\t")
    if len(fields) < 9 or not (fields[1] == "tRNA" or fields[1] == "pretRNA"):
        continue
    seqname = fields[0]
    if seqname not in set(fragstk.getseqnames()):
        continue
    
    fraglengths[seqname] = int(fields[4])
    allfrags.append(seqname)
    #print >>sys.stderr, seqname
    sprinzelstring = tuple(fields[5].split(".."))
    #print >>sys.stderr, fields[9]
    currsprinzel = sprinzelrange(sprinzelstring[0],sprinzelstring[1])
    sprinzelseqs[seqname] = currsprinzel
    aminos[seqname] = set(fields[7].split(","))
    trnatranscripts = set(fields[3].split(","))

    if fields[10] == "None" or fields[10] == "NA":
        mismatches[seqname] = set()
    else:
        mismatches[seqname] = set(fields[10].split(","))
    fragtypes[seqname] = fields[2].split("-")[1]
    if len(aminos[seqname]) > 1 and int(fields[4]) < 25 :
        badfrags.add(seqname)
        
    #print >>sys.stderr, fields
    #print >>sys.stderr, len(fields)
    #print >>sys.stderr, domainconfigs
    #print >>sys.stderr, len(domainconfigs)
    
    if fields[1] == "tRNA":
        domainconfigs[fields[11]] += 1
        trnadomains = createdomains(fields[11])
        for currtrna in trnatranscripts:
            for currdomain in trnadomains.listdomains():
                overlapgroups[currtrna][currdomain].add(seqname)
    elif fields[1] == "pretRNA":

        #print >>sys.stderr, trnatranscripts
        #overlapgroups for pretrnas
        for currtrna in trnatranscripts:
            overlapgroups[currtrna]["pretrna"].add(seqname)
#for curr in domainconfigs.iterkeys():
    #print curr+":"+str(domainconfigs[curr])
#print >>sys.stderr, fragtypes["frag1164:64"]
#print >>sys.stderr, "**||"
overlappairs = dict()
coverpairs = dict()
trnagroups = allgroup()
equalgroups = allgroup()
parentfrags = defaultdict(set)
equalfrags = defaultdict(set)
connections = defaultdict(int)
smallfrags = defaultdict(set)
smallcutoff = 25

for currtrna in overlapgroups.iterkeys():
    fragset = set()
    for currdomain in overlapgroups[currtrna].iterkeys():
        for currfrag in overlapgroups[currtrna][currdomain]:
            fragset.add(currfrag)
    #print >>sys.stderr, fraglist
    fraglist = list(fragset)
    for trnapair in itertools.combinations(fraglist, 2):
        firtrna = list(trnapair)[0]
        sectrna = list(trnapair)[1]
        #print >>sys.stderr, trnapair
        if firtrna == sectrna:
            #print >>sys.stderr, "*equals"
            continue
        if sprinzelseqs[firtrna].equals(sprinzelseqs[sectrna]):
            equalgroups.addpair(firtrna, sectrna)
            equalfrags[sectrna].add(firtrna)
            equalfrags[firtrna].add(sectrna)
            pass
            #coverpairs[firtrna] = sectrna 
        elif sprinzelseqs[firtrna].overlaps(sprinzelseqs[sectrna]):
            parentfrags[sectrna].add(firtrna)
           
        elif sprinzelseqs[sectrna].overlaps(sprinzelseqs[firtrna]):
            parentfrags[firtrna].add(sectrna)
            #trnagroups.addpair(firtrna, sectrna)

        if True: #sprinzelseqs[firtrna].coverage(sprinzelseqs[sectrna]) > 5: no longer works with new sprinzel
            
            if firtrna in badfrags:
                smallfrags[sectrna].add(firtrna)
                
            elif sectrna in badfrags:
                smallfrags[firtrna].add(sectrna)
                #trnagroups.addunlinkpair(firtrna, sectrna)
            elif classicfrags and fragtypes[firtrna] == fragtypes[sectrna]:
                trnagroups.addpair(firtrna, sectrna)
                pass
            elif not classicfrags:
                trnagroups.addpair(firtrna, sectrna)
            
            connections[firtrna] += 1
            connections[sectrna] += 1
#sys.exit() 
#print >>sys.stderr, "**|||"
for curr in sorted(connections.iterkeys(), key = lambda x: -connections[x]):
    pass
    #print curr+":"+str(connections[curr])
#print "Overlapping fragment groups"
#print >>sys.stderr, fragstk.consensus
#sys.exit()
#remove overly connecting tRNAs

orphanfrags = badfrags - reduce(set.union, (smallfrags.values()), set()) 
for currtrna in allfrags:
    if currtrna not in trnagroups.allgenes() and currtrna not in badfrags:
        trnagroups.addsingle(currtrna)
    elif currtrna in orphanfrags:
        trnagroups.addsingle(currtrna)
#print >>sys.stderr, "**||||"        
smallprinted = set()
for i, currgroup in enumerate(trnagroups.allsets()):
    if all(curr in badfrags for curr in currgroup):
        #here
        #continue
        pass
    smalltrnas = reduce(set.union, (smallfrags[curr] for curr in currgroup), set())
    smalltrnas = smalltrnas - smallprinted
    smallprinted = smallprinted | smalltrnas 
    print "Fragment Cluster "+str(i + 1)+": "+",".join(currgroup | smalltrnas)

    #print >>sys.stderr,fragtree
    #sys.exit()
    clustertrnas = set()
    subtrnas = set()
    clusteraminos = set()
    print "Fragment Coverage"
    
    smallprinted = smallprinted | smalltrnas 
    hasvarloop = False
    for currtrna in overlapgroups.iterkeys():
        trnadomains = domain()
        for currdomain in overlapgroups[currtrna].iterkeys():
            if len(currgroup & overlapgroups[currtrna][currdomain]) > 0:
                trnadomains.adddomain(currdomain)
                if currdomain == "V":
                    hasvarloop = True
        
        if len(trnadomains) > 0:
            #print >>sys.stderr, currtrna
            print currtrna+":"+trnadomains.domainstring()
            clustertrnas.add(currtrna)
    clustertrnas = list(trnasort(clustertrnas))
    print "Fragment Aminos"
    clusteraminos = reduce(set.union,(aminos[curr] for curr in currgroup), set())

    print ",".join(clusteraminos)
    if len(clusteraminos) > 1:
        pass
        #print >>sys.stderr, "****"+"Cluster "+str(i + 1)+":"+str( len(currgroup)) +str(clusteraminos)
    #print >>sys.stderr, smalltrnas
    #reduce(domain.adddomains,currgroup)
    #Removing fragment clustering for now because alignment shows it better
    #print "Fragment Tree"
    #fragtree = printchildren(currgroup, parentfrags, equalfrags)
    duplist = list(getduplicates(currgroup, equalfrags))
    if len(duplist) > 0:
        #print "Duplicate fragments: "
        for currdups in duplist:
            pass
            #print "\t"+",".join(currdups)
    #print >>sys.stderr, "||||**"    
    mismatchcounts = defaultdict(int) 
    for curr in itertools.chain.from_iterable(mismatches[curr] for curr in currgroup):
        mismatchcounts[curr] += 1
    if len(mismatchcounts.keys()) > 0:
        print "Fragment Mismatches"
        for curr in mismatchcounts.iterkeys():
            print curr +":"+str(mismatchcounts[curr])
    print "Fragment Alignment"
    #if None in currgroup | clustertrnas:
    trnaannotatation = {curr:str(i+1) for i, curr in enumerate(clustertrnas)}

    if all(len(currtrna.split(":")) > 1 and is_number(currtrna.split(":")[1]) for currtrna in (currgroup | smalltrnas)):
        trnaorder = list(sorted(currgroup | smalltrnas, key = lambda x: -int(x.split(":")[1])))
    else:
        trnaorder = list(sorted(currgroup | smalltrnas, key = lambda x: [sprinzelseqs[x].start, len(fragstk.aligns[x].replace(".",""))]))
    for currfrag in trnaorder:
        fragtrnas = set(curr for curr in overlapgroups.keys() if currfrag in set(itertools.chain.from_iterable(overlapgroups[curr].itervalues())))
        fragtrnas = list(curr for curr in clustertrnas if curr in fragtrnas )
        trnaannotatation[currfrag] = ",".join(trnaannotatation[curr] for curr in fragtrnas)
        #print >>sys.stderr, fragtrnas
    #print >>sys.stderr, clustertrnas
    curralign = fragstk.getsubset(list(clustertrnas) + list(trnaorder)) #here is the trouble line

    if not hasvarloop:
        trnanums = gettrnanums(curralign, 0)
        varlooppos = list(i for i, currnum in enumerate(trnanums) if 'e' in currnum)
        curralign = curralign.trimcolumns(varlooppos)
    curralign.printstk(ordering = list(clustertrnas) + list(trnaorder), annotation = trnaannotatation)
    pass
#print >>sys.stderr, "**||||"      
groupnum = len(list(trnagroups.allsets())) + 1
#smallprinted are the small RNAs
for currtrna in allfrags:
    if  currtrna not in smallprinted and currtrna not in trnagroups.allgenes():
        print "Fragment Singleton "+str(groupnum)+": "+currtrna
        groupnum += 1
    
    

for currtrna in coverpairs.iterkeys():
    #print currtrna+":"+coverpairs[currtrna]
    pass
            
    