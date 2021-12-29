#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import pysam
import subprocess
from tdrdbutils import *


from collections import defaultdict

'''
Here is where I need to use the tRNA ontology between mature tRNAs and chromosomes

bowtie2 -f -x /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-tRNAs -k 100 --very-sensitive --ignore-quals --np 5 --reorder -p 4 -U testfrags.fa 2>bowtieerr.txt | nametdrs.py /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-trnatable.txt /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-maturetRNAs.bed
'''
minnontrnasize = 20

maxmaps = 50
header = list(["seq_name","source_type",'tdR_name','source_trnas','length','Sprinzl_range','anticodons','isotypes','mismatches','indels','mismatch_locations','trna_domains'])
#tRFdb_ID	source_type	tdR_name	source_trnas	length	Sprinzl_range	anticodons	isotypes	mismatches	indels	mismatch_locations	trna_domains

gapchars = set("-._~")
def is_number(num):
    try:
        int(num)
        return True
    except ValueError:
        return False

def ishead(trnanum):
    return trnanum.startswith("head")
def istail(trnanum):
    return trnanum.startswith("tail")
def ismargin(trnanum):
    return ishead(trnanum) or istail(trnanum)
def ifelse(a, b, c):
    if a:
        return b
    else:
        return c
def readcoverage(chromrange, read, readrefname):
    
    readstart = read.pos
    readend = read.aend
    if chromrange.strand == ifelse(read.is_reverse, "-","+") and chromrange.chrom == readrefname:
        start = max([chromrange.start,readstart])
        end = min([chromrange.end,readend])
        
        if end - start < 0:
            return 0
        else:
            return end - start
    else:
        return 0


def varloopinfo(trnaname, start, end, stkalign, trnanums):
    seqpos = 0
    i = 0
    result = None
    vartotal = 0
    for alignpos, curr in enumerate(stkalign.aligns[trnaname]):
        if curr in gapchars:
            pass
        else:
            i += 1
            if start < i < end:
                if 'e' in trnanums[alignpos]:
                    vartotal += 1
    return vartotal                                                                                
                    

#print >>sys.stderr, len(positions)
#this gets the tRNA numbers by the sprinzel numbering system
#dompositions = {'accend' : 7,'dstart':10,"dend":25,"acstart":27,'acend':43,'tstart':49,'tend':65,'accstart':66}
dompositions = {'accend' : 7,'dstart':10,"dloopstart":13,"dloopend":21,"dend":25,"acstart":27,'acloopstart':31,'acloopend':39,'acend':43,'tstart':49,'tloopstart':54,'tloopend':60,'tend':65,'accstart':66,'tail':74}
domnames = ['accend','dstart',"dloopstart","dloopend","dend","acstart",'acloopstart','acloopend','acend','tstart','tloopstart','tloopend','tend','accstart','tail']
armnames = ['accend','dstart',"dend","acstart",'acend','tstart','tend','accstart','tail']
domabbv = {'accend' : 'AS','dstart':'DS',"dloopstart":'DLS',"dloopend":'DLE',"dend":'DE',"acstart":'CS','acloopstart':'CLS','acloopend':'CLE','acend':'CE','varloop':'VA','tstart':'TS','tloopstart':'TLS','tloopend':'TLE','tend':'TE','accstart':'AE','tail':'CCA'}

newdomabbv = {'accend' : 'AS','dstart':'DS',"dloop":'DL',"dend":'DE',"acstart":'CS','acloop':'CL','acend':'CE','varloop':'V','tstart':'TS','tloop':'TL','tend':'TE','accstart':'Z','tail':'CCA'}
newdompos = {'accend' : [0,7],'dstart':[10,12],"dloop":[13,21],"dend":[22,25],"acstart":[27,30],'acloop':[31,39],'acend':[40,43],'tstart':[49,53],'tloop':[54,60],'tend':[61,65],'accstart':[66,73],'tail':[74,76]} 
newdomnames = ['accend','dstart',"dloop","dend","acstart",'acloop','acend','tstart','tloop','tend','accstart','tail']



newdompos = {'accend' : [1,7],'dstart':[10,13],"dloop":[14,21],"dend":[22,25],"acstart":[27,31],'acloop':[32,38],'acend':[39,43],'varm':[45,48],'tstart':[49,53],'tloop':[54,60],'tend':[61,65],'accstart':[66,72],'tail':[74,76]} 
arms = {'darm':['dstart',"dloop","dend"],'acarm':["acstart",'acloop','acend'],'tarm':['tstart','tloop','tend']} 
newarmabbv = {'accend' : 'A','darm':'D',"acarm":'C','tarm':'T','varm':"V",'accstart':'Z','tail':'CCA'}
armlist = ['accend','darm',"acarm",'varm','tarm','accstart','tail']
armsuffix = {'dstart' : 'S','acstart':'S',"tstart":'S','dloop' : 'L','acloop':'L',"tloop":'L','dend' : 'E','acend':'E',"tend":'E'}

'''
def getdomains(trnastart, trnaend, varlooplength = 0):
    ['accend','darm','acarm','tarm','acend']
    fragelements = list()
    if is_number(trnastart):
        start = int(trnastart)
        
    elif trnastart.startswith("e"):
        start = 45
    elif trnastart in sprinzeladditional:
        start = int(trnastart[:-1])
    else:
        print >>sys.stderr, "Invalid start position "+str(trnastart)
        sys.exit(1)
    if is_number(trnaend):
        end = int(trnaend)
        
    elif trnaend.startswith("e"):
        end = 46
    elif trnaend in sprinzeladditional:
        end = int(trnaend[:-1])
    else:
        print >>sys.stderr, "Invalid end position "+str(trnaend)
        sys.exit(1)
        
    for currdom in newdomnames:
        if start < newdompos[currdom][1] or end > newdompos[currdom][0]:
            fragelements.append(currdom)
        if currdom == 'acend': # and varlooplength > 1:
            fragelements.append('varloop')
    return fragelements
'''
def domainalign(trnanums ):
    currdom = None
    domlist = list()
    for currpos in trnanums:
        if is_number(currpos):
            for dom in newdompos.iterkeys():
                if int(currpos) == newdompos[dom][0]:
                    currdom = dom
                elif int(currpos) == newdompos[dom][1]:
                    currdom = None
                
        domlist.append(currdom)
        
    return domlist
        
def getdomainlist(domset):       
    domstring = list()
    for currarm in armlist:
        
        
        if currarm in arms:
            suff = ''
            for currseg in arms[currarm]:
                if currseg in domset:
                    suff += armsuffix[currseg]
            if len(suff) > 0:        
                domstring.append(newarmabbv[currarm] + suff)
        elif currarm in domset:
            domstring.append(newarmabbv[currarm])
    #print >>sys.stderr, domstring
    return domstring
       
'''
def getarms(trnastart, trnaend, varlooplength = 0):
    if is_number(trnastart):
        start = int(trnastart)
        
    elif trnastart.startswith("e"):
        start = 45
    elif trnastart in sprinzeladditional:
        start = int(trnastart[:-1])
    else:
        print >>sys.stderr, "Invalid start position " +str(trnastart)
        sys.exit(1)
    if is_number(trnaend):
        end = int(trnaend)
    elif trnaend.startswith("e"):
        end = 46
    elif trnaend in sprinzeladditional:
        end = int(trnaend[:-1])
    else:
        print >>sys.stderr, "Invalid end position "+str(trnaend)
        sys.exit(1)                                                                                                           
        
    fragelements = list()
    for currdom in armlist:
        if currdom in arms:
            subdoms = ''
            #print >>sys.stderr, start
            #print >>sys.stderr, newdompos[arms[currdom][0]][1]
            if start <= newdompos[arms[currdom][0]][1] and end >= newdompos[arms[currdom][0]][0]:
                subdoms += 'S'
            if start <= newdompos[arms[currdom][1]][1] and end >= newdompos[arms[currdom][1]][0]:
                subdoms += 'L'
            if start <= newdompos[arms[currdom][2]][1] and end >= newdompos[arms[currdom][2]][0]:
                subdoms += 'E'
            if len(subdoms) > 0:
                fragelements.append(newarmabbv[currdom]+subdoms)
        elif currdom == 'vloop':
            if trnastart.startswith("e") or trnaend.startswith('e'):
                fragelements.append(newarmabbv[currdom])
            elif start <= 48  and end >= 45: #and varlooplength > 1:  #
                fragelements.append(newarmabbv[currdom])
        else:
            if start < newdompos[currdom][1] and end > newdompos[currdom][0]:
                fragelements.append(newarmabbv[currdom])
        #if currdom == 'acarm' and varlooplength > 1 and end >= 46:
            #fragelements.append('V')
    return fragelements
'''

def findtag(read, tag):
    for currtag in read.tags:
        if currtag[0] == tag:
            return currtag[1]
def readmdtag(currread):
    digitre = re.compile(r"^(\d+)(.*)")
    mdtag = findtag(currread, "MD")
    currtag = mdtag
    currpos = currread.pos
    mismatches = dict()
    #print >>sys.stderr, currtag
    while len(currtag) > 0:
        digitmatch = digitre.search(currtag)
        if digitmatch:
            currpos += int(digitmatch.group(1))
            currtag = digitmatch.group(2)
        else:
            if currtag[0] == '^':
                mismatches[currpos + 1] = currtag[0]+currtag[1]
                currtag = currtag[2:]
                
            else:
                mismatches[currpos + 1] = currtag[0]
                currtag = currtag[1:]
                currpos += 1
    return mismatches
    
def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)

        
        #print >>sys.stderr, ",".join(bamfile.getrname(curr.tid) for curr in newset)
def processbatch(newset, bamfile, trnaalign,trnatranscripts, stkmargin = 0): #frag1471:29

    for currread in newset:
        chromname = bamfile.getrname(currread.tid)
        if chromname in trnatranscripts:
            #locusnames.append(currpretrna)
            trnaalign = addreadalign(trnaalign, chromname, currread, trnatranscripts[chromname].start, stkmargin)
        else:
            for currtrna in trnatranscripts.iterkeys():
                if readcoverage(trnatranscripts[currtrna], currread,chromname) > 10:
                    #locusnames.append(currpretrna)
                    if trnatranscripts[currtrna].strand == "-":
                        
                        trnaalign = addreadalign(trnaalign, currtrna, currread, trnatranscripts[currtrna].end, stkmargin,  reverse =  trnatranscripts[currtrna].strand == "-")
                    else:
                        trnaalign = addreadalign(trnaalign, currtrna, currread, trnatranscripts[currtrna].start, stkmargin,  reverse =  trnatranscripts[currtrna].strand == "-")
                    #trnastart = gettrnaposition(locusname, readend, locialign, locuspositionnums, returnpos = 'prev',start = pretrnalocs[locusname].end, margin = locimargin,includemargins = True, reverse = pretrnalocs[locusname].strand == '-')            
                    #trnaend = gettrnaposition(locusname,readstart+1  , locialign, locuspositionnums, returnpos = 'next',  start = pretrnalocs[locusname].end, margin = locimargin, includemargins = True,reverse = pretrnalocs[locusname].strand == '-')
                
                    break
    #print >>sys.stderr, "||"
    return trnaalign
    
    
    
'''
trnalign at this point should include flanking sequence as well
'''
def addreadalign(trnaalign, chrom, read, trnastart, margin, reverse = False):
    readstart = read.pos
    readend = read.aend
    
    #if reverse:
    #    resultpos = getalignposition(trnaname,position - start+ margin, stkalign )
    #    print >>sys.stderr, "||"
    #else:
    #    resultpos = len(trnanums) - getalignposition(trnaname,position - start + margin, stkalign )
    if reverse: 
        alignend = getalignposition(chrom,trnastart - readstart + margin, trnaalign )
        alignstart = getalignposition(chrom,trnastart - readend  + margin, trnaalign)
        #print >>sys.stderr, "||"
    else:
        #print >>sys.stderr, chrom
        #print >>sys.stderr, readstart
        #print >>sys.stderr, trnastart 
        #print >>sys.stderr, readstart - trnastart+ margin
        #print >>sys.stderr, "||"
        alignstart = getalignposition(chrom,readstart - trnastart+ margin, trnaalign )
        alignend = getalignposition(chrom,readend - trnastart + margin, trnaalign)
    #print >>sys.stderr, alignstart
    #print "-"*alignstart + read.seq + "-"*(trnaalign.alignlength - alignend)
    alnstring = "."*alignstart
    #print >>sys.stderr, alnstring
    currpos = alignstart
    #flatten out CIGAR string into a list
    cigar = list(itertools.chain.from_iterable(itertools.repeat(operation,length ) for operation, length in read.cigartuples))
    cigarpoint = 0
    readpos = 0
    newaligns = dict()
    newstruct = ''
    newconsensus = ''
    for curr in trnaalign.aligns.iterkeys():
        newaligns[curr] = trnaalign.aligns[curr][:alignstart]
    newstruct = trnaalign.currstruct[:alignstart]
    newconsensus = trnaalign.consensus[:alignstart]
    #print >>sys.stderr, chrom
    #print >>sys.stderr, currpos
    #print >>sys.stderr, readstart
    #print >>sys.stderr, len(trnaalign.aligns[curr])
    #print >>sys.stderr, len(list(cigar))
    #print >>sys.stderr, cigar
    #print >>sys.stderr, trnaalign.aligns[chrom]
    #print >>sys.stderr, len(trnaalign.aligns[chrom])
    readseq = read.seq.replace("T","U")
    if reverse:
        #print >>sys.stderr, "**||"
        readseq = revcom(readseq,rnamode = True)
    #need to add things here for dealing with shit before and after the alignment
    for currtype in cigar:
        RnaAlignment(newaligns, newstruct, newconsensus)
        if currtype == 0:  
            #print >>sys.stderr, currpos
            while trnaalign.aligns[chrom][currpos] in gapchars:
                alnstring += trnaalign.aligns[chrom][currpos]
                for curr in newaligns.iterkeys():
                    newaligns[curr] += trnaalign.aligns[curr][currpos]
                newstruct += trnaalign.currstruct[currpos]
                newconsensus += trnaalign.consensus[currpos]
                currpos += 1
            for curr in newaligns.iterkeys():
                newaligns[curr] += trnaalign.aligns[curr][currpos]
            newstruct += trnaalign.currstruct[currpos]
            newconsensus += trnaalign.consensus[currpos]
            alnstring += readseq[readpos]
            currpos += 1
            readpos += 1
            
        elif currtype == 1: #insertion to the reference
            alnstring += readseq[readpos]

            if trnaalign.aligns[chrom][currpos] not in gapchars:
                for curr in trnaalign.aligns.iterkeys():
                    newaligns[curr] += '-'
                newstruct += '.'
                newconsensus += '.'
            else:
                for curr in trnaalign.aligns.iterkeys():
                    newaligns[curr] += trnaalign.aligns[curr][currpos]
                newstruct += trnaalign.currstruct[currpos]
                newconsensus += trnaalign.consensus[currpos]
                currpos += 1
                
            readpos += 1
            #alnstring += '-'
        elif currtype == 2: #deletion from the reference
            while trnaalign.aligns[chrom][currpos] in gapchars:
                alnstring += trnaalign.aligns[chrom][currpos]
                for curr in newaligns.iterkeys():
                    newaligns[curr] += trnaalign.aligns[curr][currpos]
                newstruct += trnaalign.currstruct[currpos]
                newconsensus += trnaalign.consensus[currpos]
                currpos += 1
            for curr in trnaalign.aligns.iterkeys():
                newaligns[curr] += trnaalign.aligns[curr][currpos]
            newstruct += trnaalign.currstruct[currpos]
            newconsensus += trnaalign.consensus[currpos]
            alnstring += '-'
            currpos += 1
        else:
            print >>sys.stderr, "nonstandard CIGAR type"
    while currpos < len(trnaalign.aligns[chrom]):
        for curr in trnaalign.aligns.iterkeys():
            newaligns[curr] += trnaalign.aligns[curr][currpos]
        newstruct += trnaalign.currstruct[currpos]
        newconsensus += trnaalign.consensus[currpos]
        alnstring += '.'
        currpos += 1
        #print >>sys.stderr, alnstring
        

    #print >>sys.stderr, trnaalign.aligns[chrom]
    #print >>sys.stderr, trnaalign.currstruct
    newseqs = newaligns
    newseqs[read.qname] = alnstring
    #print >>sys.stderr, alnstring
    return RnaAlignment(newseqs, newstruct, newconsensus)
'''        
old gapless form
def addreadalign(trnaalign, chrom, read):
    readstart = read.pos
    readend = read.aend
    

    alignstart = getalignposition(chrom,readstart, trnaalign)
    alignend = getalignposition(chrom,readend, trnaalign)
    
    #print "-"*alignstart + read.seq + "-"*(trnaalign.alignlength - alignend)
    alnstring = "."*alignstart
    currpos = alignstart
    for currbase in read.seq.replace("T","U"):
        
        while trnaalign.aligns[chrom][currpos] in gapchars:
            alnstring += trnaalign.aligns[chrom][currpos]
            currpos += 1

        alnstring += currbase
        currpos += 1
    alnstring =  alnstring +"."*(trnaalign.alignlength - currpos)
    #print trnaalign.aligns[chrom]
    #print "**********"
    newseqs = trnaalign.aligns.copy()
    newseqs[read.qname] = alnstring
    return RnaAlignment(newseqs, trnaalign.currstruct, trnaalign.consensus)
'''    
def oldgetalignposition(trnaname, position, stkalign):         
    seqpos = 0
    i = 0
    result = None
    resultpos = None
    
    for alignpos, curr in enumerate(stkalign.aligns[trnaname]): #for each place in the alignment
        if curr in gapchars:
            pass
        else:
            if i == position: #if it's reached the number of bases in the read
                #result = trnanums[alignpos]
                resultpos = alignpos
                trnabase = curr
            i += 1
    if resultpos is None:
        resultpos = alignpos + (position - i)
    return resultpos


def getalignposition(trnaname, position, stkalign):         
    seqpos = 0
    i = 0
    result = None
    resultpos = None
    #print >>sys.stderr, position
    #print >>sys.stderr, stkalign.aligns[trnaname]
    for alignpos, curr in enumerate(stkalign.aligns[trnaname]): #for each place in the alignment
        if curr in gapchars:
            pass
        else:
            if i == position: #if it's reached the number of bases in the read
                #result = trnanums[alignpos]
                resultpos = alignpos
                trnabase = curr
            i += 1
    if resultpos is None:
        resultpos = alignpos + (position - i)
        #print >>sys.stderr, "|(|**"

    return resultpos
    
def gettrnaposition(trnaname, position, stkalign, trnanums, returnpos = 'actual', getbase = False, includemargins = False, start = 0, margin = 0, reverse = False):
    
    seqpos = 0
    i = 0                                                    
    result = None
    resultpos = None
    
    #for alignpos, curr in enumerate(stkalign.aligns[trnaname]):
    #    if curr in gapchars:
    #        pass
    #    else:
    #        i += 1
    #        if i == position - start + margin:
    #            #result = trnanums[alignpos]
    #            resultpos = alignpos
    #            trnabase = curr
    ##print >>sys.stderr, resultpos                
    #if position  - start + margin < 0:
    #    resultpos = position  - start + margin
    #if resultpos is None:
    #    #print >>sys.stderr, str(position - start)
    #    #print >>sys.stderr, str(i)
    #    resultpos = alignpos + (position - start - i)
    #print >>sys.stderr, trnanums[resultpos]
    #print >>sys.stderr, resultpos
    #print >>sys.stderr, position - start + margin
    #print >>sys.stderr, stkalign.currstruct
    
    #resultpos = getalignposition(trnaname, position, stkalign)
    if reverse:
        resultpos = getalignposition(trnaname,start - position + margin, stkalign )
        #print >>sys.stderr, resultpos
        #print >>sys.stderr, start - position + margin
        #print >>sys.stderr, trnanums[resultpos]

        
        #print >>sys.stderr, "|*|"
    else:
        resultpos = getalignposition(trnaname,position - start + margin, stkalign )
        #print >>sys.stderr, "||"

        
    #print >>sys.stderr, resultpos
    #print >>sys.stderr, position - start+ margin
    #
    #print >>sys.stderr, "**"

    #print >>sys.stderr, len(trnanums)
    if resultpos < 0:
        return "head" + str(margin- (resultpos))
    elif resultpos >= len(trnanums):
        
        return "tail" + str(margin+ (resultpos - len(trnanums)))
    if getbase:
        return trnabase
    if trnanums is None:
        return resultpos
    if returnpos == 'actual' and is_number(trnanums[resultpos]):
        return trnanums[resultpos]
    if ishead(trnanums[resultpos]) and includemargins:
        return trnanums[resultpos]
    elif ishead(trnanums[resultpos]):
        return "0"
    elif istail(trnanums[resultpos]) and includemargins:
        return trnanums[resultpos]
    elif istail(trnanums[resultpos]):
        return "76"

    if returnpos == 'fixed':
      currpos = resultpos
      additional = 0
      while(not is_number(trnanums[currpos]) and currpos >= 0 and currpos < len(trnanums)):

        currpos -= 1
        additional += 1

                
        if is_number(trnanums[currpos]):
            
            return trnanums[currpos]+":i"+str(additional)
        elif currpos <= 0:
            return  "0" #'Start'
        elif currpos < len(trnanums):
            return  "76" #'End'
        else:
            return "NA"
    else:
        currpos = resultpos
        #print >>sys.stderr, trnanums[currpos]
        #print >>sys.stderr, "||"
        while(not is_number(trnanums[currpos]) and currpos >= 0 and currpos < len(trnanums)):
            if returnpos == 'next':
                currpos += 1
            elif returnpos == 'prev':
                currpos -= 1
            
                
            else:
                #maybe throw an exception here later
                return None
            
        if is_number(trnanums[currpos]):
            return trnanums[currpos]
        elif currpos <= 0:
            return "0" #'Start'
        elif currpos < len(trnanums):
            return "76" #'End'
        else:
            return "NA"

'''
orgmode	arch
forcecca	False
trna_margin	3
loci_margin	3
'''
class trnadbinfo:
    def __init__(self, dbfilename):
        
        self.orgmode = "euk"
        self.trnamargin = 20
        self.locimargin = 100
        self.forcecca = False  
        if not os.path.exists(dbfilename):
            print >>sys.stderr, "No database file "+dbfilename+" exists, using defaults"
        else:
            dbfile = open(dbfilename, "r")
            for currline in dbfile:
                fields = currline.rstrip().split("\t")
                #print >>sys.stderr, fields
                if fields[0] == "trna_margin":
                    self.trnamargin = int(fields[1])
                if fields[0] == "orgmode":
                    self.orgmode = fields[1]
                if fields[0] == "loci_margin":
                    self.locimargin = int(fields[1])
                if fields[0] == "forcecca":
                    self.forcecca = fields[1]
            
def gettrnafrags( trnafile,trnabed,stkfile,locibed,locistkfile,genomefile, dbinfoname,stkout = None, locistkout = None):
    trnastk = list(readrnastk(open(stkfile, "r")))[0]
    locistk = list(readrnastk(open(locistkfile, "r")))[0]                                                  
    
    dbinfo = trnadbinfo(dbinfoname)
    
    edgemargin = dbinfo.trnamargin
    #locusmargin = 100
    locimargin = dbinfo.locimargin
    trnamargin = dbinfo.trnamargin
    
    
    orgtype = dbinfo.orgmode

    positionnums = gettrnanums(trnastk, margin = trnamargin, orgtype = orgtype)
    locuspositionnums = gettrnanums(locistk, margin  = locimargin, orgtype = orgtype)
    #print positionnums
    #trnastk = trnastk.addmargin(edgemargin)
    # locistk = locistk.addmargin(locusmargin) 
    #print >>sys.stderr, "::"
    #print >>sys.stderr, trnamargin
    #print >>sys.stderr, positionnums
    trnadata = transcriptfile(trnafile)
    trnalocs = getnamedict(readbed(trnabed,orgdb = "genome",seqfile=genomefile))
    pretrnalocs = getnamedict(readbed(locibed,orgdb = "genome",seqfile=genomefile))
    trnatranscripts = set(trnadata.gettranscripts())
    trnaloci = set(trnadata.getloci())
    readtargets = set()

    trnaalign = trnastk.addmargin(trnamargin)
    
    if os.path.getmtime(genomefile+".fai") < os.path.getmtime(genomefile): 
        print >>sys.stderr, "fasta index file "+genomefile+".fai is older than fasta file and may not match"
        print >>sys.stderr, 'use command "samtools faidx '+genomefile+'.fai" to fix'
        sys.exit(1)
    lociflanks = getseqdict(list(curr.addmargin(locimargin) for curr in readbed(locibed,orgdb = "genome",seqfile=genomefile)), faifiles = {"genome":genomefile+".fai"})
    locistarts = {curr: lociflanks[curr][0:locimargin].replace("T","U") for curr in lociflanks}
    lociends = {curr: lociflanks[curr][-locimargin:].replace("T","U") for curr in lociflanks}
    locialign = locistk.addupstream(locistarts)
    locialign = locialign.adddownstream(lociends)
    #checking if grabbing the sequence failed
    #print >>sys.stderr, min(len(curr) for curr in lociflanks.values())
    if min(len(curr) for curr in lociflanks.values()) < locimargin:
        print >>sys.stderr, "Grabbing flanking sequence failed, cannot continue"
        print >>sys.stderr, "Check if fasta file"+genomefile+"is present and if index file "+genomefile+".fai matches it"
        sys.exit(1)
    #print >>sys.stderr, genomefile
    transcriptseqalign = trnaalign
    lociseqalign = locialign
    #locialign = locistk.addmargin(locimargin)
    
    hitscores = dict()
    curreadname = None
    readlength = None
    totalmatch = re.compile(r"^(?:(?P<startclip>\d+)S)?(?P<matchlength>\d+)M(?:(?P<endclip>\d+)S)?$")
    totalreads = 0
    skippedfrags = 0
    bamfile = pysam.Samfile("-", "r" )
    acsets = defaultdict(int)
    aminosets = defaultdict(int)
    
    print "\t".join(header)
    for pairedname, allmaps in itertools.groupby(bamfile,lambda x: x.qname):

        trnarange = "NA"
        allmaps = list(allmaps)
        #if unmapped
        if sum(curr.flag & 0x004 > 0 for curr in allmaps):
            #print pairedname +"\t"+"No_match"
            skippedfrags += 1
            continue
        totalreads += 1
        #print >>sys.stderr, "**"+pairedname
        readlength = None
        hitscores = dict()
        readtargets = set()
        
        mappings = 0
        currscore = None
        newset = set()
        readlength = None
        allaminos = set()
        allanticodons = set()
        #iterate through all mappings of the current read

        #if locistkfile is None:
        #    allmaps =  set(curr for curr in trnatranscripts for curr in allmaps if bamfile.getrname(curr.tid) not in trnaloci)
            
        newset = getbestmappings(allmaps, bamfile)

        if sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset) > 0 and sum(bamfile.getrname(curr.tid) in trnaloci for curr in newset) > 0:
            newset =  set(curr for curr in trnatranscripts for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
        transcriptseqalign = processbatch(newset, bamfile, transcriptseqalign,trnalocs, stkmargin = trnamargin)
        lociseqalign = processbatch(newset, bamfile, lociseqalign,pretrnalocs, stkmargin = locimargin)
        #print >>sys.stderr, len(lociseqalign.aligns.keys())
        trnastart = None
        trnaend = None
        #print >>sys.stderr, "||**"+",".join(bamfile.getrname(curr.tid) for curr in newset)
        #print >>sys.stderr, "||**"+",".join(trnaloci)
        
        if sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset) > 0:
            

            newset = list(curr for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            diff = len(newset) - sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset)
            anticodons = frozenset(trnadata.getanticodon(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            aminos = frozenset(trnadata.getamino(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            trnamappings = list(curr for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            trnanames = list(bamfile.getrname(curr.tid) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            #print >>sys.stderr, trnanames
            #print >>sys.stderr, anticodons
            #print >>sys.stderr, aminos
            trnaname = "tdR"
            fragtypes = set()
            mismatch = set()
            indel = set()
            mismatchpos = None
            mismatchlocs = set()
            for currread in newset:
                mismatch = set()
                chromname = bamfile.getrname(currread.tid)
                readlength = len(currread.seq)
                readstart = currread.pos
                readend = currread.aend
                currindel = 0
                for currtag in currread.tags:
                    if currtag[0] == "XM":
                        mismatch.add(currtag[1])
                    if currtag[0] == "XO":
                        #currindel += int(currtag[1])
                        pass
                    if currtag[0] == "XG":
                        currindel += int(currtag[1])
                indel.add(str(currindel))
                mismatchpos = readmdtag(currread)
                
                
                #trnastart = gettrnaposition(chromname, readstart+ 1, trnaalign, positionnums, returnpos = 'next', margin = trnamargin, start = trnalocs[chromname].start)  #dealing with sprinzel +1
                
                #print >>sys.stderr, readstart
                #trnaend = gettrnaposition(chromname, readend, trnaalign, positionnums, returnpos = 'prev',margin = trnamargin,  start = trnalocs[chromname].start)
                #varlooplen = varloopinfo(chromname, readstart+1, readend, trnastk, positionnums)
                #trnarange = trnastart+".."+trnaend
                #hasvarloop(chromname, trnastk)







                trnanums = gettrnanums(transcriptseqalign, margin = 0, orgtype = orgtype)
                if transcriptseqalign.aligns[pairedname][trnamargin] in gapchars:
                    firstpos = "1"
                else:
                    firstpos = "-1"
                finalpos = trnanums[-trnamargin - 1]
                
                for i in range(0,trnamargin):
                    #old naming
                    #trnanums[i] = "-"+str(trnamargin - i)
                    #trnanums[-(i+1)] = "+"+str(trnamargin - i)
                    
                        
                    trnanums[i] = "i"+str(trnamargin - i)+str(firstpos)
                    trnanums[-(i+1)] = str(finalpos)+"i"+str(trnamargin - i)
                
                
                #print >>sys.stderr, transcriptseqalign.consensus
                currdomalign = domainalign(trnanums)
                #print >>sys.stderr, currdomalign
                #print >>sys.stderr, "end"
                #print >>sys.stderr, len(trnanums)
                #print >>sys.stderr, trnaalign.alignlength
                #if pairedname == "frag10503:12":
                #bash testcommands.bash testfrag10503  alkb2-frags.fa /projects/lowelab/users/holmes/pythonsource/tdrtest/hg19/hg19
              
                trnastart = None
                trnaend = None
                domstart = None
                varlooplen = 0
                alldoms =  set()
                haszeroposition = False
                for i in range(transcriptseqalign.alignlength):
                    if str(trnanums[i]) == "0":
                        haszeroposition = True
                    if transcriptseqalign.aligns[pairedname][i] not in gapchars: #and (is_number(trnanums[i]) or trnanums[i].startswith("e") or trnanums[i] in sprinzeladditional):
                        if trnastart is None:
                            trnastart = trnanums[i]
                        trnaend = trnanums[i]
                        if currdomalign[i] is not None:
                            alldoms.add(currdomalign[i])
                    if  transcriptseqalign.aligns[pairedname][i] not in gapchars and trnanums[i].startswith("e"):
                        varlooplen += 1
                    if transcriptseqalign.aligns[pairedname][i] in gapchars or transcriptseqalign.aligns[chromname][i] in gapchars:
                        continue
                    if transcriptseqalign.aligns[pairedname][i].upper() != transcriptseqalign.aligns[chromname][i].upper():
                        mismatchlocs.add(transcriptseqalign.aligns[pairedname][i] + ":"+trnanums[i] + ':'+ transcriptseqalign.aligns[chromname][i])
                #special for handling histidines
                if str(trnastart) == "0":
                    trnastart = "-1"
                    
                elif haszeroposition and trnastart.startswith("-"):
                    newpos = int(trnastart[1:])
                    trnastart = "-" + str(newpos + 1)
                    #print >>sys.stderr, trnanums
                    #print >>sys.stderr, currdomalign
                        
                #this is going to get triggered in long gaps        
                if trnastart is None:
                    #print >>sys.stderr, "**"+pairedname
                    trnastart = "1"
                if trnaend is None:
                    #print >>sys.stderr, "**"+pairedname
                    trnaend = "76"
                if trnastart is None or trnaend is None:
                    print >>sys.stderr, "Cannot find Sprinzel coordinates for "+pairedname
                    sys.exit(1)
                trnarange = trnastart+".."+trnaend
                #print >>sys.stderr, alldoms

                    
                armdomains = getdomainlist(alldoms)
                #print >>sys.stderr, pairedname+"**"+str(varlooplen)
                
                for currpos in mismatchpos.iterkeys():
                    pass
                    #mismatchpos[currpos]
                    #mismatchlocs.add(gettrnaposition(chromname, currpos - 1, trnaalign, positionnums, getbase = True)+":"+gettrnaposition(chromname, currpos - 1, trnaalign, positionnums) +':'+ mismatchpos[currpos])
                    
                    #print >>sys.stderr, findtag(currread, "MD")
                    #print >>sys.stderr, currpos
                    #print >>sys.stderr, currread.seq
                    #trnanums = gettrnanums(trnaalign, margin = edgemargin)
                    #mismatchlocs.add(trnaalign.aligns[pairedname][getalignposition(chromname,readstart, trnaalign)]+":"+gettrnaposition(chromname, currpos - 1, trnaalign, positionnums) +':'+ mismatchpos[currpos])
                    
                currtrna = trnalocs[chromname]
                if currtrna.start - 3 < readstart <  currtrna.start + 3 and currtrna.end - 5 < readend <  currtrna.end +3:
                    fragtypes.add("W")
                elif currtrna.start - 3 < readstart <  currtrna.start + 3:
                    fragtypes.add("5P")
                elif currtrna.end - 5 < readend <  currtrna.end + 3:
                    fragtypes.add("3P")
                else:
                    fragtypes.add("OF")
            if len(fragtypes) == 1:
                fragtype = list(fragtypes)[0]
            else:
                fragtype = "OF"
            if len(aminos - frozenset(['Und'])) > 1:
                trnatype = "Mult"
            elif len(anticodons - frozenset(['NNN'])) > 1:
                trnatype = list(aminos)[0]
            elif len(trnamappings) > 1:
                trnatype = list(aminos)[0]+"-"+list(anticodons)[0]
            else:
                trnatype = list(aminos)[0]+"-"+list(anticodons)[0] +"-"+ list(trnanames)[0].split("-")[-1] #+"-"+"uniq"

            #tags 
            tags = [("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            for currtrnamap in trnamappings:
                currtrnamap.tags = currtrnamap.tags + [("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            finalset = trnamappings
            if len(mismatchlocs) < 1:
                mismatchlocs.add("None")
            trnaname = "tdR-"+fragtype+"-"+trnatype
            #print pairedname +"\ttRNA\t"+trnaname+"\t"+",".join(trnanames)+"\t"+"\t"+",".join(anticodons)+"\t"+",".join(aminos)+"\t"+",".join(str(curr) for curr in mismatch) +"\t"+trnarange+ "\t"+",".join(mismatchlocs)
            #print >>sys.stderr, trnastart
            #print >>sys.stderr, trnaend
            

            
            #armdomains = getarms(trnastart,trnaend,varlooplen)
            #print >>sys.stderr, pairedname+str(varlooplen)
            #trnadomains = list(newdomabbv[currdom] for currdom in getdomains(trnastart,trnaend,varlooplen))
            #
            #domrange = list()
            #domrange.append(trnadomains[0] +'_'+ trnadomains[-1])
            #domrange.append("_".join(trnadomains))
            #domrange.append(str(trnastart)+":"+trnadomains[0] +'_'+ trnadomains[-1] +":"+str(trnaend))
            if list(fragtypes)[0] == 'W':
                trnaname = list(trnanames)[0]
            #list(["tRFdb_ID","source_type",'tdR_name','source_trnas','length','Sprinzl_range','anticodons','isotypes','mismatches','indels','mismatch_locations','trna_domains'])
            


            if len(armdomains) == 0:
                armdomains.append("None")
                
            print "\t".join([pairedname,"tRNA",trnaname,",".join(trnasort(trnanames)),str(readlength),trnarange,",".join(anticodons),",".join(aminos),",".join(str(curr) for curr in mismatch),",".join(str(curr) for curr in indel),",".join(mismatchlocs),"_".join(armdomains)])
            if len(mismatchlocs) != int(",".join(str(curr) for curr in mismatch)):
                #print >>sys.stderr, pairedname +" mismatch count does not match"
                pass
        #sum(bamfile.getrname(curr.tid) in trnaloci for curr in newset) > 0 or True:
        #check for loci
        elif locistkout is not None:
            currread = list(newset)[0]
            #here is where I allow for either db type
            chromname = bamfile.getrname(list(newset)[0].tid)
            locusnames = list()
            for currread in  newset:
                chromname = bamfile.getrname(currread.tid)
                if chromname in pretrnalocs:
                    locusnames.append(chromname)
                else:
                    for currpretrna in pretrnalocs.iterkeys():
                        if readcoverage(pretrnalocs[currpretrna], currread,chromname) > 10:
                            locusnames.append(currpretrna)
                            break
            #trnanames = set(bamfile.getrname(currread.tid) for currread in newset)
            if len(locusnames) < 1:
                continue
            if len(locusnames) > 1:
                #print >>sys.stderr, "Multiple pre-tRNA matches to "+pairedname +" "+",".join(locusnames)
                pass
                
            locusname = locusnames[0]
            trnanames = locusnames
            if pretrnalocs[locusname].strand == '+':
                pass
                #continue
            transcripts = list(trnadata.getlocustranscript(locusname) for curr in newset)
            anticodons = frozenset(trnadata.getanticodon(curr) for curr in transcripts)
            aminos = frozenset(trnadata.getamino(curr) for curr in transcripts)
            mismatchlocs = set(["NA"])
            mismatch = set()
            currindel = 0
            indel = set()
            
            for currtag in currread.tags:
                
                if currtag[0] == "XM":
                    mismatch.add(currtag[1])
                if currtag[0] == "XO":
                    pass
                    #currindel += int(currtag[1])
                if currtag[0] == "XG":
                    currindel += int(currtag[1])
            indel.add(str(currindel))
            chromname = list(newset)[0]
            currtrna = pretrnalocs[locusname]
            readstart = currread.pos
            readend = currread.aend
            #print >>sys.stderr, "||***"
            if readend < currtrna.start - 5 or readstart > currtrna.end + 5:
                continue
                skippedfrags += 1
            elif currtrna.end - 2 < readstart < currtrna.end + 2:
                fragtype = "T"
            elif currtrna.start - 2 < readend < currtrna.start + 2:
                fragtype = "L"
            elif readstart < currtrna.start  and readend >  currtrna.end:
                fragtype = "LWT"
            elif readstart < currtrna.start:
                fragtype = "L5P"
            elif readend >  currtrna.end:
                fragtype = "3PT"
            else:
                fragtype = "I"
            if len(aminos - frozenset(['Und'])) > 1:
                trnatype = "Mult"
            elif len(anticodons - frozenset(['NNN'])) > 1:
                trnatype = list(aminos)[0]
            elif len(set(transcripts)) > 1:
                trnatype = list(aminos)[0]+"-"+list(anticodons)[0]
            elif len(trnanames) > 1:
                trnatype = list(aminos)[0]+"-"+list(anticodons)[0] +"-"+ list(trnanames)[0].split("-")[-2]
            else:
                trnatype = list(aminos)[0]+"-"+list(anticodons)[0] +"-"+ list(trnanames)[0].split("-")[-2] +"-"+list(trnanames)[0].split("-")[-1] #+"-"+"uniq"
            #locusre = re.compile(r'.*tRNA\-(.*)')
            #locusmatch = locusre.search(locusname)
            #print >>sys.stderr, "||***"
            #if locusmatch:
            #    shortlocusname = locusmatch.group(1)
            #else:
            #    print >>sys.stderr, "Cannot find tRNA name for"+ locusname
            #    sys.exit(1)
            
            #print >>sys.stderr, readstart - pretrnalocs[locusname].start
            
            #gotta deal with reverse strand ugh
            #print >>sys.stderr, "||"
            #
            #print >>sys.stderr,pretrnalocs[locusname].strand  
            #print >>sys.stderr,  pretrnalocs[locusname].start
            #print >>sys.stderr,  pretrnalocs[locusname].end
            #print >>sys.stderr,  readstart+1
            #print >>sys.stderr,  locuspositionnums
            mismatchlocs = set()
            trnastart = None
            trnaend = None
            #print >>sys.stderr, lociseqalign.aligns[pairedname]
            
            locitrnanums = gettrnanums(lociseqalign, margin = 0, orgtype = orgtype)
            locitrnanums = replacepadding(locitrnanums, margin = locimargin)
            
            
            #print >>sys.stderr, str(len(locuspositionnums)) + "**"
            #print >>sys.stderr, lociseqalign.aligns[pairedname]
            #print >>sys.stderr, lociseqalign.currstruct
            #print >>sys.stderr, "".join("*" if curr.startswith("head") or curr.startswith("tail") else "-" for curr in locitrnanums)
            #print >>sys.stderr, lociseqalign.consensus
            #
            #sys.exit()

            
            for i in range(lociseqalign.alignlength): #frag5569:22
                if  i > len(locuspositionnums):
                    #print >>sys.stderr, pairedname + "**"
                    #print >>sys.stderr, str(len(locitrnanums)) + "**"
                    #print >>sys.stderr, str(len(locitrnanums)) + "**"
                    #
                    #print >>sys.stderr, str(lociseqalign.alignlength) + "**"
                    #print >>sys.stderr, locitrnanums[i] + "**"
                    pass
                elif lociseqalign.aligns[locusname][i]:
                    #print >>sys.stderr, pairedname + "||"
                    pass
                
                if lociseqalign.aligns[pairedname][i] in gapchars or lociseqalign.aligns[locusname][i] in gapchars: # or not (is_number(locitrnanums[i]) or locitrnanums[i].startswith('head') or locitrnanums[i].startswith('tail')):
                    continue
                if locitrnanums[i].startswith('head') or locitrnanums[i].startswith('tail'):
                    #print >>sys.stderr, "**"
                    pass
                if lociseqalign.aligns[pairedname][i] not in gapchars:# and (is_number(locitrnanums[i]) or locitrnanums[i].startswith('head') or locitrnanums[i].startswith('tail')):
                    if trnastart is None:
                        trnastart = locitrnanums[i]
                    trnaend = locitrnanums[i]
                if lociseqalign.aligns[pairedname][i].upper() != lociseqalign.aligns[locusname][i].upper():
                    mismatchlocs.add(lociseqalign.aligns[pairedname][i] + ":"+locitrnanums[i].replace("head","-").replace("tail","+") + ':'+ lociseqalign.aligns[locusname][i])
            if trnastart is None or trnaend is None:
                print >>sys.stderr, str(len(locuspositionnums)) + "**|*|"
                print >>sys.stderr, lociseqalign.aligns[pairedname]
                print >>sys.stderr, lociseqalign.currstruct
                print >>sys.stderr, "".join("*" if curr.startswith("head") or curr.startswith("tail") else "-" for curr in locitrnanums)
                print >>sys.stderr, trnastart
                print >>sys.stderr, trnaend
            if len(mismatchlocs) < 1:
                mismatchlocs.add("None")
                
                
                
            '''
            if pretrnalocs[locusname].strand == '-':
                #print >>sys.stderr, "||||" 
                
                trnastart = gettrnaposition(locusname, readend, locialign, locuspositionnums, returnpos = 'prev',start = pretrnalocs[locusname].end, margin = locimargin,includemargins = True, reverse = pretrnalocs[locusname].strand == '-')            
                trnaend = gettrnaposition(locusname,readstart+1  , locialign, locuspositionnums, returnpos = 'next',  start = pretrnalocs[locusname].end, margin = locimargin, includemargins = True,reverse = pretrnalocs[locusname].strand == '-')

            else:
                trnastart = gettrnaposition(locusname, readstart+1, locialign, locuspositionnums, returnpos = 'prev',start = pretrnalocs[locusname].start, margin = locimargin,includemargins = True, reverse = pretrnalocs[locusname].strand == '-')            
                trnaend = gettrnaposition(locusname, readend, locialign, locuspositionnums, returnpos = 'next',  start = pretrnalocs[locusname].start, margin = locimargin, includemargins = True,reverse = pretrnalocs[locusname].strand == '-')
            '''
            #varlooplen = varloopinfo(chromname, readstart, readend, trnastk, positionnums)
            if trnastart is None or trnaend is None:
                print >>sys.stderr, pairedname
                sys.exit(1)
            trnastart = trnastart.replace("head","-")
            trnastart = trnastart.replace("tail","+")
            trnaend = trnaend.replace("head","-")
            trnaend = trnaend.replace("tail","+")
            trnarange = trnastart+".."+trnaend
            ##########
            
            #trnaname = "tdR-"+fragtype+"-"+shortlocusname
            #trnaname = "tdR-"+fragtype+"-"+shortlocusname
            trnaname = "tdR-"+fragtype+'-'+trnatype
            readlength = len(currread.seq)
            #list(["tRFdb_ID","source_type",'tdR_name','source_trnas','length','Sprinzl_range','anticodons','isotypes','mismatches','indels','mismatch_locations','trna_domains'])
            print "\t".join([pairedname,"pretRNA",trnaname,",".join(trnasort(trnanames)),str(readlength),trnarange,",".join(anticodons),",".join(aminos),",".join(str(curr) for curr in mismatch),",".join(str(curr) for curr in indel),",".join(mismatchlocs),'None'])
            #print ",".join(transcripts)
            pass
    #print >>sys.stderr, "skipped: "+str(skippedfrags)
    if stkout is not None:
        stkfile = open(stkout, "w")
        
        transcriptseqalign.pruneNbases().printstk(output = stkfile,ordering= trnasort(transcriptseqalign.aligns.keys()))
    if locistkout is not None:
        stkfile = open(locistkout, "w")
        lociseqalign.printstk(output = stkfile, ordering= trnasort(lociseqalign.aligns.keys()))


#nametdrs.py /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-trnatable.txt /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-maturetRNAs.bed trnadbs/hg19-trnaalign.stk
 #bowtie2 -f -x /projects/lowelab/users/holmes/pythonsource/tdrnamer/hgdbtest/hgtest-tRNAgenome -k 100 --very-sensitive --ignore-quals --np 5 --reorder -p 4 -U testfrags.fa 2>bowtieerr.txt | nametdrs.py hgdbtest/hgtest
dbname = sys.argv[1]
stkout = sys.argv[2]
locistkout = None
if len(sys.argv) > 3:
    locistkout = sys.argv[3]
#bowtie2 -f -x /projects/lowelab/users/holmes/pythonsource/tdrnamer/trnadbs/hg19-tRNAs -k 100 --very-sensitive --ignore-quals --np 5 --reorder -p 4 -U testfrags.fa 2>bowtieerr.txt | nametdrs.py hgdbtest/hgtest
trnatable = dbname + '-trnatable.txt'
trnabed = dbname + '-maturetRNAs.bed'
locibed = dbname + '-trnaloci.bed'
stkfile = dbname + '-trnaalign.stk'
#stkfile = dbname + '-trnaconvert.stk'
locistk = dbname + '-trnaloci.stk'
sequences = dbname + '-tRNAgenome.fa'

dbinfo = dbname + '-dbinfo.txt'

#locistkout = None
gettrnafrags( trnatable,trnabed,stkfile,locibed, locistk ,sequences,dbinfo, stkout = stkout, locistkout = locistkout)
