#!/usr/bin/env python


import re
import os

from collections import defaultdict
import argparse

import time
import sys
import argparse
import os.path
from collections import defaultdict
from tdrdbutils import *
import itertools
import subprocess
from distutils.spawn import find_executable



#This program gets the mature tRNA sequences

#need to move this somewhere else later
allintronseqs = defaultdict(list)
locustranscripts = dict()
intronsplices = defaultdict(list)


def getgithash(scriptdir):
    gitloc = get_location("git")
    
    if gitloc is None:
        print >>sys.stderr, "Cannot find git in path"
        print >>sys.stderr, "Recording of versioning not possible"
    gitjob = subprocess.Popen([gitloc,"--git-dir="+scriptdir+"/.git","rev-parse","HEAD"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    githash = gitjob.communicate()[0].rstrip()
    if gitjob.returncode != 0:
        githash = "Cannot find git hash"
    gitjob = subprocess.Popen([gitloc,"--git-dir="+scriptdir+"/.git","describe"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    gitversion = gitjob.communicate()[0].rstrip()
    if gitjob.returncode != 0:
        gitversion = "Cannot find git version"
    return gitversion, githash 

class tRNAlocus:
    def __init__(self, loc, seq, score, amino, anticodon, intronseq, rawseq = None):
        self.name = loc.name
        self.loc = loc
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.intronseq = intronseq
        self.rawseq = rawseq

class tRNAtranscript:
    def __init__(self, seq, score, amino, anticodon, loci, intronseq, name = None, rawseq = None):
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.loci = tuple(loci)
        self.intronseqs = intronseq
        self.name = name

        self.rawseq = rawseq
    def getmatureseq(self, prokmode = False):
        prefix = ""
        
        #print >>sys.stderr, self.amino
        if self.amino == "His":
            prefix = "G"
        end = "CCA"
        if prokmode:
            end = ""
        return prefix + self.seq + end 
    
            

def getuniquetRNAs(trnalist):
    sequencedict = defaultdict(list)
    scoredict = defaultdict(set)
    for curr in trnalist:
        sequencedict[curr.seq].append(curr)
    for currtrans in sequencedict.iterkeys():
        scores = set(curr.score for curr in sequencedict[currtrans])
        anticodon = set(curr.anticodon for curr in sequencedict[currtrans])
        amino = set(curr.amino for curr in sequencedict[currtrans])
        introns = set(curr.intronseq for curr in sequencedict[currtrans])
        #remove psuedo if there's a real one somewheres
        if len(anticodon) > 1:
            anticodon.discard('Xxx')
            amino.discard('X')
        if introns == set([""]):
            introns = set()
        loci = list(curr for curr in sequencedict[currtrans])
        if len(scores) > 1:
            #print >>sys.stderr, "Multiple scores"
            pass
        if len(anticodon) > 1:
            print >>sys.stderr, "tRNA file contains identical tRNAs with seperate anticodons, cannot continue"
            sys.exit()
        yield tRNAtranscript(currtrans, scores,list(amino)[0],list(anticodon)[0],loci, introns)
        
        
'''
def readrnacentral(scanfile,chromnames, mode = 'locus'):
    reftochrom = dict()
    convertfile = open(chromnames)
    
    for line in convertfile:
        if line.startswith("#"):
            continue
        fields = line.split()
        reftochrom[fields[1]] = fields[0]
    orgname = "genome"
   
    #print reftochrom
    trnafile = open(scanfile)
    transcriptinfo = defaultdict(list)
    for line in trnafile:
        fields = line.split(",")
        #print len(fields)
        #print "**"
        if fields[0] == 'Entry number':
            continue
        #print fields[1]
        if True:
            name = fields[2]
            amino = fields[11]
            anticodon = fields[12]
            sequence = fields[21]
            score = float(fields[18])
            #print >>sys.stderr, score
            #print name+":"+amino+":"+anticodon
            #print sequence
            gbchrom = re.sub(r'\.\d+$', '', fields[7])
            if gbchrom not in reftochrom:
                continue
            chrom = reftochrom[gbchrom]
            #print >>sys.stderr, chrom
            start = int(fields[8].split('-')[0]) - 1
            end = fields[8].split('-')[1]
            transcriptname = re.sub(r'-\d+$', '', fields[2])
            if fields[9] == "yes":
                strand = "-"
                
            elif fields[9] == "no":
                strand = "+"
            currtRNA = GenomeRange(orgname, chrom,start,end, name = name,strand = strand,orderstrand = True)
            #print >>sys.stderr, chrom
            intronnums = set()
            intronseqs = ""
            for currintron in fields[14:16]:
                intronmatch  = re.search(r'(\d+)\.\.(\d+)',currintron)
                if intronmatch:
                    intronstart = int(intronmatch.group(1))-1
                    intronend = int(intronmatch.group(2))
                    intronnums |= set(range(intronstart, intronend))
                    intronseqs = sequence[intronstart:intronend]
            
            newseq = ''
            for i in range(len(sequence)):
                if i in set(intronnums):
                    newseq += '-'
                else:
                    newseq += sequence[i]
                    
            rawseq = newseq.replace('-','')
            currlocus = tRNAlocus(currtRNA,rawseq, score,amino,anticodon,intronseqs)
            transcriptinfo[transcriptname].append(currlocus)
            if mode == 'locus':
                yield currlocus
        allseqs = dict()
    print >>sys.stderr, len(transcriptinfo.keys())
    for currtrans in transcriptinfo.iterkeys():
        if len(set(curr.seq for curr in transcriptinfo[currtrans])) > 1:
            print >>sys.stderr, "multiple"
        #print transcriptinfo[currtrans][0].seq
        if transcriptinfo[currtrans][0].seq in allseqs:
            print  >>sys.stderr, "duplicate:" + currtrans + ":"+allseqs[transcriptinfo[currtrans][0].seq]
        allseqs[transcriptinfo[currtrans][0].seq] = currtrans
        if mode == 'transcript':
            #print >>sys.stderr,currtrans 
            yield tRNAtranscript( transcriptinfo[currtrans][0].seq, set(curr.score for curr in transcriptinfo[currtrans]), transcriptinfo[currtrans][0].amino, transcriptinfo[currtrans][0].anticodon, set(transcriptinfo[currtrans]), transcriptinfo[currtrans][0].intronseq, name = currtrans)
        

def readtRNAscan(scanfile, genomefile, mode = None):
    #mode = 'gtRNAdb'

    trnalist = list()
    orgname = "genome"
    if  hasattr(scanfile ,'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile)
    trnascore = dict()
    trnaanticodon = dict()
    trnaamino = dict()
    tRNAintron = dict()
    trnas = dict()

    for currline in trnascan:
        
        if  currline.startswith("Sequence") or  currline.startswith("Name") or  currline.startswith("------"):
            continue
        fields = currline.split()
        #print >>sys.stderr, fields[0]
        if mode == "gtRNAdb":
            #print >>sys.stderr, fields[6:8]
            del fields[6:8]
        curramino = fields[4]
        currac = fields[5]
        
        if currac == "???":
            currac = 'Xxx'
        if fields[2] > fields[3]:
            end = int(fields[3]) - 1
            start = int(fields[2])
            
        else:
            end = int(fields[3])
            start = int(fields[2]) - 1
        currchrom = fields[0]
        trnanum = fields[1]
        currtRNA = GenomeRange(orgname, currchrom,start,end, name = currchrom+"."+"tRNA"+trnanum+"-"+curramino+currac,strand = "+",orderstrand = True)
        currtrans = currtRNA

        trnaamino[currtrans.name] = curramino
        trnaanticodon[currtrans.name] = currac
        #print >>sys.stderr, "**".join(fields)#currline
        trnascore[currtrans.name] =  float(fields[8])
        trnas[currtrans.name] =  currtrans
    
    
        
        currtRNA.fastafile = genomefile
        #print >>sys.stderr, len(currtrans.name)
        trnalist.append(currtRNA)
        if int(fields[6]) != 0:
            if currtRNA.strand ==  "-":
                intronstart = int(fields[2]) - int(fields[6]) 
                intronend = int(fields[2]) - int(fields[7]) +1
            else:
                intronstart = int(fields[6]) - int(fields[2]) - 1
                intronend = int(fields[7]) - int(fields[2])
            tRNAintron[currtRNA.name] = tuple([intronstart, intronend])
    trnaseqs = getseqdict(trnalist, faifiles = {orgname:genomefile+".fai"})
    intronseq = defaultdict(str)
    print >>sys.stderr, "***"
    for curr in trnaseqs.iterkeys():
        allintronseqs[curr] = ''
        print >>sys.stderr, curr
        if curr in tRNAintron:
            start = tRNAintron[curr][0]
            end = tRNAintron[curr][1]
            intronseq[curr] = trnaseqs[curr][start:end] 
            allintronseqs[curr] =  trnaseqs[curr][start-1:end-1].lower().replace("t","u")
            intronsplices[curr] = [start,end]
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
            if curr == "tRNA-Leu-CAA-1-2":
                #print >>sys.stderr, trnaseqs[curr]
                pass
        locustranscripts[transname] = set(locus.name for locus in currlocus)
        
        yield tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr],trnaamino[curr],trnaanticodon[curr],intronseq[curr])
        
'''
def striplocus(trnaname):
    return re.sub(r"\-\d+$", "",trnaname)
#Sequence		tRNA     	Bounds   	tRNA	Anti	Intron Bounds	Inf	HMM	2'Str	Hit	      	Isotype	Isotype	Type
#Name    	tRNA #	Begin    	End      	Type	Codon	Begin	End	Score	Score	Score	Origin	Pseudo	CM	Score	
#--------



def padintrons(intronseqs):
    
    maxsize = max(len(intronseqs[currname]) for currname in intronseqs.iterkeys())
    intronal = dict()
    for currname in intronseqs.iterkeys():
        #for i in range(len(intronseqs[currname])):
        intronal[currname] = intronseqs[currname] + ("." * (maxsize - len(intronseqs[currname])))
    consensus = "." * maxsize
    structure = "." * maxsize
    
    #print >>sys.stderr, intronal
    return RnaAlignment(intronal,structure , consensus = consensus)
    
#=======++====+++++++++++====+=====++***+
def oldconvertalignment(trnaalignment, intronseqs, locusdict):
    alignment =  list(readrnastk(open(trnaalignment, "r")))[0]
    
    startcolumns = list()
    endcolumns = list()
    realcols = list()
    intronstart = None
    currpos = 0
    for i, currcons in enumerate(alignment.consensus):
        if currcons == '+' and len(startcolumns) == 0:
            realcols.append(i)
        if currcons in set('+=*'):
            if currpos == 40:
                intronstart = i
            currpos += 1
            
    for currname in intronseqs.iterkeys():
        for i in range(len(intronseqs[currname])):
            intronseqs[currname][i] = intronseqs[currname][i].lower().replace("t","u")
            
    endcolumns = range(realcols[-3]-1,len(alignment.consensus))
    startcolumns = range(realcols[0] + 1)
    #print >>sys.stderr, alignment.consensus
    alignment = alignment.trimcolumns(startcolumns)
    alignment = alignment.trimcolumns(endcolumns)
    
    intronstart = intronstart - realcols[0]
    intronalignment = padintrons(intronseqs)
    
    locusalignment = alignment.convertnames(locusdict)
    
    #print >>sys.stderr, set(locusalignment.aligns.keys()) - set(intronalignment.aligns.keys())
    #print >>sys.stderr, set(locusalignment.aligns.keys()) 
    #print >>sys.stderr, set(intronalignment.aligns.keys())
    
    locusalignment = locusalignment.insertcolumns(intronstart, intronalignment)

    
    return locusalignment
    
#still does not handle non-canonical introns correctly    
def convertalignment(trnaalignment, intronseqs, locusdict):
    alignment =  list(readrnastk(open(trnaalignment, "r")))[0]
    
    startcolumns = list()
    endcolumns = list()
    realcols = list()
    intronstart = None
    currpos = 0
    
    for i, currcons in enumerate(alignment.consensus):
        if currcons == '+' and len(startcolumns) == 0:
            realcols.append(i)
        if currcons in set('+=*'):
            if currpos == 40:
                intronstart = i
            currpos += 1
            
    #for currname in intronseqs.iterkeys():
    #    for i in range(len(intronseqs[currname])):
    #        intronseqs[currname][i] = intronseqs[currname][i].lower().replace("t","u")
            
    endcolumns = range(realcols[-3]-1,len(alignment.consensus))
    startcolumns = range(realcols[0] + 1)
    #print >>sys.stderr, alignment.consensus
    alignment = alignment.trimcolumns(startcolumns)
    alignment = alignment.trimcolumns(endcolumns)
    
    intronstart = intronstart - realcols[0]
    #print >>sys.stderr, intronseqs
    #this is not correct but only matters for multiintronic tRNAs
    introncat = defaultdict(str)
    #print >>sys.stderr, intronseqs
    
    for currname in itertools.chain.from_iterable(locusdict.itervalues()):
        introncat[currname] = "".join(intronseqs[currname])
    #print >>sys.stderr, introncat    tRNA-Ser-TGA-1-1
    #print >>sys.stderr, introncat["tRNA-Ser-TGA-1-1"]
    #print >>sys.stderr, "**"
    #print >>sys.stderr, introncat
    intronalignment = padintrons(introncat)
    #print >>sys.stderr, "**||"
    #print >>sys.stderr, intronalignment.aligns
    locusalignment = alignment.convertnames(locusdict)
    
    #print >>sys.stderr, set(locusalignment.aligns.keys()) - set(intronalignment.aligns.keys())
    #print >>sys.stderr, set(locusalignment.aligns.keys()) 
    #print >>sys.stderr, intronalignment.aligns.keys()
    
    locusalignment = locusalignment.insertcolumns(intronstart, intronalignment)

    return locusalignment 
    
def splicealignment(trnaalignment, intronlocs, locusdict, prokmode = False):
    alignment =  list(readrnastk(open(trnaalignment, "r")))[0]
    #print >>sys.stderr, intronlocs["tRNA-Leu-CAA-1-2"]
    #print >>sys.stderr, intronlocs
    #locusalignment = alignment.splice(intronlocs)
    locusalignment = alignment.multisplice(intronlocs)
    transdict = dict()
    for currtranscript in locusdict.iterkeys():
        currlocus = list(locusdict[currtranscript])[0]
        transdict[currlocus] = [currtranscript]
        if currtranscript == 'tRNA-Leu-CAA-1':
            #print >>sys.stderr, currlocus #tRNA-Leu-CAA-1-2
            #sys.exit(1)
            pass
        
            
    #print >>sys.stderr, transdict 
    
    transalignment = locusalignment.convertnames(transdict)
    
    #print >>sys.stderr, set(locusalignment.aligns.keys()) - set(intronalignment.aligns.keys())
    #print >>sys.stderr, set(locusalignment.aligns.keys()) 
    #print >>sys.stderr, set(intronalignment.aligns.keys())
    

    downstream = {curr[0]:"CCA" for curr in transdict.values()}
    upstream = dict()
    for curr in transdict.values():
        if curr[0].split("-")[1] == "His":
            upstream[curr[0]] = "G"
        else:
            upstream[curr[0]] = "-"
            
    #upstream = {curr:"G"  else curr:"G" }
    transalignment = transalignment.addupstream(upstream, consensus = "+")
    if not prokmode:
        transalignment = transalignment.adddownstream(downstream, consensus = "+++")    
    return transalignment 
    
def readtRNAdb(scanfile, genomefile, trnamap):
    
    #mode = 'gtRNAdb'
    trnalist = list()
    orgname = "genome"
    if  hasattr(scanfile ,'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile)
    trnascore = dict()
    trnaanticodon = dict()
    trnaamino = dict()
    tRNAintron = defaultdict(list)
    trnas = dict()
    for currline in trnascan:
        if  currline.startswith("Sequence") or  currline.startswith("Name") or  currline.startswith("------"):
            continue
        fields = currline.split()
        
        curramino = fields[4]
        currac = fields[5]
        
        if currac == "???":
            pass
            #currac = 'Xxx'
        if fields[2] > fields[3]:
            end = int(fields[3]) - 1
            start = int(fields[2])
            
        else:
            end = int(fields[3])
            start = int(fields[2]) - 1
        currchrom = fields[0]
        trnanum = fields[1]
        
        trnascanname = currchrom+"."+"trna"+trnanum+"-"+curramino+currac
        #print >>sys.stderr, trnamap.keys()
        #print >>sys.stderr, trnascanname
        shorttrnascanname = currchrom+"."+"trna"+trnanum

        
        if trnascanname in trnamap:
            currtRNA = GenomeRange(orgname, currchrom,start,end, name = trnamap[trnascanname],strand = "+",orderstrand = True)
        elif shorttrnascanname in trnamap:
            currtRNA = GenomeRange(orgname, currchrom,start,end, name = trnamap[shorttrnascanname],strand = "+",orderstrand = True)
        else:
            #print >>sys.stderr, "Skipping "+trnascanname+", has no transcript name"
            continue
        currtrans = currtRNA

        trnaamino[currtrans.name] = curramino
        trnaanticodon[currtrans.name] = currac
        #print >>sys.stderr, "**".join(fields)#currline
        trnascore[currtrans.name] =  float(fields[8])
        trnas[currtrans.name] =  currtrans
    
    
        
        currtRNA.fastafile = genomefile
        trnalist.append(currtRNA)
        
        startintronnums = fields[6].split(",")
        endintronnums = fields[7].split(",")
        for i in range(len(startintronnums)):
            if int(startintronnums[i]) != 0:
                if currtRNA.strand ==  "-":
                    intronstart = int(fields[2]) - int(startintronnums[i])  
                    intronend = int(fields[2]) - int(endintronnums[i]) + 1
                else:
                    intronstart = int(startintronnums[i]) - int(fields[2]) 
                    intronend = int(endintronnums[i]) - int(fields[2]) + 1
                tRNAintron[currtRNA.name].append(tuple([intronstart, intronend]))
    trnaseqs = getseqdict(trnalist, faifiles = {orgname:genomefile+".fai"})
    intronseq = defaultdict(list)
    trnaloci = list()
    
    
    for curr in trnaseqs.iterkeys():
        #allintronseqs[curr] = ''
        #print >>sys.stderr, curr
        for currintron in reversed(tRNAintron[curr]):
            
            start = currintron[0]
            end = currintron[1]
            intronseq[curr].append(trnaseqs[curr][start:end])
            allintronseqs[curr].append( trnaseqs[curr][start:end].lower().replace("t","u"))
            intronsplices[curr].append([start,end])
            if curr in set(["tRNA-Leu-CAA-1-1","tRNA-Leu-CAA-1-2"]):
                #print >>sys.stderr, curr
                #print >>sys.stderr, [start,end]
                #print >>sys.stderr, trnaseqs[curr]
                #print >>sys.stderr, trnaseqs[curr][:start] + trnaseqs[curr][end:]
                pass
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
            
        trnaloci.append( tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr],trnaamino[curr],trnaanticodon[curr],intronseq[curr]))
    trnaloci.sort(key = lambda x: striplocus(x.name))
    for transname, currloci in itertools.groupby(trnaloci, lambda x: striplocus(x.name)):
        currlocus = list(currloci)
        locustranscripts[transname] = set(locus.name for locus in currlocus)
        yield tRNAtranscript( currlocus[0].seq, set(curr.score for curr in currlocus), currlocus[0].amino,currlocus[0].anticodon, set(currlocus), currlocus[0].intronseq, name = transname)
                                                                                                                                                                                    


def getmaturetrnas(**args):
    args = defaultdict(lambda: None, args)
    
    prokmode = args["prokmode"]
    margin = args["margin"]

    #>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-2 (tRNAscan-SE ID: chrVI.trna6) chrVI:204924-204996 (-) Ala (AGC) 73 bp Sc: 69.6
                                #nmt-tRNA-Gln-TTG-13-1
    trnanamere = re.compile(r"^\>\S+_((?:\w+\-)?\w+\-\w+\-[\w\?]+\-\d+\-\d+)\s+\((?:tRNAscan\-SE\s+ID:\s+)?([^\-]+)(\S+)\)")
    #trnanamere = re.compile(r"^\>\S+_(tRNA\-\w+\-\w+\-\d+\-\d+)\s+")#\((S+)\)")
    if "maturetrnafa" in args:
        maturetrnafa = open(args["maturetrnafa"], "w")
    else:
        maturetrnafa=sys.stdout
    gtrnatrans = None
    
    if args["namemap"]:
        gtrnafa = open(args["namemap"], "r")
        gtrnatrans = dict()
        for currline in gtrnafa:
            fields = currline.rstrip().split()


            if "tRNAscan-SE_id" == fields[0] or len(fields) < 2 or len(fields[1].split("-")) < 4:
                continue
            shortname = fields[0].split("-")[0]

            gtrnatrans[shortname] = fields[1] 

    elif args["gtrnafa"]:
        gtrnafa = open(args["gtrnafa"], "r")
        gtrnatrans = dict()
        for currline in gtrnafa:
            trnamatch = trnanamere.match(currline)
            #print >>sys.stderr, currline
            if trnamatch:

                #print >>sys.stderr, trnamatch.group(1)+":"+trnamatch.group(2)
                gtrnatrans[trnamatch.group(2)] = trnamatch.group(1)
            elif currline.startswith(">"):
                pass
                #print >>sys.stderr, currline
        if len(gtrnatrans.keys()) == 0:
            print >>sys.stderr, "Could not extract names from gtrnadb fasta file"
            print >>sys.stderr, "must have names in format '>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-10 (tRNAscan-SE ID: chrXIII.trna9)'"
            sys.exit()
            
            
    alltrnas = list()
    trnascantrnas = list()
    trnadbtrnas = list()
    trnacentraltrnas = list()
    for currfile in args["trnascan"]:
        if gtrnatrans:
            trnadbtrnas.extend(readtRNAdb(currfile, args["genome"], gtrnatrans))
        else:
            trnascantrnas.extend(readtRNAscan(currfile, args["genome"]))
            #print >>sys.stderr, len(trnascantrnas)
            
        
    '''
    for currfile in args["trnascan"]:
        trnascantrnas.extend(readtRNAscan(currfile, args["genome"]))
    trnacentraltrnas = list()
    for currfile in args["rnacentral"]:
        trnacentraltrnas.extend(readrnacentral(currfile,args.chromtranslate,mode = 'transcript'))
    '''    
    alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas + trnadbtrnas
    locustrnas = trnadbtrnas
    mitomode = args["mitomode"]
    trnabed = None
    if args["bedfile"]:
        trnabed = open(args["bedfile"], "w")
    
    trnatable = None
    if args["maturetrnatable"]:
        trnatable = open(args["maturetrnatable"], "w")
    
    
    
    
    def readmultistk(struct):
        currrecord = ""
        structs = list()
        for line in struct.split("\n"):
            currrecord += line+"\n"
            if line == "//":
                yield currrecord
                currrecord = ""
                
    
    
    anticodoncount = defaultdict(int)
    trnanames = dict()
    trnalist = list()
    for currtrans in alltrnas:
        if currtrans.name is None:
            name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
            currtrans.name = name
            anticodoncount[currtrans.anticodon] += 1
            
            
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    
    #cp /projects/lowelab/users/pchan/data/tRNA/stdModels/infernal-1.1/TRNAinf.cm ./
    

    trnacmfile = args["cmmodel"] 
    stkfile = args["trnaalignment"]
    if args["trnaalignment"]:
        devnull = open(os.devnull, 'w')
        seqfile = tempmultifasta(((currtrans.name, currtrans.getmatureseq(prokmode = prokmode)) for currtrans in alltrnas))

        cmcommand = ['cmalign', "-o", stkfile,"--nonbanded","--notrunc", "-g",trnacmfile,seqfile.name]
        #print >>sys.stderr, " ".join(cmcommand)
        #while True:
            #time.sleep(100)
        cmrun = subprocess.Popen(cmcommand, stdout = devnull)
        result = cmrun.wait()
        if result:
            print >>sys.stderr, "Failure to align tRNAs"
            sys.exit(1)
        #stkout = cmrun.communicate()[0]
        #trnaalign = readrnastk(stkout.split("\n"))[0]
        seqfile.close()
 
    if args["locibed"]:
        locibed = open(args["locibed"],"w")
  
    for currtrans in alltrnas:
        name = currtrans.name
        #trnanames[name] = currtrans
        trnalist.append(name)
        #print >>sys.stderr, name
        print >>maturetrnafa, ">"+name
        print >>maturetrnafa, str("N" * margin) +currtrans.getmatureseq(prokmode = prokmode)+str("N" * margin)
        if trnatable is not None:
            print >>trnatable, "\t".join([name,",".join(currlocus.name for currlocus in currtrans.loci),currtrans.amino,currtrans.anticodon])
        if trnabed is not None:
            transcriptrange = GenomeRange("genome", name, margin, margin + len(currtrans.getmatureseq(prokmode = prokmode)), strand = "+", name = name)
            print >>trnabed, transcriptrange.bedstring()
        if args["locibed"]:
            for currlocus in currtrans.loci:
                pass
                print >>locibed, currlocus.loc.bedstring()
        
def gettrnalocus(**args):
    #margin = 100
    margin = args["margin"]
    args = defaultdict(lambda: None, args)
    stkfile = args["stkfile"]
    genomefile = args["genomefile"]
    mitomode = args["mitomode"]
    trnalocifile = args["trnaloci"]
    trnalociseqs = args["trnalociseqs"]
    trnalocipos = args["trnalocipos"]
    trnacmfile = args["trnamodel"]                                              
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"        
    #if mitomode:
    #    trnacmfile = scriptdir+'TRNAinf.cm'
    #else:
    #    trnacmfile = scriptdir+'TRNAinf-euk.cm'
    
    trnaloci = list(readbed(trnalocifile, orgdb = "genome", seqfile=genomefile))
    lociseqs = getseqdict(trnaloci, faifiles = {"genome":genomefile+".fai"})
    #print lociseqs
    #lociseqfile = tempmultifasta(lociseqs)
    upstreamdist = dict()
    downstreamdist = dict()
    locisort = sorted(trnaloci, key = lambda x: int(x.start))
    for i in range(len(trnaloci) - 1):
        firgene = locisort[i]
        secgene = locisort[i+1]
        

        #if firgene.start > secgene.end and firgene.start - margin < secgene.end + margin:
        #    pass
        if firgene.end < secgene.start and firgene.end + margin > secgene.start - margin:

            
            distance = secgene.start - firgene.end
            downstreamdist[firgene.name] = margin - distance / 2
            upstreamdist[secgene.name] = margin - distance / 2
            #print >>sys.stderr, firgene.name
            #print >>sys.stderr, secgene.name
            #print >>sys.stderr, margin - distance / 2
        
    
    lociflanks = getseqdict(list(curr.addmargin(margin) for curr in trnaloci), faifiles = {"genome":genomefile+".fai"})
    #print >>sys.stderr, lociflanks.keys()
    for currgene in lociflanks.iterkeys():
        if currgene in downstreamdist:
            #print >>sys.stderr, lociflanks[currgene]
            lociflanks[currgene] = lociflanks[currgene][:-downstreamdist[currgene]] + downstreamdist[currgene] * "N"
            #print >>sys.stderr, lociflanks[currgene]
        if currgene in upstreamdist:
            lociflanks[currgene] = upstreamdist[currgene] * "N" + lociflanks[currgene][upstreamdist[currgene]:] 
        
    #sys.exit()
    lociseqfile = open(trnalociseqs, "w") 
    #print >>sys.stderr, trnalociseqs
    for name, seq in lociflanks.iteritems():
        
        print >>lociseqfile, ">"+name
        print >>lociseqfile, seq
    lociseqfile.close()
    
    locistarts = {curr: lociflanks[curr][0:margin].replace("T","U") for curr in lociflanks}
    lociends = {curr: lociflanks[curr][-margin:].replace("T","U") for curr in lociflanks}
    
    lociseqs = getseqdict(list(trnaloci), faifiles = {"genome":genomefile+".fai"})
    #creating loci stk file here
    #need to add flanking regions after alignment
    
    devnull = open(os.devnull, 'w')
    seqfile = tempmultifasta(lociseqs.iteritems())
    stktemp = tempfile.NamedTemporaryFile()
    cmcommand = ['cmalign', "-o",stktemp.name ,"--nonbanded","--notrunc",trnacmfile,seqfile.name] #"-g"
    #print >>sys.stderr, " ".join(cmcommand)
    cmrun = subprocess.Popen(cmcommand, stdout = devnull)
    result = cmrun.wait()
    if result:
        print >>sys.stderr, "Failure to align tRNAs"
        sys.exit(1)
    #sys.exit()

    devnull.close()    #lociseqs = getseqdict(trnaloci, faifiles = {"genome":genomefile+".fai"})
    locusstk = list(readrnastk(stktemp))[0]
    #locusstk = locusstk.addupstream(locistarts)
    #locusstk = locusstk.adddownstream(lociends)
    
    locistkfile = open(stkfile, "w")
    locusstk.printstk(output = locistkfile)
    locistkfile.close()
    trnalocifile = open(trnalocipos, "w")
    for currloci in trnaloci:
        transcriptrange = GenomeRange("genome", currloci.name ,margin, margin + currloci.length(), strand = "+", name = currloci.name)
        print >>trnalocifile, transcriptrange.bedstring()

    trnalocifile.close()

def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
        


parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--databasename',required=True,
                   help='database name to be used')
parser.add_argument('--genomefile',required=True,
                   help='Fasta file containing genome sequence')
parser.add_argument('--trnascanfile',required=True,
                   help='output from tRNAscan-SE run')
parser.add_argument('--gtrnafafile',
                   help='Fasta file of tRNA sequences from gtRNAdb')
parser.add_argument('--namemap',
                   help='gtrnadb file containing pairs of tRNA names')
parser.add_argument('--orgmode',
                   help='organism mode (euk/arch/bact/mito, default euk)')
parser.add_argument('--forcecca', action="store_true", default=False,
                   help='force the addition of a CCA tail')


def shellcall(shellcommand,failquit = False):
    retcode = subprocess.call(shellcommand, shell=True)
    if retcode > 0 and failquit:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
        print >>sys.stderr, "quiting program..."
        sys.exit(1)
    elif retcode > 0:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
    return retcode
        
        
args = parser.parse_args()
dbname = args.databasename
orgmode = args.orgmode
scanfile = os.path.expanduser(args.trnascanfile)
genomefile = os.path.expanduser(args.genomefile)
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
forcecca = args.forcecca


prokmode = False
if orgmode is None:
    orgmode = "euk"






if orgmode == "euk":
    maturemodel =  scriptdir+'trnamature-euk.cm'
    maturemargin = 20
    trnamodel =  scriptdir+'TRNAinf-euk.cm'
    locimargin = 50
elif orgmode == "arch":
    maturemodel =  scriptdir+'trnamature-arch.cm'
    maturemargin = 20
    if forcecca:
        trnamodel =  scriptdir+"TRNAinfnocca-arch.cm"
    else:
        trnamodel =  scriptdir+'TRNAinf-arch.cm'
    locimargin = 50
    prokmode = True
elif orgmode == "mito":
    maturemodel =  scriptdir+'TRNAMatureMitoinf.cm'
    maturemargin = 20
    if forcecca:
        trnamodel =  scriptdir+"trna1415G-nocca.cm"
    else:
        trnamodel =  scriptdir+'TRNAinf.cm'
    locimargin = 50
    prokmode = False
elif orgmode == "bact":
    maturemodel =  scriptdir+'bactmature.cm'
    maturemargin = 20
    if forcecca:
        trnamodel =  scriptdir+"bactnocca-030216.cm"
    else:
        trnamodel =  scriptdir+'bact-030216.cm'
    

    locimargin = 50
    prokmode = True
    
if forcecca:
    prokmode = False


gtrnafafile = None
namemap = None


if args.gtrnafafile:
    gtrnafafile = os.path.expanduser(args.gtrnafafile)
if args.namemap:

    namemap = os.path.expanduser(args.namemap)


    
#$1 is database name
#trnascanfile is trnascan file
#genomefile is fasta file of genome

#test command line programs

runtime = time.time()
loctime = time.localtime(runtime)
    
get_location("samtools")
get_location("bowtie2-build")

if not os.path.isfile(genomefile+".fai"):
    shellcall("samtools faidx "+genomefile)

getmaturetrnas(trnascan=[scanfile], genome=genomefile,gtrnafa=gtrnafafile,bedfile=dbname+"-maturetRNAs.bed",maturetrnatable=dbname+"-trnatable.txt",trnaalignment=dbname+"-trnaalign.stk",locibed=dbname+"-trnaloci.bed",maturetrnafa=dbname+"-maturetRNAs.fa", namemap = namemap,cmmodel = maturemodel, prokmode = prokmode, margin = maturemargin)
gettrnalocus(genomefile=genomefile,stkfile=dbname+"-trnaloci.stk",trnaloci=dbname+"-trnaloci.bed",maturetrnafa=dbname+"-maturetRNAs.fa", trnalocipos = dbname+"-trnaloci.bed",trnalociseqs = dbname+"-tRNAloci.fa", margin = locimargin, trnamodel = trnamodel)

#getmaturetrnas(trnascan=[scanfile], genome=genomefile,gtrnafa=gtrnafafile,namemap=namemapfile, bedfile=dbdirectory+dbname+"-maturetRNAs.bed",maturetrnatable=dbdirectory+dbname+"-trnatable.txt",trnaalignment=transcriptstk,locibed=dbdirectory+dbname+"-trnaloci.bed",maturetrnafa=dbdirectory+dbname+"-maturetRNAs.fa",cmmodel = maturemodel, prokmode = prokmode)
#aligntrnalocus(genomefile=genomefile,stkfile=locusstk,trnaloci=dbdirectory+dbname+"-trnaloci.bed", cmmodel = trnamodel)


gitversion, githash = getgithash(scriptdir)

intronseqs = allintronseqs

#print >>sys.stderr, allintronseqs

conv = convertalignment(dbname+"-trnaalign.stk",intronseqs, locustranscripts)
#conv.printstk(output = open(dbname+"-convertedloci.stk", "w"))

revconv = splicealignment(dbname+"-trnaloci.stk",intronsplices,locustranscripts, prokmode = prokmode)
revconv.printstk(output = open(dbname+"-trnaconvert.stk", "w"))


combinedseqs = open(dbname+"-tRNAgenome.fa",'w')
catcall = subprocess.Popen([find_executable("cat"),dbname+"-maturetRNAs.fa",dbname+"-tRNAloci.fa"],stdout = combinedseqs )
catcall.wait()
combinedseqs.close()

shellcall("samtools faidx "+dbname+"-tRNAgenome.fa")
dbdirectory = os.path.dirname(dbname) + "/"
if dbdirectory == "/":
    dbdirectory = ""
print >>sys.stderr, dbdirectory
#shellcall("cat "+dbname+"-maturetRNAs.fa "+dbname+"-tRNAlocipositions.fa >"+dbname+"-tRNAgenome.fa", failquit = True)
bowtie2buildcall = subprocess.Popen([find_executable("bowtie2-build"),dbname+"-tRNAgenome.fa",dbname+"-tRNAgenome"])
bowtie2buildcall.wait()

#bowtie2buildcall = subprocess.Popen([find_executable("bowtie2-build"),dbname+"-maturetRNAs.fa",dbname+"-maturetRNAs"])
#bowtie2buildcall.wait()
#print >>sys.stderr, " ".join([find_executable("bowtie2-build"),dbname+"-tRNAgenome.fa",dbname+"-tRNAgenome"])

#shellcall("bowtie2-build "+dbname+"-tRNAgenome.fa "+dbname+"-tRNAgenome", failquit = True)
dbinfo = open(dbname+ "-dbinfo.txt","w")
print >>dbinfo, "time\t"+str(runtime)+" ("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "creation\t"+" ".join(sys.argv)
print >>dbinfo, "genomefile\t"+str(genomefile)
print >>dbinfo, "trnascanfile\t"+str(scanfile)
print >>dbinfo, "orgmode\t"+str(orgmode)
print >>dbinfo, "forcecca\t"+str(forcecca)
print >>dbinfo, "trna_margin\t"+str(maturemargin)
print >>dbinfo, "loci_margin\t"+str(locimargin)
print >>dbinfo, "git_version\t"+str(gitversion)
dbinfo.close()
