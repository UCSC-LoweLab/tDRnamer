#!/usr/bin/env python

import pysam
import sys
import tempfile
import re
import gzip
from collections import defaultdict

gapchars = set("-.~")


def is_number(num):
    try:
        int(num)
        return True
    except ValueError:
        return False


def trnatranscriptkey(trnaname):
    fields = trnaname.split("-")

        
    return fields[1],fields[2],int(fields[3])

def trnalocuskey(trnaname):
    fields = trnaname.split("-")

        
    return fields[1],fields[2],int(fields[3]),int(fields[4])

def is_locusname(trnaname):
     fields = trnaname.split("-")
     if len(fields) != 5:
         return False
         
     if not is_number(fields[3]) or not is_number(fields[4]):
         return False
     return True
def is_transcriptname(trnaname):
     fields = trnaname.split("-")
     if len(fields) != 4:
         return False
     if not is_number(fields[3]):
         return False
     return True
def trnasort(trnanames):
    trnatranscripts = sorted((curr for curr in trnanames if is_transcriptname(curr)), key = trnatranscriptkey)
    trnaloci = sorted((curr for curr in trnanames if is_locusname(curr)), key = trnalocuskey)
    other = list(curr for curr in trnanames if not is_transcriptname(curr) and not is_locusname(curr))
    return trnatranscripts+trnaloci+other

def readmultifasta(fafile):
    #print chrom+":"+ chromstart+"-"+ chromend
    if fafile == "stdin":
        fafile = sys.stdin
    else:
        fafile = open(fafile, "r")
    currloc = 0
    currseq = None
    sequence = ""
    reheader = re.compile(r"\>([^\s\,]+)")
    for line in fafile:
        line = line.rstrip("\n")
        currheader = reheader.match(line)
        if currheader and sequence != "" and currseq is not None:
            yield currseq, sequence
            currseq = currheader.groups(1)[0]
            sequence = ""
        elif currheader:
            currseq = currheader.groups(1)[0]
        else:
            sequence += line
    if currseq is not None:
        yield currseq, sequence
            
def fastadict(fafile):
    seqdict = dict()
    for name, seq in readmultifasta(fafile):
        seqdict[name] = seq
    return seqdict
        
def tempmultifasta(allseqs):
    fafile = tempfile.NamedTemporaryFile(suffix=".fa")


    for seqname, seq in allseqs:
        fafile.write(">"+seqname+"\n")
        fafile.write(seq+"\n")
    fafile.flush()
    return fafile
def invertstrand(strand):
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"

class alignment:
    def __init__(self, sequences):
        if len(sequences.keys()) < 1:
            raise EmptyAlignException() 
        self.aligns = sequences

        if max(len(curr) for curr in self.aligns.itervalues()) != min(len(curr) for curr in self.aligns.itervalues()):
            print >>sys.stderr, "Non-matching sequence lengths in multiple alignment"
        #self.alignlength = len(self.aligns[self.aligns.keys()[0]])
        self.alignlength = max(len(curr) for curr in self.aligns.itervalues())
    def getseqlength(self, seqname):
        return len(self.aligns[seqname].replace("-",""))
    def numseqs(self):
        return len(self.aligns.keys())
    def getseqnames(self):
        return self.aligns.keys()
    def getalignseq(self ,seqname):
        return self.aligns[seqname]
    def toseqdict(self):
        seqs = dict()
        for currname, currseq in self.aligns.iteritems():
            seqs[currname] = string.translate(currseq, None, "-.~")
        return seqs
    def getsubset(self, subset):
        newaligns = dict()
        for curr in subset:
            newaligns[curr] = ""
        for i in range(self.alignlength):
            if any(self.aligns[curr][i] not in gapchars for curr in subset):
                for curr in subset:
                    newaligns[curr] += self.aligns[curr][i]
        return alignment(newaligns)
    def prottocodonalignment(self, nucseqs):
        #I don't do any checking to make sure that the nucleotide sequence is equivalent here
        newalign = defaultdict(str)
        currpoint = defaultdict(int)
        for i in range(self.alignlength):
            for currgene in nucseqs.iterkeys():
                if self.aligns[currgene][i] != "-":
                    newalign[currgene] += nucseqs[currgene][currpoint[currgene]:currpoint[currgene] + 3]
                    currpoint[currgene] += 3
                else:
                    newalign[currgene] += "---"
         
        #print "|".join(str(curr) for curr in currpoint.values())
        #print "|".join(str(len(curr)) for curr in nucseqs.values())
        #print "|".join(curr for curr in newalign.values())
        
        #print str(self.alignlength) pyroIsla1  
        return alignment(newalign, aligned = True)
    def tempmultifasta(self):
        return tempmultifasta(self.aligns)
    def fastaformat(self, prefix = ""):
        output = ""
        for seqname, seq in self.aligns.iteritems():
            output += ">"+prefix+seqname+"\n"+seq+"\n"
        return output
    def getsegment(self,start, end):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            newalign[name] = seq[start:end]
        return alignment(newalign)
    def getsubsets(self,windowsize, stepsize):
        for curr in range(0, self.alignlength - windowsize, stepsize):
            start = curr
            end = curr + windowsize
            yield self.getsegment(start,end)
    def removeemptyseqs(self):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            if len(string.translate(seq, None, "-.")) != 0:
                newalign[name] = seq
        return alignment(newalign)
    #http://bioinformatics.oxfordjournals.org/content/25/5/668.full

    def phylipformat(self):
        #print ",".join(str(curr)for curr in self.aligns.keys())
        output = str(len(self.aligns.keys()))+ " "+str(self.alignlength)+"\n"
        for currgene in self.aligns.iterkeys():
            #sys.stderr.write("**"+str(currgene).ljust( 11)+"**\n")
            output += str(currgene).ljust( 14)
            output +=  self.aligns[currgene]+"\n"
        #print output
        return output


    def nexusformat(self):
        output = "#NEXUS\nBegin data;\nDimensions ntax="+str(len(self.aligns.keys()))+" nchar="+str(self.alignlength)+";\n"
        output += "Format datatype=dna symbols=\""+"ATCG"+"\" missing=? gap=-;\n"
        output += "Matrix\n"
        for currgene in self.aligns.iterkeys():
            output += currgene.ljust( 11)
            output +=  self.aligns[currgene]+"\n"
        output+=";\nEnd;"    
    def getrealcoordinate(self, alignname, coord):
        currreal = 0
        #print self.aligns[alignname]
        #print coord
        for i in range(0, coord):
            if self.aligns[alignname][i] not in  set("-.~"):
                currreal += 1
        return min([currreal, len(list(curr for curr in self.aligns[alignname] if curr not in set("-.~"))) - 1])
    def getseqrange(self, alignname, start, end):
        return self.getsegment(self.getaligncoordinate(alignname, start),self.getaligncoordinate(alignname, end))

            
    def printstk(self, name = None, annotation = None):
        print "# STOCKHOLM 1.0"
        if name is not None:
            print "#=GF ID "+name
        padlength = max(len(curr) for curr in self.aligns.iterkeys()) + 4
        for currname, currseq in self.aligns.iteritems():
            tail = ""
            if annotation is not None:
                tail = "    ["+annotation[currname]+']'
            print string.ljust(currname,padlength ) +currseq+tail
        print "//"
    def printhtml(self, name = None):
        print "<CODE>"
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"<BR/>"
        print "</CODE>"
    def clustalwformat(self):
        output = "CLUSTAL W 2.1 multiple sequence alignment\n\n"
        conservestring = ""
        for i in range(0, self.alignlength):
            conservestring += ":"                     
        for currpos in range(0, self.alignlength, 60):
            for seqname, seq in self.aligns.iteritems():
                output += seqname+"\t"+seq[currpos:min([currpos+60,self.alignlength])]+"\n"
            output += "\t"+conservestring[currpos:min([currpos+60,self.alignlength])] + "\n\n"
            
        return output

def getbestmappings(allmaps, bamfile):
    currscore = None
    newset = set()
    for currmap in allmaps:
        tagdict = dict()
        for curr in currmap.tags:
            tagdict[curr[0]] = curr[1]
        if currmap.tid is -1:
            continue
        #print >>sys.stderr, bamfile.getrname(currmap.tid)
        chromname = bamfile.getrname(currmap.tid)
        #sys.exit()
        readlength = len(currmap.seq)
        #if this is the best mapping of this read, discard any previous mappings            
        #print >>sys.stderr, tagdict["AS"]
        #if the current score is worse than the new one
        if currscore is None or currscore < tagdict["AS"]:
            #print >>sys.stderr, str(currscore) +"<"+str( tagdict["AS"])
            newset = set()
            newset.add(currmap)
            currscore = tagdict["AS"]
        #if this mappings is the same as the previous best mapping, add it to the list
        elif currscore == tagdict["AS"]:
            newset.add(currmap)
            
        else:
            pass
        #print >>sys.stderr, currscore

    return newset

class RnaAlignment(alignment):
    def __init__(self,alignseqs, structure, consensus = None,energies = None):
        self.aligns = alignseqs
        self.currstruct = structure
        self.energies = energies 
        self.consensus = consensus
        self.alignlength = max(len(curr) for curr in alignseqs.values())
        if self.alignlength != len(structure):
            for curr in alignseqs.iterkeys():
                print >>sys.stderr, alignseqs[curr]
            print >>sys.stderr, structure
            
            raise AlignError('Alignment structure length not equal to sequence lengths')
            #print >>sys.stderr, "Cannot create alignment, " 
            #sys.exit(1)
    def addupstream(self, seqs, struct = None, consensus = None):
        newseqs = dict()
        #print >>sys.stderr, self.aligns.keys()
        for currname in self.aligns.iterkeys():
            newseqs[currname] = seqs[currname] + self.aligns[currname]
        if struct is None:
            newstruct = (max(len(curr) for curr in seqs.itervalues()) * ":") + self.currstruct
        else:
            newstruct = struct + self.currstruct
        if self.consensus is None:
            pass
        elif consensus is not None:
            newconsensus = consensus + self.consensus
            
        else:
            newconsensus = (max(len(curr) for curr in seqs.itervalues()) * ".") + self.consensus
            

        return RnaAlignment(newseqs, newstruct, consensus = newconsensus)
    def adddownstream(self, seqs, struct = None, consensus = None):
        newseqs = dict()
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  self.aligns[curr]+  seqs[curr]
        if struct is None:
            newstruct = self.currstruct + (max(len(curr) for curr in seqs.itervalues()) * ":")
        else:
            newstruct = self.currstruct + struct
        if self.consensus is None:
            pass
        elif consensus is not None:
            newconsensus = self.consensus + consensus
            
        else:
            newconsensus = self.consensus + (max(len(curr) for curr in seqs.itervalues()) * ".")
            
            
        return RnaAlignment(newseqs, newstruct, consensus = newconsensus)
    def addmargin(self, length):
        newseqs = dict()
        newconsensus = None
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  length*"N" + self.aligns[curr]+ length*"N"
        newstruct = length * "-" + self.currstruct + length * "-"
        if self.consensus:
            newconsensus = length * "-" + self.consensus + length * "-"

        return RnaAlignment(newseqs, newstruct, newconsensus)
        
    def splice(self, splices):
        newaligns = dict()
        currpoint = defaultdict(int)
        for curr in self.aligns.keys():
            newaligns[curr] = ''
        for i in range(self.alignlength):
            for curr in self.aligns.keys():
                
                if curr in splices and splices[curr][0] <= currpoint[curr] < splices[curr][1]:
                    newaligns[curr] += '-'
                else:
                    newaligns[curr] += self.aligns[curr][i]
                if self.aligns[curr][i] not in gapchars:
                    currpoint[curr] += 1
        #print >>sys.stderr, "||||||||||||||||||||||||||||||||"            
        for currname in splices.keys():
            newseq = "".join(curr for curr in newaligns[currname] if curr not in gapchars)
            oldseq = "".join(curr for curr in self.aligns[currname] if curr not in gapchars)
            oldranged = (oldseq[:splices[currname][0]] + oldseq[splices[currname][1]:])
            if newseq.upper() != oldranged.upper():
                #print >>sys.stderr, currname+"|"
                #print >>sys.stderr, newseq
                #print >>sys.stderr, oldranged
                pass
                    
            if currname in set(["tRNA-Leu-CAA-1-1","tRNA-Leu-CAA-1-2"]):
                #print >>sys.stderr, [splices[currname][0],splices[currname][1]]
                #print >>sys.stderr, oldseq
                #print >>sys.stderr, newseq
                pass
        return RnaAlignment(newaligns,self.currstruct , consensus = self.consensus) 
        
    def multisplice(self, splices):
        splpoint = 0
        newaligns = dict()
        currpoint = defaultdict(int)
        for curr in self.aligns.keys():
            newaligns[curr] = ''
        for i in range(self.alignlength):
            for curr in self.aligns.keys():
                
                if curr in splices and splices[curr][splpoint][0] <= currpoint[curr] < splices[curr][splpoint][1]:
                    newaligns[curr] += '-'
                else:
                    newaligns[curr] += self.aligns[curr][i]
                if self.aligns[curr][i] not in gapchars:
                    currpoint[curr] += 1
        #print >>sys.stderr, "||||||||||||||||||||||||||||||||"            
        for currname in splices.keys():
            newseq = "".join(curr for curr in newaligns[currname] if curr not in gapchars)
            oldseq = "".join(curr for curr in self.aligns[currname] if curr not in gapchars)
            oldranged = (oldseq[:splices[currname][splpoint][0]] + oldseq[splices[currname][splpoint][1]:])
            if newseq.upper() != oldranged.upper():
                #print >>sys.stderr, currname+"|"
                #print >>sys.stderr, newseq
                #print >>sys.stderr, oldranged
                pass
                    
            if currname in set(["tRNA-Leu-CAA-1-1","tRNA-Leu-CAA-1-2"]):
                #print >>sys.stderr, [splices[currname][0],splices[currname][1]]
                #print >>sys.stderr, oldseq
                #print >>sys.stderr, newseq
                pass
        return RnaAlignment(newaligns,self.currstruct , consensus = self.consensus) 
        
    def getsubset(self, subset):
        newaligns = dict()
        newstruct = ""
        if self.consensus is not None:
            newconsensus = ''
        for curr in subset:
            newaligns[curr] = ""
            if curr not in self.aligns:
                print >>sys.stderr, "No such sequence "+curr
                sys.exit(1)
        for i in range(self.alignlength):
            if any(self.aligns[curr][i] not in gapchars for curr in subset) or self.currstruct[i] not in gapchars:
                for curr in subset:
                    newaligns[curr] += self.aligns[curr][i]
                newstruct += self.currstruct[i]
                if self.consensus is not None:
                    newconsensus += self.consensus[i]
                    
        return RnaAlignment(newaligns,newstruct, newconsensus)
    def trimcolumns(self, columns):
        newaligns = dict()
        newstruct = ""
        if self.consensus is not None:
            newconsensus = ''
        else:
            newconsensus = None
        for curr in self.aligns.keys():
            newaligns[curr] = ''
        for i in range(self.alignlength):
            if i not in columns:
                for curr in self.aligns.keys():
                    newaligns[curr] += self.aligns[curr][i]
                newstruct += self.currstruct[i]
                if self.consensus is not None:
                    newconsensus += self.consensus[i]
                    
        return RnaAlignment(newaligns,newstruct , consensus = newconsensus) 
    def convertnames(self, newnames):
        newaligns = dict()
        newstruct = self.currstruct
        newconsensus = self.consensus
        for currkey in newnames.keys():
            for newname in newnames[currkey]:
                newaligns[newname] = self.aligns[currkey]
                
        return RnaAlignment(newaligns,newstruct, consensus = newconsensus) 

    def insertcolumns(self, column, insertalign):
        newaligns = dict()
        newstruct = ""
        if self.consensus is not None:
            newconsensus = ''
        else:
            newconsensus = None
        for curr in self.aligns.keys():
            newaligns[curr] = ''
        for i in range(self.alignlength):
            if i == column:
                for curr in newaligns.keys():
                    #print >>sys.stderr, curr
                    #print >>sys.stderr, insertalign.aligns[curr]
                    newaligns[curr] += insertalign.aligns[curr]
                newstruct += insertalign.currstruct
                if self.consensus is not None:
                    newconsensus += insertalign.consensus
                    
            for curr in self.aligns.keys():
                newaligns[curr] += self.aligns[curr][i]
            newstruct += self.currstruct[i]
            if self.consensus is not None:
                newconsensus += self.consensus[i]
        return RnaAlignment(newaligns,newstruct, consensus = newconsensus) 
        
        
    def viennaformat(self):
        output = ""
        for currseq in self.aligns.iterkeys():
            output += ">"+currseq+"\n"
            output += self.aligns[currseq]+"\n"
            output += self.currstruct+"\n"
        return output
    def viennatempfile(self):
        viennafile = tempfile.NamedTemporaryFile()
        viennafile.write(self.viennaformat())
        viennafile.flush()
        return viennafile
    def printstk(self, name = None, ordering = None, annotation = None, output = sys.stdout):
        maxnamelength = max(len(name) for name in self.aligns.keys())
        maxnamelength = max([maxnamelength,len("#=GC SS_cons")])
        print >>output, "# STOCKHOLM 1.0"
        if name is not None:
            print  >>output, "#=GF ID "+name
        if ordering is None:
            outlist = self.aligns.keys() 
        else:
            outlist = ordering
        for currname in outlist:
            tail = ''
            if annotation is not None:
                tail = "    ["+annotation[currname]+']'
            currseq = self.aligns[currname]
            print  >>output, currname.ljust(maxnamelength+4)+currseq+tail
        structline = ""
        structpos = 0
        secpos = 0
        print  >>output, "#=GC SS_cons".ljust(maxnamelength+4)+self.currstruct
        if self.consensus is not None:
            print  >>output, "#=GC RF".ljust(maxnamelength+4)+self.consensus
        print >>output, "//"
    def pruneNbases(self):
        newaligns = dict()
        newstruct = ""
        newconsensus = None
        if self.consensus is not None:
            newconsensus = ''
        for curr in self.aligns.iterkeys():
            newaligns[curr] = ""
            if curr not in self.aligns:
                print >>sys.stderr, "No such sequence "+curr
                sys.exit(1)
        for i in range(self.alignlength):
            for curr in self.aligns.iterkeys():
                if self.aligns[curr][i] == "N":
                    newaligns[curr] += "-"
                else:
                    newaligns[curr] += self.aligns[curr][i] 
            newstruct += self.currstruct[i]
            if self.consensus is not None:
                newconsensus += self.consensus[i]
        return RnaAlignment(newaligns,newstruct , consensus = newconsensus)
        
    def printhtml(self, name = None):
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"</BR>"
        structline = ""
        structpos = 0
        secpos = 0
        print "#=GC SS_cons\t"+self.currstruct+"</BR>"
        print "//"+"</BR>"
sprinzeleuk = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e5','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
sprinzelarch = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])

sprinzeladditional = set(['17a','20a','20b'])

sprinzelbact = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e27','e26','e25','e24','e23','e22','e21',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])

sprinzelmito = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e11','e12','e13','e14','e15','e16','e17','e1','e2','e3','e4','e27','e26','e25','e24','e23','e22','e21',47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
def gettrnanums(trnaalign, margin = 20, orgtype = "euk"):
    
    trnanum = list()
    currcount = 0
    enum = 11
    gapnum = 1
    intronnum = 1
    sprinzelpositions = sprinzeleuk
    if orgtype == "bact":
        sprinzelpositions = sprinzelbact
    elif orgtype == "arch":
        sprinzelpositions = sprinzelarch
    elif orgtype == "mito":
        sprinzelpositions = sprinzelmito
        
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    insertlength = 0    
    lastpos = '0'
    taillength = 0
    for i, struct in enumerate(trnaalign.consensus):
        
        if currcount >= len(sprinzelpositions):
            
            trnanum.append(str(lastpos)+"i"+str(insertlength+1))
            #trnanum.append("tail"+str(insertlength+1))
            gapnum += 1
            currcount += 1
            insertlength += 1
            
        elif struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=':
                currcount = 1
                insertlength = 0
            if sprinzelpositions[currcount] == 'e':
                trnanum.append(sprinzelpositions[currcount])
                insertlength = 0
                lastpos = sprinzelpositions[currcount]
                #if enum == 17:
                #    enum = 1
                #elif enum == 5:
                #    enum = 27
                #elif enum > 20:
                #    enum -= 1
                #else:
                #    enum += 1
                currcount += 1
            elif sprinzelpositions[currcount] == '-':
                insertlength += 1
                #trnanum.append('gap'+str(gapnum))
                trnanum.append(str(lastpos)+"i"+str(insertlength+1))
                gapnum += 1
                currcount += 1
            else:
                trnanum.append(str(sprinzelpositions[currcount]))
                lastpos = sprinzelpositions[currcount]
                currcount += 1
                insertlength = 0
        else:
            if sprinzelpositions[currcount] == 38:
                trnanum.append('intron'+str(intronnum))
                intronnum += 1
            else:
                
                #trnanum.append('gap'+str(gapnum))
                trnanum.append(str(lastpos)+"i"+str(insertlength+1))
                
                insertlength += 1
                gapnum += 1
        #print >>sys.stderr, i + margin + 1
        #print >>sys.stderr, trnanum
        #print >>sys.stderr, i
        #print >>sys.stderr, "".join("-" if curr != "76" else "*" for curr in trnanums)
        #print >>sys.stderr, "".join("*" if curr == i else "-" for curr in range(len(trnaalign.consensus)))
        #print >>sys.stderr, "".join(str(curr) for curr in trnaalign.consensus)
        #print >>sys.stderr, "".join(str(curr)[-1] if is_number(curr) else "-" for curr in trnanum)
        #print >>sys.stderr, "**"
        if i + margin + 1 != len(trnanum):
            #print >>sys.stderr, trnaalign.consensus[:i]
            #print >>sys.stderr, trnanum[-1]
            #print >>sys.stderr, i
            #print >>sys.stderr, len(trnanum)
            #print >>sys.stderr, ",".join(trnanum)
            #print >>sys.stderr, "**||"
            #sys.exit(1)
            pass
    for i in range(margin):
        trnanum.append('tail'+str(i+1))
    if len(trnanum) != len(trnaalign.consensus):
        #print >>sys.stderr, "Could not calculate tRNA numbering"
        #sys.exit(1)
        pass
        
    
    return trnanum
def replacepadding(trnanums, margin = 20):
    newnum = list()
    for i in range(margin):
        newnum.append('head'+str(margin - i))
    for i in range(margin, len(trnanums) - margin):
        newnum.append(trnanums[i])
    for i in range(margin):
        newnum.append('tail'+str(i +1))
    return newnum
    
def convertmaturealign(rnaalign):
    newseqs = dict()
    newalign = rnaalign.struct
    for name, seq in rnaalign.aligns.iteritems():
        newseqs[name]
    
def readrnastk(stk):
    seqs = defaultdict(str)
    struct = ""
    consensus = ""
    energyscore = None
    for line in stk:
        line = line.rstrip()
        if line.startswith("//"):
            if consensus == "":
                consensus = None
            yield RnaAlignment(seqs, struct, consensus = consensus)
            seqs = defaultdict(str)
            struct = ""
            consensus = ""
            energyscore = None
        elif not line.startswith("#") and len(line.split()) > 1:
            currname = line.split()[0]
            currseq = line.split()[1]
            seqs[currname] += currseq
        elif line.startswith("#=GC SS_cons"):
            struct += line.split()[2]
        elif line.startswith("#=GC RF"):
            consensus += line.split()[2]
            

class transcriptfile:
    def __init__(self, trnafilename):
        trnafile = open(trnafilename)
        locustranscript = dict()
        trnatranscripts = list()
        loci = list()
        amino = dict()
        anticodon = dict()
        for i, line in enumerate(trnafile):
            fields = line.split()
            if len(fields) < 2:
                continue
            trnatranscripts.append(fields[0])
            amino[fields[0]] = fields[2]
            anticodon[fields[0]] = fields[3]
            for currlocus in fields[1].split(','):
                locustranscript[currlocus] = fields[0]
                loci.append(currlocus)

        
        self.locustranscript = locustranscript
        self.transcripts = trnatranscripts
        self.amino = amino
        self.anticodon = anticodon
        self.loci = loci
    def gettranscripts(self):
        return set(self.transcripts)
    def getlocustranscript(self, locus):
        return  self.locustranscript[locus]
    def getloci(self):
        return  self.loci
    def getamino(self, trna):
        return  self.amino[trna]
    def getanticodon(self, trna):
        return  self.anticodon[trna]
        
    def allaminos(self):
        return  set(self.amino.values())
    def allanticodons(self):
        return  set(self.anticodon.values())
        


class samplefile:
    def __init__(self, samplefilename):
        try:
            samplefile = open(samplefilename)
            samplelist = list()
            samplefiles = dict()
            replicatename = dict()
            
            replicatelist = list()
            for i, line in enumerate(samplefile):
                fields = line.split()
                if len(fields) < 2:
                    continue
                samplefiles[fields[0]] = fields[2]
                replicatename[fields[0]] = fields[1]
                
                samplelist.append(fields[0])
                if fields[1] not in set(replicatelist):
                    replicatelist.append(fields[1])
            
            #bamlist = list(curr + "_sort.bam" for curr in samplefiles.iterkeys())
            samplenames = list(curr  for curr in samplefiles.iterkeys())
            self.samplelist = samplelist
            self.samplefiles = samplefiles
            self.replicatename = replicatename
            self.replicatelist = replicatelist
            #self.bamlist = list(curr+ "_sort.bam" for curr in samplelist)
        except IOError:
            print >>sys.stderr, "Cannot read sample file "+samplefilename
            print >>sys.stderr, "exiting..."
            sys.exit(1)
    def getsamples(self):
        return self.samplelist
    def getbamlist(self):
        return list(curr+ ".bam" for curr in samplelist)
    def getbam(self, sample):
        return sample+ ".bam"
    def getfastq(self, sample):
        return self.samplefiles[sample]
    def getreplicatename(self, sample):
        return self.replicatename[sample]
    def allreplicates(self):
        return self.replicatelist
    def getrepsamples(self, replicate):
        return set(currsample for currsample in self.replicatename.iterkeys() if self.replicatename[currsample] == replicate)
        
def getsizefactors( sizefactorfile):
    sizefactorfile = open(sizefactorfile)
    sizefactors = dict()
    bamheaders =list()
    sizes = list()
    for i, line in enumerate(sizefactorfile):
        if i == 0:
            bamheaders = list(curr.strip("\"\n") for curr in line.split())
        elif i == 1:
            sizes = list(float(curr.strip("\"\n")) for curr in line.split())
    for i in range(0, len(bamheaders)):
        sizefactors[bamheaders[i]] = sizes[i]
        #print >>sys.stderr, bamheaders[i]+":"+ str(sizes[i])
    return sizefactors

#special class that uses the read indentifier for hashing in sets

class GenomeRange:
    __slots__ = "dbname", "chrom", "strand","name", "fastafile"
    def __eq__(self, other):
        return self.strand == other.strand and self.chrom == other.chrom and self.start == other.start and self.end == other.end
    def __hash__(self):
        return  self.start + self.end + hash(self.chrom) + hash(self.strand)
    def __init__(self, dbname, chrom, start, end, strand = None,name = None, orderstrand = False, data = None, fastafile = None):
        self.dbname =dbname
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.fastafile = fastafile
        if orderstrand and self.start > self.end:
            temp = self.start
            self.start = self.end
            self.end = temp
            self.strand = invertstrand(strand)
        self.data = data
        self.name = name        
    def coverage(self, other):
        if self.strand == other.strand and self.chrom == other.chrom:
            start = max([self.start,other.start])
            end = min([self.end,other.end])
            
            if end - start < 0:
                return 0
            else:
                return end - start
        else:
            return 0
            
    def bedstring(self, name = None,score = 1000):
        if self.strand == None:
            self.strand = "+"
        if name is None and self.name is None:
            name = "FEAT"
        elif name is None:
            name = self.name
        return "\t".join([self.chrom,str(self.start),str(self.end),name,str(score), self.strand])
    def length(self):
        return self.end - self.start
    def addmargin(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist,self.end + dist,self.strand, name = newname, fastafile = self.fastafile)
    def getupstream(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist,self.start,self.strand, name = newname, fastafile = self.fastafile)
    def getdownstream(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand != '-':
            return GenomeRange(self.dbname, self.chrom, self.end,self.end + dist,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start - dist,self.start,self.strand, name = newname, fastafile = self.fastafile)
    def getfirst(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - dist,self.end ,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start,self.start + dist,self.strand, name = newname, fastafile = self.fastafile)
    def getlast(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.start,self.start + dist,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.end - dist,self.end ,self.strand, name = newname, fastafile = self.fastafile)
    def getbase(self, basenum, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - basenum - 1 ,self.end -basenum,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start +basenum ,self.start + basenum + 1,self.strand, name = newname, fastafile = self.fastafile)
    def antisense(self):
        if self.strand == "+":
            newstrand = "-"
        else:
            newstrand = "+"
        return GenomeRange(self.dbname, self.chrom, self.start,self.end,newstrand, name = self.name, fastafile = self.fastafile)
    def bamseq(self):
        if self.strand == "+":
            return self.data["seq"]
        else:
            return revcom(self.data["seq"])

class GenomeRead(GenomeRange):
    def __eq__(self, other):
        return self.name == other.name
    def __hash__(self):
        return  hash(self.name)
    def __init__(*args, **nargs):
        GenomeRange.__init__(*args, **nargs)
'''
Still need to add trailer fragment code, both here and elsewhere
'''
def getfragtype(currfeat, currread, maxoffset = 10):
    if currread.start < currfeat.start + maxoffset and currread.end > currfeat.end - maxoffset:
        return "Whole"
    elif currread.start < currfeat.start + maxoffset:
        if currfeat.strand == "+":
            return "Fiveprime"
        else:
            return "Threeprime"
    elif currread.end > currfeat.end - maxoffset:
        if currfeat.strand == "+":
            return "Threeprime"
        else:
            return "Fiveprime"
            
smallrnatypes = set([])            
def readfeatures(filename, orgdb="genome", seqfile= None, removepseudo = False):
    if filename.endswith(".bed") or filename.endswith(".bed.gz"):
        return readbed(filename, orgdb, seqfile)
    elif filename.endswith(".gtf") or filename.endswith(".gtf.gz") or filename.endswith(".gff") or filename.endswith(".gff.gz"):
        #print >>sys.stderr, removepseudo
        return (curr for curr in readgtf(filename, orgdb, seqfile, filterpsuedo = removepseudo, filtertypes =set(['retained_intron','antisense','lincRNA']) ))
    else:
        print >>sys.stderr, filename+" not valid feature file"
        sys.exit()


def readgtf(filename, orgdb="genome", seqfile= None, filterpsuedo = False, filtertypes = set(['retained_intron','antisense','lincRNA'])):
    bedfile = None
    #print >>sys.stderr, "****"
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split("\t")
        if len(fields) > 2:
            biotype = None
            featname = None
            genename = None
            #print >>sys.stderr, len(fields)
            genesource = fields[1]  
            #retained introns are often other things as well, so I skip em
            if fields[2] != "transcript" or genesource in filtertypes:
                continue

              
            for currattr in fields[8].rstrip(";").split(";"):
                #print >>sys.stderr,  currattr
                currname = currattr.strip().split()[0]
                currvalue = currattr.strip().split()[1]
                if currname == "transcript_name":
                    featname = currvalue.strip('"')
                elif currname == "name" or currname == "gene_id" and featname is None:
                    featname = currvalue.strip('"')

                elif currname == "gene_biotype":
                    biotype = currvalue.strip('"')
                elif currname == "gene_name":
                    genename = currvalue.strip('"')
                #print >>sys.stderr, "***||"
            if genename is None:
                #print >>sys.stderr, "No name for gtf entry "+featname
                genename = featname
                pass
            if filterpsuedo and biotype == "pseudogene":
                #print >>sys.stderr, "*******"
                continue
                
            if genesource == 'ensembl':
                #print >>sys.stderr, biotype
                genesource = biotype
            if not (fields[6] == "+" or fields[6] == "-"):
                print >>sys.stderr, "strand error in "+filename
                skippedlines += 1
            elif not (fields[3].isdigit() and fields[4].isdigit()):
                print >>sys.stderr, "non-number coordinates in "+filename
                skippedlines += 1
            else:
                                    
                yield GenomeRange( orgdb, fields[0],fields[3],fields[4],fields[6], name = featname, fastafile = seqfile, data = {"biotype":biotype, "source":genesource, "genename":genename})
            
def readbed(filename, orgdb="genome", seqfile= None):
    bedfile = None
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    skippedlines = 0
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split()
        if len(fields) > 2:
            if len(fields) < 5:
                strand = "+"
            else:
                strand = fields[5]
            if not (strand == "+" or strand == "-"):
                print >>sys.stderr, "strand error in "+filename
                skippedlines += 1
            elif not (fields[1].isdigit() and fields[2].isdigit()):
                print >>sys.stderr, "non-number coordinates in "+filename
                skippedlines += 1
            else:
                yield GenomeRange( orgdb, fields[0],fields[1],fields[2],strand, name = fields[3], fastafile = seqfile)
    
    if skippedlines > 0:
        print >>sys.stderr, "skipped "+str(skippedlines)+" in "+filename
def ifelse(arg, trueres,falseres):
    if arg:
        return trueres
    else:
        return falseres
        
        
def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)        
def issinglemapping(mapping):
    return mapping.mapq > 2        

def getbamrange(bamfile, chromrange = None, primaryonly = False, singleonly = False, maxmismatches = None, allowindels=True):
    bamiter = None
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            rname = bamfile.getrname(currline.rname)
            #print rname
            #need to fix this with cigar stuff
            #len(currline.pos)
            strand = "+"
            strand = ifelse(currline.is_reverse, '-','+')
            #not giving the reverse complement for now
            seq = currline.seq
            #print currline.cigar
            #print >>sys.stderr, currline.qname
            if 'SRR1508404.892272' ==  currline.qname:
                pass
                #print >>sys.stderr, "***"
                #print >>sys.stderr, currline.mapq
                #print >>sys.stderr, issinglemapped(currline)
            if primaryonly and not isprimarymapping(currline):
                continue
            if singleonly and not issinglemapping(currline):
                continue
            if strand == "-":
                pass
            #[("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            uniqueac = True
            uniqueamino = True
            uniquetrna = True
            #print >>sys.stderr, dir(currline)
            mismatches = None
            alignscore = None
            secondbestscore = None
            uniquemapping = False
            for currtag in currline.tags:
                if currtag[0] == "YA" and currtag[1] > 1:
                    uniqueac = False
                if currtag[0] == "YM" and currtag[1] > 1:
                    uniqueamino = False
                if currtag[0] == "YR" and currtag[1] > 1:
                    uniquetrna = False   
                if currtag[0] == "XM":
                    mismatches = currtag[1]
                if currtag[0] == "XS":
                    secondbestscore = float(currtag[1])
                if currtag[0] == "AS":
                    alignscore = float(currtag[1])
            if  secondbestscore is None or alignscore > secondbestscore:
                uniquemapping = True
            if maxmismatches is not None and currtag[1] > maxmismatches:
                continue  
            if not allowindels and len(currline.cigar) > 1:
                continue
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname , data = {"score":currline.mapq, "CIGAR":currline.cigar,"CIGARstring":currline.cigarstring, "seq":seq, "flags": currline.flag, "qual":currline.qual,"bamline":currline,'uniqueac':uniqueac,"uniqueamino":uniqueamino,"uniquetrna":uniquetrna,"uniquemapping":uniquemapping})
    except ValueError as err:
        #print>>sys.stderr, err
        #print>>sys.stderr, bamfile.name
        if chromrange is not None:
            #print >>sys.stderr, chromrange.chrom+":"+ str(chromrange.start)+"-"+str(chromrange.end) +" failed"
            pass
        
#'uniqueac':uniqueac,"uniqueamino":uniqueamino,"uniquetrna":uniquetrna})
def isuniquetrnamapping(read):
    return read.data["uniquetrna"]
def isuniqueaminomapping(read):
    return read.data["uniqueamino"]
def isuniqueacmapping(read):
    return read.data["uniqueac"]
def issinglemapped(read):
    return (read.data["score"] >= 2)        
def getpileuprange(bamfile, chromrange = None):
    bamiter = None
    if chromrange is not None:
        bamiter = bamfile.pileup(chromrange.chrom, chromrange.start, chromrange.end)
    else:
        bamiter = bamfile.pileup()
    
    for currpos in bamiter:
        readcounts = defaultdict(int)
        reference = 'N'
        for read in currpos.pileups:
            aln = read.alignment
            
            
            if not read.indel and not read.is_del:
                readcounts[aln.seq[read.qpos]] += 1
        yield read.qpos, readcounts
                
        
def getseqdict(genelist, faifiles = None):
    namedict = getnamedict(genelist)
    allorgs = set(currgene.dbname for currgene in genelist)
    dbdict = dict()
    fastafiles = dict()
    for currorg in allorgs:
        dbdict[currorg] = dict()
        
    for currgene in genelist:
        dbdict[currgene.dbname][currgene.name] = currgene
        if currgene.fastafile is not None:
            fastafiles[currgene.dbname] =  currgene.fastafile
        else:
            print >>sys.stderr, "no fasta file"
            sys.exit(1)
    seqdict = dict()
    
    for currorg in allorgs:
        
        if faifiles is not None:
            
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg], faindex = faifiles[currorg])
        else:
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg])
        seqdict.update(currseqs)
    return seqdict
    
def getnamedict(genelist):
    namedict = dict()
    for currgene in genelist:
        namedict[currgene.name] = currgene
    return namedict
        
        
def getseqs(fafile,rangedict, faindex = None):
    #print >>sys.stderr, rangedict
    
    if faindex is not None:
        try:
            faifile = fastaindex(fafile, faindex)
        except IOError as e:
            print >>sys.stderr, "Cannot read fasta file "+fafile
            print >>sys.stderr, "Ensure that file "+fafile +" exits and generate fastaindex "+faindex+" with samtools faidx"
            sys.exit(1)
        return faifile.getseqs(rangedict)
    genomefile = open(fafile, "r")
    reheader = re.compile( r"\>([^\s\,]+)")
    
    allseqs = defaultdict(str)
    currloc = 0
    
    for line in genomefile:
        line = line.rstrip("\n")
        currheader = reheader.match(line)
        if currheader: #sequence += line[chromstart - currloc:chromend - currloc]
            currseq = currheader.groups(1)[0]
            #print >>sys.stderr, currseq
            currloc = 0
        else:
            for currname, location in rangedict.iteritems():
                if currseq == location.chrom:
                    
                    chromstart = location.start
                    chromend = location.end
                    if location.dbname == 'eschColi_K12':
                        pass
                        #print >>sys.stderr, genomefile
                        #print >>sys.stderr, currseq+":"+str(chromstart)+"-"+str(chromend)+":"+str(currloc)
                    if currloc <= chromstart <= currloc + len(line) and currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:chromend - currloc]
                    elif currloc <= chromstart <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:]
                    elif currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[:chromend - currloc]
                    elif  currloc < chromstart <chromend < currloc + len(line):
                        pass
                    elif chromstart <= currloc < currloc + len(line) < chromend:
                        allseqs[currname] += line
            currloc += len(line)
    genomefile.close()
    finalseqs = dict()
    for currname in allseqs.iterkeys():
        #allseqs[currname] = allseqs[currname].upper()
        if (rangedict[currname].strand == "-"):
            seq = list(allseqs[currname].upper())
            seq.reverse()
            comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
            finalseqs[currname]  = ''.join(comp[base] for base in seq)
        else:
            finalseqs[currname] = allseqs[currname].upper()
    for currseq in rangedict.iterkeys():
        if currseq not in finalseqs:
            print >>sys.stderr, "No sequence extracted for "+rangedict[currseq].dbname+"."+rangedict[currseq].chrom+":"+str(rangedict[currseq].start)+"-"+str(rangedict[currseq].end)
    return finalseqs        
    
class fastaindex:
    def __init__(self, fafile, faifile):
        self.fafile = fafile
        fai = open(faifile)
        self.chromsize = dict()
        self.chromoffset = dict()
        self.seqlinesize = dict()
        self.seqlinebytes = dict()
        
        for line in fai:
            fields = line.split("\t")
            self.chromsize[fields[0]] = int(fields[1])
            self.chromoffset[fields[0]] = int(fields[2])
            self.seqlinesize[fields[0]] = int(fields[3])
            self.seqlinebytes[fields[0]] = int(fields[4])
    def getchrombed(self, dbname = 'genome'):
        for curr in self.chromsize.iterkeys():
            yield GenomeRange(dbname,curr,0,self.chromsize[curr],name=curr, strand = "+")
    def getseek(self, currchrom,loc):
        #print >>sys.stderr, (self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
        return self.chromoffset[currchrom] + loc + int(loc/(self.seqlinesize[currchrom]))*(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
    def getfullseqs(self, names):
        genomefile = open(self.fafile, "r")
        for currchrom in names:
            #print >>sys.stderr, currchrom+":"+str(self.chromsize[currchrom])
            genomefile.seek(self.getseek(currchrom,0))
            #seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]))
            seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]) - self.getseek(currchrom,0))
            #seq = seq.replace("\n","")
            yield currchrom, seq 
    def getseqs(self,  rangedict):
        genomefile = open(self.fafile, "r")
        #print >>sys.stderr, "(())"
        #print >>sys.stderr, self.chromsize
        allseqs = dict()
        for currname, currregion in rangedict.iteritems():
            currchrom = currregion.chrom
            #faskip = 
            #print >>sys.stderr, int(currregion.start/(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom]))
            genomefile.seek(self.getseek(currchrom,currregion.start))
            seq = genomefile.read(self.getseek(currchrom,currregion.end) - self.getseek(currchrom,currregion.start))
            seq = seq.replace("\n","")
            allseqs[currname] = seq
            #print >>sys.stderr, len(seq)
            #print >>sys.stderr, str(currregion.end - currregion.start)
            
        genomefile.close()
        finalseqs = dict()
        for currname in allseqs.iterkeys():
            #allseqs[currname] = allseqs[currname].upper()
            if (rangedict[currname].strand == "-"):
                seq = list(allseqs[currname].upper())
                seq.reverse()
                comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
                finalseqs[currname]  = ''.join(comp[base] for base in seq)
            else:
                finalseqs[currname] = allseqs[currname].upper()
        return finalseqs    
        
        
        
#the object I use for storing groups of genomeranges.  Useful for finding overlaps
class RangeBin:     
    def __init__(self,rangelist, binfactor = 10000):
        self.binfactor = binfactor
        self.bins = []
        self.length = 0
        for curr in rangelist:
            self.additem(curr)
    def __len__(self):
        return self.length
    def __iter__(self):
        for currbin in self.bins:
            #print currbin
            for currgene in currbin:
                yield currgene
    def additem(self, item):
        binstart = int(item.start / self.binfactor)
        #print "**"+str(binstart)
        binend = int(item.end / self.binfactor) + 1
        while (binstart + 2 >= len(self.bins)):
            self.bins.append(set())
        self.bins[binstart].add(item)
        self.length += 1
        #print self.bins[binstart]
    def getrange(self, item):
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins):
                for currrange in self.bins[i]:
                    if currrange.start >= item.start and currrange.end <= item.end:
                            yield currrange
                            
    def getbin(self, item):
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    yield currrange
    def getbinpos(self, item):
        for i in range(int(item / self.binfactor)-1,int(item / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    yield currrange
      
def revcom(sequence, rnamode = False):
    seq = list(sequence)
    seq.reverse()
    if rnamode:
        comp = {"A":"U","U":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K","-":"-",".":"."}
    else:
        comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K","-":"-",".":"."}
    return ''.join(comp[base] for base in seq)            