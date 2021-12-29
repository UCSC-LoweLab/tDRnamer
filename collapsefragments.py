#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import subprocess
from tdrdbutils import *
from distutils.spawn import find_executable
import time

import argparse



from collections import defaultdict

def readmultifastq(fqfile, fullname = False):
    #print chrom+":"+ chromstart+"-"+ chromend
    if fqfile == "stdin":
        fqfile = sys.stdin
    elif fqfile.endswith(".gz"):
        fqfile = gzip.open(fqfile, "rb")
    else:
        fqfile = open(fqfile, "r")
    currloc = 0
    currseq = None
    sequence = ""
    quality = ""
    if fullname:
        reheader = re.compile(r"\@(.+)")
    else:
        reheader = re.compile(r"\@([^\s\,]+)")
    qualheader = re.compile(r"\+([^\s\,]*)")
    readqual = False
    for line in fqfile:
        #print line
        line = line.rstrip("\n")
        seqheader = reheader.match(line)
        qheader = qualheader.match(line)
        if readqual:
            quality += line
            if len(quality) == len(sequence):
                yield currseq, sequence, quality
                readqual = False
        elif seqheader:
            currseq = seqheader.groups(1)[0]
            sequence = ""
            quality = ""
        elif qheader and readqual == False:
            readqual = True
            pass
        else:
            sequence += line
            
class prunedict:
    def __init__(self, trimcutoff = 10):
        self.counts = defaultdict(int) 
        self.trimcutoff = trimcutoff
        self.totalkeys = 0
        self.trimmed = 0
    def trimold(self):
        #print >>sys.stderr, "**"
        if self.trimcutoff is not None:
            newdict = defaultdict(int) 
            trimmed = 0
            for curr in self.counts.iterkeys():
                if self.counts[curr] > self.trimcutoff:
                    newdict[curr] = self.counts[curr]
                else:
                    trimmed += 1
            self.trimmed += trimmed
            print >>sys.stderr, str(trimmed)+"/"+str(self.totalkeys)+" at "+str(self.trimcutoff)
            self.totalkeys = len(self.counts.keys())
            self.counts = newdict
    def trim(self):
        #print >>sys.stderr, "**"
        if self.trimcutoff is not None:
            newdict = defaultdict(int) 
            trimmed = 0
            for curr in self.counts.iterkeys():
                if self.counts[curr] > self.trimcutoff:
                    newdict[curr] = self.counts[curr]
                else:
                    trimmed += 1
            self.trimmed += trimmed
            self.totalkeys = len(self.counts.keys())
            print >>sys.stderr, str(trimmed)+"/"+str(self.totalkeys)+" at "+str(self.trimcutoff)
            
            self.counts = newdict
        
    def __getitem__(self, key):
        return self.counts[key]
    def __setitem__(self, key, count):
        if key not in self.counts:
            self.totalkeys += 1
        self.counts[key] = count

    def iterkeys(self):
        return self.counts.keys()
            


parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')

parser.add_argument('--skipprune', action="store_true", default=False,
                   help='Sum samples that are replicates')
parser.add_argument('--trimcutoff', type=int, default=10,
                   help='margin to add to feature coordinates')
parser.add_argument("fastqfile")

args = parser.parse_args()
argdict = vars(args)

trimcutoff = None
skipprune = argdict["skipprune"]
fastqfile = argdict["fastqfile"]
if not skipprune:
    trimcutoff = int(argdict["trimcutoff"])

#skipprune = True
mincounts = 0
#maxseqs = 100000 
seqcount = None

seqcount = prunedict(trimcutoff = trimcutoff)
total = 0
for name, seq, qual in readmultifastq(fastqfile):
    seqcount[seq] += 1
    total += 1
    #if total % 100000 == 0:
        #print >>sys.stderr, str(total)

if not skipprune:        
    seqcount.trim()
#print >>sys.stderr, "results: "+str(len(seqcount.counts.keys()))+"/"+str(total)+":"+str(1.*len(seqcount.counts.keys())/total)
#print >>sys.stderr, "totaltrimmed: "+str(seqcount.trimmed)
#print >>sys.stderr, "finalcutoff: "+str(seqcount.trimcutoff)
#print >>sys.stderr, "maxkeys: "+str(seqcount.maxkeys)

for i, curr in enumerate(seqcount.iterkeys()):
    if seqcount[curr] > mincounts:
        print ">frag"+str(i+1)+":"+str(seqcount[curr])
        print curr