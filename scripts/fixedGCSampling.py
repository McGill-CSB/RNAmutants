#!/usr/bin/python 
import sys
import os
import numpy
import random
import math
import re

DEBUG = False

NUM_SAMPLES = 10000
GC_WEIGHT = 1.0
TARGET_GC = 0.5
LIMIT_MUTS = sys.maxint
SAMPLES_HEADER = ">> Sampling "
SAMPLES_NEW_MUTS_HEADER = "> sampling "
# Indicate the location of RNAmutants here
path_to_RNAmutants = "."

COMMAND_LINE = "%s/RNAmutants -l %s/lib/ -M %s -C -f %s -m %s --sample-file %s > %s"


def buildSampleFile(path, n, num=NUM_SAMPLES,kmin=0,kmax=sys.maxint):
  sfile = open(path,"w")
  for i in range(n+1):
    if kmin<=i<=kmax:
      sfile.write("%s\t%s\n"%(i,num))
  sfile.close()
  return path

def buildInputFile(path,seq,mutMask=None,structConstraints=None):
  if mutMask is None:
    mutMask = "o"*len(seq)
  if structConstraints is None:
    structConstraints = "?"*len(seq)
  sfile = open(path,"w")
  sfile.write("%s\n"%(seq))
  sfile.write("%s\n"%(mutMask))
  sfile.write("%s\n"%(structConstraints))
  sfile.close()
  return path

def buildWeightFile(path,gcWeight=1.0):
  sfile = open(path,"w")
  bases = ['A','C','G','U']
  for r in bases:
    for c in bases:
      coeff = 1.0
      if r in ['C','G']:
        coeff /= gcWeight
      if c in ['C','G']:
        coeff *= gcWeight
      sfile.write("%s "%coeff)
    sfile.write("\n")
  sfile.close()
  return path

  
def runRNAMutants(seq,mutMask=None,structConstraints=None,gcWeight=GC_WEIGHT,numSamples=NUM_SAMPLES,kmin=0,kmax=sys.maxint):
  inputPath = "test.faa"
  samplePath = "samples.txt"
  outputPath = "out-%.2f.txt"%gcWeight
  weightPath = "weights.txt"
  kmax = min(kmax,len(seq))
  buildInputFile(inputPath,seq,mutMask,structConstraints)
  buildSampleFile(samplePath,len(seq),numSamples,kmin,kmax)
  buildWeightFile(weightPath,gcWeight)
  cmd = COMMAND_LINE%(path_to_RNAmutants,path_to_RNAmutants,weightPath,inputPath,kmax,samplePath,outputPath)
  os.system(cmd)
  results = parseSampleSet(outputPath,kmin,kmax)
  return results

def parseSampleSet(path,kmin,kmax):
  result = {}
  currentSeq = ""
  inSampling = False
  seqLine = True
  muts = 0
  mutre = re.compile(r"""(\d+)""",re.VERBOSE)
  try:
    for l in open(path):
      if l.startswith(SAMPLES_HEADER):
        inSampling = True
      elif inSampling:
        if l.startswith(SAMPLES_NEW_MUTS_HEADER):
          m = mutre.findall(l)
          muts = int(m[1])
        else:
          if seqLine:
            s = l[:-1].strip()
            data = s.split()
            currentSeq = data[0]
            energy = float(data[1][1:-1])
          else:
            str = l[:-1].strip()
            if muts not in result:
              result[muts] = []
            result[muts].append((currentSeq,str,energy))
          seqLine = not seqLine
  except Exception, e:
    print e,l
  return result

def countGC(seq):
  totGC = 0 
  for c in seq:
    if c in ['G','C','g','c']:
      totGC += 1
  return totGC

def analyzeGCResults(results):
  freqs = {}
  for muts in results:
    totLength = 0
    totGC = 0
    for (seq,str,energy) in results[muts]:
      totGC += countGC(seq)
      totLength += len(seq)
    freqs[muts] = float(totGC)/float(totLength)
  return freqs

def addFilteredSampleSet(samples,finalSamples,lowerBound=-sys.maxint,upperBound=+sys.maxint,kmin=0,kmax=sys.maxint):
  for muts in samples:
    if kmin<=muts<=kmax:
      if muts not in finalSamples:
        finalSamples[muts] = []
      for (seq,str,energy) in samples[muts]:
        if lowerBound <= countGC(seq) <= upperBound:
          finalSamples[muts].append((seq,str,energy))

MAGIC_RATIO = 3

def displayExps(finalSamples,weightsToFreqs,kmin,kmax):
  weights = weightsToFreqs.keys()
  weights.sort()
  for muts in finalSamples:
    if kmin<=muts<=kmax:      
      s = ""
      for w in weights:
        s += "%s.2f->%.2f "%(w,weightsToFreqs[w][muts])
      print "  B[%s]=%s\t[%.2f-%.2f]\t%s"%(muts,len(finalSamples[muts]),weightsToFreqs[weights[0]][muts],weightsToFreqs[weights[-1]][muts],s)  

def checkBucketsFull(finalSamples,targetNumSamples,initSeq,lowerBoundGC,upperBoundGC,kmin,kmax):
  result = True
  for muts in  finalSamples:
    if (len(finalSamples[muts]) < targetNumSamples) and (kmin<=muts<=kmax) and isBucketFillable(initSeq,muts,lowerBoundGC,upperBoundGC):
      result = False
  return result

def isBucketFillable(seq,muts,lbGC,ubGC):
  lbMuts = countGC(seq)-muts
  ubMuts = countGC(seq)+muts
  overlap =  not ( (lbMuts<=ubMuts<lbGC<=ubGC)
              or (lbGC<=ubGC<lbMuts<=ubMuts))
  return overlap


def displayWarningUnfillableBuckets(initSeq,lowerBoundGC,upperBoundGC,kmin,kmax):
  for invMuts in  range(len(initSeq)+1):
    muts = len(initSeq)+1-invMuts
    if kmin<= muts<=kmax:
      if not isBucketFillable(initSeq,muts,lowerBoundGC,upperBoundGC):
        sys.stderr.write("Warning: Targetted GC%% unreachable within <=%s mutation(s)\n"%muts)
        return

def newWeight(seq,w,finalSamples,targetGC, targetNumSamples,weightsToFreqs, lowerBoundGC,upperBoundGC,kmin,kmax):
  # Figuring out the worst "bucket"
  worst = -1
  worstOccupancy = sys.maxint
  for muts in  finalSamples:
    if kmin<=muts<=kmax:
      occupancy = float(len(finalSamples[muts]))/float(targetNumSamples)
      if (occupancy<worstOccupancy) and isBucketFillable(seq,muts,lowerBoundGC,upperBoundGC):
        worst = muts
        worstOccupancy = occupancy
  if DEBUG: print "Worst=",worst,worstOccupancy
  # Building (weight->Freqs) list for worst, sorted on frequency
  exps = []
  for w in weightsToFreqs:
    numGC = weightsToFreqs[w][worst]
    exps.append((numGC,w))
  exps.sort()
  # No experiment yet, try current weight
  if len(exps)==0:
    return w
  if DEBUG: print "  Wanted:",targetGC
  # Targetted GC frequency smaller than smallest GC frequency achieved so far, half smallest weight
  if targetGC<exps[0][0]:
    if DEBUG: print "  Smaller than smallest"
    return exps[0][1]/2.0
  # Targetted GC frequency greater than greatest GC frequency achieved so far, double biggest weight
  if targetGC>exps[-1][0]:
    if DEBUG: print "  Greater than greatest"
    return exps[-1][1]*2.0
  # Targetted GC frequency 'squeezed' between two values, find them and go dichotomically
  if DEBUG: print "    exps"
  for i in range(len(exps)-1):
    if DEBUG: print "  Greater than greatest"
    (numGCLow,wLow) = exps[i]
    (numGCHigh,wHigh) = exps[i+1]
    if numGCLow <= targetGC <= numGCHigh:
      return (wLow+wHigh)/2.0
  return None


def fillBuckets(seq,targetGC,numSamples, tolerance = 0.0, kmin=0,kmax=sys.maxint):
  """Samples [numSamples] elements within each class of mutations (buckets) for
  a sequence [seq] having a proportion [targetGC] of GC. A [tolerance] can be specified,
  in which case the targetted proportion is (targetGC-tolerance,targetGC+tolerance)
  """
  finalSamples = {}
  n = len(seq)
  lowerBoundGC = math.floor(n*(targetGC-tolerance))
  upperBoundGC = math.floor(n*(targetGC+tolerance))
  w = 1.0
  bucketsFull = False
  weightToFreqs = {}
  displayWarningUnfillableBuckets(seq,lowerBoundGC,upperBoundGC,kmin,kmax)
  while not bucketsFull:
    samples = runRNAMutants(seq,None,None,w,MAGIC_RATIO*numSamples,kmin,kmax)
    weightToFreqs[w] = analyzeGCResults(samples)
    addFilteredSampleSet(samples,finalSamples,lowerBoundGC,upperBoundGC,kmin,kmax)
    bucketsFull = checkBucketsFull(finalSamples,numSamples,seq,lowerBoundGC,upperBoundGC,kmin,kmax)
    displayExps(finalSamples,weightToFreqs,kmin,kmax)
    if (not bucketsFull):
      wb = w
      w = newWeight(seq,w,finalSamples,targetGC,numSamples,weightToFreqs,lowerBoundGC,upperBoundGC,kmin,kmax)
      print "New weight:\n  %.3f \t ->\t %.3f"%(wb,w)
  return finalSamples

def printSamples(samples,kmin,kmax):
  for i in samples:
    if kmin<=i<=kmax:
      print "> Sampled %s sequence(s) and secondary structure(s) with %s mutations"%(len(samples[i]),i)
      for (seq,str,E) in samples[i]:
        print "  %s\t%s"%(seq,E)
        print "  %s"%str

def printUsage(cmd):
  print """\nCalls RNAMutants repeatedly to estimate suitable weights.\n
Usage: %s seq [opts]
Options:
  -n <n>   - Samples n structures each time
  -g <x>   - Targets a given C+G-content x (0<x<1)
  -e <y>   - Allows for a relative tolerance y on C+G-content
  -k n1 n2 - Only fills buckets from n1 to n2
  """ % (cmd)
  sys.exit(1)


if __name__ == '__main__':
  if len(sys.argv)==1:
    printUsage(sys.argv[0]);
    sys.exit(1);
  seq = sys.argv[1]
  i = 2
  numSamples = NUM_SAMPLES
  gcWeight = GC_WEIGHT
  targetGC = TARGET_GC
  tolerance = 0.0
  kmin = 0
  kmax = sys.maxint
  while i<len(sys.argv):
    if sys.argv[i]=="-k":
      i+=1
      kmin = int(sys.argv[i])
      i+=1
      kmax = int(sys.argv[i])
      print "a..."
    if sys.argv[i]=="-d":
      i+=1
      gcWeight = float(sys.argv[i])
    elif sys.argv[i]=="-n":
      i+=1
      numSamples = int(sys.argv[i])
    elif sys.argv[i] in ["-g"]:
      i+=1
      targetGC = float(sys.argv[i])
    elif sys.argv[i] in ["-e"]:
      i+=1
      tolerance = float(sys.argv[i])
    i += 1
  printSamples(fillBuckets(seq,targetGC,numSamples,tolerance,kmin,kmax),kmin,kmax)
