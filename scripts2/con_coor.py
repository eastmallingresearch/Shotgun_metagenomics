#!/home/deakig/usr/local/bin/python
# coding: utf-8

from hirbin import *
from hirbin.parsers import *
import os
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys


def parseArgs(argv):
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Perform functional annotation of assembled metagenomics sequences, using TIGRFAM or PFAM.')
  parser.add_argument('-p','--maxOverlap',dest='max_acceptable_overlap',default=0.1,help='Max percentage of acceptable sequence overlap for bins at the same contig (if overlapping more than p percent, take the best scoring HMM-profile), default=%(default)s')
  parser.add_argument('-o',dest='output_file')
  parser.add_argument('-d',dest='protein_file')
  parser.add_argument('-m',dest='hmm_file')
  arguments=parser.parse_args(argv)
  return arguments


def main(protein_file,hmm_file,output_file,max_acceptable_overlap):
  convert_coordinates(protein_file,hmm_file,output_file,"",max_acceptable_overlap)
  print "Done"

if __name__=='__main__':
  arguments=parseArgs(sys.argv[1:])
  main(arguments.protein_file,arguments.hmm_file,arguments.output_file,arguments.max_acceptable_overlap)
