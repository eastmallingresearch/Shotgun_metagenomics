#!/home/deakig/usr/local/bin/python
from hirbin import *
from hirbin.parsers import *
import os
from os import mkdir
import sys

cov_file=str(sys.argv[1])
tab_file=str(sys.argv[2])

print "cov file: " + cov_file  +" tab file: " + tab_file

parseCoverageBed(cov_file,tab_file)