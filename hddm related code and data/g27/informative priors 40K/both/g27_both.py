
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import kabuki
from kabuki.analyze import gelman_rubin
import hddm
print(hddm.__version__)
# Python version:   
import sys
print (sys.version)
print (pd.__version__)

# Adding a monkey in order to save the models dbs:
import pickle
# Get around a problem with saving regression outputs in Python 3
def savePatch(self, fname):
    with open(fname, 'wb') as f:
        pickle.dump(self, f)
hddm.HDDM.savePatch = savePatch


data = hddm.load_csv("/ems/elsc-labs/loewenstein-y/lior.lebovich/DD/g27/data_csv_g27_dl0_3columns.csv")
# data = data[data.subj_idx!=125]
docdir = "/ems/elsc-labs/loewenstein-y/lior.lebovich/DD/g27/both"
os.chdir(docdir)

i = sys.argv[1]
print(i)
m = hddm.HDDM(data, include=('v','z'), 
			  p_outlier=.05, informative=True)
m.find_starting_values()
m.sample(40000, burn=20000, thin=5, dbname='db'+str(i), db='pickle')
m.savePatch('submodel'+str(i))

print('Thats it!!')
