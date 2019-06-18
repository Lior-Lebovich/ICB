import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
import kabuki
from kabuki.analyze import gelman_rubin
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
import hddm
import sys

# Adding a monkey in order to save the models dbs:
import pickle
# Get around a problem with saving regression outputs in Python 3
def savePatch(self, fname):
    with open(fname, 'wb') as f:
        pickle.dump(self, f)
hddm.HDDM.savePatch = savePatch

docdir = 'D:\\all hddm analysis\\BL\\informative priors 40K\\ic'
saved_model_name = 'BL_ic_no_tttv'
os.chdir(docdir)

models = []

for i in range(12):
    print(i)
    m = pickle.load(open('submodel'+str(i), 'rb'))
    models.append(m)
print('Thats it!!')

gelman_rubin_stats = gelman_rubin(models)
with open('gelman_rubin_stats.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in gelman_rubin_stats.items():
       writer.writerow([key, value])
 
      
# The combined model:
running_model = kabuki.utils.concat_models(models)
# Saving the model:
running_model.savePatch(saved_model_name)
# Saving the traces:
running_model.get_traces().to_csv('traces.csv') 
# Saving the model's stats:
running_model_stats = running_model.gen_stats()
running_model_stats.to_csv('stats_csv_'+saved_model_name+'.csv')


# For model comparison:
running_model.dic_info
running_model.logp


# simulate data:
ppc_data = hddm.utils.post_pred_gen(running_model, append_data=True, 
                                    progress_bar=True, samples=100)
ppc_data.to_csv('ppc_data.csv')

# load raw data:
data = hddm.load_csv('D:\\all hddm analysis\\verticalBL_imp.csv')
data = data[data.rt<=3]

# compare simulated and raw data:
ppc_compare = hddm.utils.post_pred_stats(data, ppc_data)
ppc_compare.to_csv('ppc_stats.csv')

ppc_compare_no_compare = hddm.utils.post_pred_stats(data, ppc_data, call_compare=False)
ppc_compare_no_compare.to_csv('ppc_stats_no_compare.csv')


# Plotting posterior-predictive:
os.chdir(docdir)
fig = plt.figure()
running_model.plot_posterior_predictive(columns=10,
                                           figsize=(100, 100), save=True)
plt.savefig('posteriors predictive.pdf')

