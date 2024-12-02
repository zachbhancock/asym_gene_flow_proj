#python bs

import os
import numpy as np
import msprime, pyslim
import tskit
import sys
import pandas as pd
import math

ts = tskit.load("ss_center_sym_50_to_1_10.trees")
ts = ts.simplify()
recap_ts = pyslim.recapitate(ts, recombination_rate=1e-8, ancestral_Ne=10000)
ts = msprime.sim_mutations(recap_ts, rate=1e-7, random_seed=12345)

pop1_inds = pyslim.individuals_alive_at(ts, 0, population=0)
pop2_inds = pyslim.individuals_alive_at(ts, 0, population=1)
pop3_inds = pyslim.individuals_alive_at(ts, 0, population=2)
pop4_inds = pyslim.individuals_alive_at(ts, 0, population=3)
pop5_inds = pyslim.individuals_alive_at(ts, 0, population=4)
pop6_inds = pyslim.individuals_alive_at(ts, 0, population=5)
pop7_inds = pyslim.individuals_alive_at(ts, 0, population=6)
pop8_inds = pyslim.individuals_alive_at(ts, 0, population=7)
pop9_inds = pyslim.individuals_alive_at(ts, 0, population=8)
pop10_inds = pyslim.individuals_alive_at(ts, 0, population=9)
    
groups = {
    'pop1' : np.random.choice(pop1_inds, size=5, replace=False),
    'pop2' : np.random.choice(pop2_inds, size=5, replace=False),
    'pop3' : np.random.choice(pop3_inds, size=5, replace=False),
    'pop4' : np.random.choice(pop4_inds, size=5, replace=False),
    'pop5' : np.random.choice(pop5_inds, size=5, replace=False),
    'pop6' : np.random.choice(pop6_inds, size=5, replace=False),
    'pop7' : np.random.choice(pop7_inds, size=5, replace=False),
    'pop8' : np.random.choice(pop8_inds, size=5, replace=False),
    'pop9' : np.random.choice(pop9_inds, size=5, replace=False),
    'pop10' : np.random.choice(pop10_inds, size=5, replace=False)
}

group_order = ['pop1', 'pop2', 'pop3', 'pop4', 'pop5', 'pop6', 'pop7', 'pop8', 'pop9', 'pop10']

sampled_nodes = [[] for _ in groups]
for j, k in enumerate(group_order):
   for ind in groups[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)

ind_nodes = []
ind_group = []
ind_ids = []
for j, group in enumerate(group_order):
   for ind in groups[group]:
      ind_ids.append(ind)
      ind_nodes.append(ts.individual(ind).nodes)
      ind_group.append(group_order[j])

nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(nind)]
ind_div = ts.divergence(ind_nodes, indexes=pairs)

x = []
for i in ind_ids:
   ind = ts.individual(i)
   label = f"tsk{ind.id}pop"
   x.append(label)

b = np.reshape(ind_div, (50,50))
panda_df = pd.DataFrame(data = b, columns = x)
panda_df.to_csv("ss_center_sym_50_to_1_10-pi.csv", sep=' ', index=True)

ind_Fst = ts.Fst(ind_nodes, indexes=pairs)

c = np.reshape(ind_Fst, (50,50))
panda2_df = pd.DataFrame(data = c, columns = x)
panda2_df.to_csv("ss_center_sym_50_to_1_10-Fst.csv", sep=" ", index=True)

indivlist = []
indivnames = []
with open("ss_center_sym_50_to_1_10.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for group in group_order:
     for i in groups[group]:
        indivlist.append(i)
        ind = ts.individual(i)
        vcf_label = f"tsk{ind.id}pop"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

with open("ss_center_sym_50_to_1_10.vcf", "w") as vcffile:
   ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

