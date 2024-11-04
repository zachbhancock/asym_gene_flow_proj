###FEEMS#####

# base
import numpy as np
import pandas as pd
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz

data_path = pkg_resources.resource_filename("feems", "data/")
(bim, fam, G) = read_plink("{}feems_bed".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

coord = np.loadtxt("{}/haust_locations.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/haust_locations.outer".format(data_path))  # outer coordinates
grid_path = "{}/grid_100.shp".format(data_path)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

projection = ccrs.EquidistantConic(central_longitude=-77.07186765, central_latitude=35.22513235)

sp_graph.fit(lamb = 2.0)

feems_nodes = sp_graph.nodes
pd.DataFrame(feems_nodes).to_csv('feems_nodes.csv', header=False, index=False)

feems_node_pos = sp_graph.node_pos
pd.DataFrame(feems_node_pos).to_csv('feems_node_pos.csv', header=False, index=False)

feems_edges = sp_graph.edges
pd.DataFrame(feems_edges).to_csv('feems_edges.csv', header=False, index=False)

feems_w = sp_graph.w
pd.DataFrame(feems_w).to_csv('feems_w.csv', header=False, index=False)
