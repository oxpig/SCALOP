import numpy as np
import optics
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.cluster import DBSCAN
from DTW_data import DTWdata
from optparse import OptionParser
from sklearn.manifold import MDS

from matplotlib.widgets import Cursor
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scalop.utils import cutoffs
import os

class Dump(object):

	def __init__(self):

		self.db = None

		self.cluster_centers = None
		self.fname = ''

	def write(self):
		print('start writing cluster data.txt', self.fname)
		with open(os.path.join(resultsdir,"Cluster_data%s.txt" %(self.fname)), "w") as clusterfile:
			
			clusterfile.write("Loop_name\tSequence_with_anchors\tCluster\tCluster_center\n")

			for i in range(ran_points.shape[0]):

				loop_name = dtwdata.loop_names[i]

				Sequence = dtwdata.sequence_dict[loop_name]

				cluster_center_ind = "No"

				current_label = self.db.labels_[i]
				
				if self.cluster_centers[current_label] == i:

					cluster_center_ind = "Yes"

				clusterfile.write(loop_name + "\t" + Sequence + "\t" + str(current_label) + "\t" + cluster_center_ind + "\n")

def find_cluster_centers(db):

	labels = db.labels_

	unique_labels = set(labels)

	for label in unique_labels:

		if label == -1:

			continue

		label_indices = np.where(labels == label)[0]

		trimmed_matrix = D[label_indices,:][:,label_indices]

		medians = np.mean(trimmed_matrix, axis = 1)

		cluster_center_internal = np.argmin(medians)

		cluster_center_global = label_indices[cluster_center_internal]

		cluster_centers[label] = cluster_center_global

	Dump_callback.cluster_centers = cluster_centers

def cluster_data(threshold):

	db = DBSCAN(eps = threshold, min_samples = 5, metric = "precomputed")

	db.fit(D)

	Dump_callback.db = db

	return db

parser = OptionParser()

parser.add_option("-d", "--distfile",
                  action="store",
                  dest="distfile",
                  default=None,
                  help="File containing distance matrix")
parser.add_option("-f", "--filelist",
                  action="store", # optional because action defaults to "store"
                  dest="filelist",
                  default=None,
                  help="File containing loop file names",)
parser.add_option("-i", "--directory",
                  action="store", # optional because action defaults to "store"
                  dest="directory",
                  default=False,
                  help="Directory containing loop structures",)
parser.add_option("-t", "--thresholds",
                  action="store", # optional because action defaults to "store"
                  dest="thresholds",
                  default=cutoffs,
                  help="CDR-specific clustering thresholds",)

(options, args) = parser.parse_args()

print("Loading distance matrix")

dtwdata = DTWdata(distmat_filename = options.distfile, loop_names_filename = options.filelist, pdb_directory = options.directory)

D = dtwdata.distance_matrix

print("Projecting distance matrix")

np.random.seed(4718) # fixed for consistent visualisation of the OPTICS plot - kept for publication plots

mds = MDS(n_components = 2, dissimilarity = "precomputed")

ran_points = mds.fit_transform(D)

cluster_centers = {-1: None}

print("Plotting")

Dump_callback = Dump()

threshold = options.thresholds[options.directory[-2:]]
db = cluster_data(threshold)
find_cluster_centers(db)
resultsdir=os.path.dirname(os.path.abspath(*options.distfile))
Dump_callback.fname = options.directory[-2:]+'_'+str(threshold)+'A'
Dump_callback.write()



