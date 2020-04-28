import numpy as np
import os

class DTWdata(object):

	res_dict = {
			 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
			 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
			 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
			 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'TYS': 'Y', 
			 'MET': 'M', 'MSE': 'M'
	}

	def __init__(self, distmat_filename, loop_names_filename, pdb_directory):

		self.distmat_filename = distmat_filename

		self.loop_names_filename = loop_names_filename

		self.distance_matrix = np.loadtxt(distmat_filename)

		self.pdb_directory = pdb_directory

		loop_names_fileobject = open(loop_names_filename, "r")

		self.loop_names = [line.strip().split('/')[-1] for line in loop_names_fileobject]

		loop_names_fileobject.close()

		self.loop_indices = dict( list(zip ( self.loop_names, np.arange(len(self.loop_names)) )) )

		self.read_sequences()

	def read_sequence(self, pdb_filename):

		Sequence = ''

		with open(pdb_filename ,'r') as pdb_file:
			
			prev_index = ''
			
			for line in pdb_file:
				
				if (line[:6] != "REMARK") and line[:3] != "TER" and line[:3] != "END":

					cur_index = line[22:27]
					
					if(cur_index != prev_index):
					
						residue_name = line[17:20]
						
						Sequence += self.res_dict[residue_name] if residue_name in self.res_dict else ''
						
						prev_index = cur_index

		return Sequence

	def read_sequences(self):

		self.sequence_dict = {}

		#preloaded_data = list(os.walk(self.pdb_directory))

		#number_of_files = len([filename for root in preloaded_data for filename in root[2] if filename[-4:] == ".pdb"])

		number_of_files = len( os.listdir(self.pdb_directory))

		i = 0
	
		file_list = os.listdir(self.pdb_directory)
		
		for filename in file_list:

			seq = self.read_sequence(self.pdb_directory + '/' + filename)

			self.sequence_dict[filename] = seq
			
			i+=1

			if i % 100 == 0:

				print("Reading sequences from pdbs %.2f complete\r" % (i * 100.0 / number_of_files))

	def __repr__(self):

		return "DTW data from files " +  self.distmat_filename + " and " + self.loop_names_filename

	def __getitem__(self, loop_names):

		if len(loop_names) == 2:

			loop_name_1 = loop_names[0]

			loop_name_2 = loop_names[1]

			if loop_name_1 not in self.loop_indices:

				raise TypeError("Loop name " + loop_name_1 + " not found")

			if loop_name_2 not in self.loop_indices:

				raise TypeError("Loop name " + loop_name_2 + " not found")

			index_1 = self.loop_indices[loop_name_1]

			index_2 = self.loop_indices[loop_name_2]

			return self.distance_matrix[index_1, index_2]

		else:

			raise NameError("This function takes exactly 2 arguments")

	def __iter__(self):

		for loop_name in self.loop_names:

			yield loop_name



