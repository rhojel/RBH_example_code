"""
This file conducts a molecular dynamics simulation for an input rvc file and outputs an rvc file and an erg file 
every 10 time steps. The output rvc file can be used to visualize the molecular dynamics simulation once it is converted 
to a crd format. The output erg file tracks the bond energy, nonbond energy, kinetic energy, and total energy for every 10 timesteps.

Usage: python MyMD.py --iF input.rvc <other params>

Richard Hojel
"""

import numpy as np 
import argparse

class MDSim(object):
	def __init__(self, **kwargs):
		"""
		Mothership class for the MD simulation. 

		Uses argparse to read in all the arguments from the command line. For all missing arguments
		the argument is set to a default (except for the input filename which is required)
		"""
		parser = argparse.ArgumentParser()
		parser.add_argument('--iF', required=True,  help='input filename') # input file
		parser.add_argument('--kB', type=float, nargs='?', default=40000.0, help='kb value') # same as kb
		parser.add_argument('--kN', type=float, nargs='?', default=400.0 , help='kn value') # same as kn
		parser.add_argument('--nbCutoff', type=float, nargs='?', default=0.50, help='non-bonded atome distance cutoff') # distance within atoms to consider nonbonded interactions
		parser.add_argument('--m', type=float, nargs='?', default=12.0, help='atom mass') # atom mass (constant for all)
		parser.add_argument('--dt', type=float, nargs='?', default=0.001, help='length of time step')  # length of time step
		parser.add_argument('--n', type=int, nargs='?', default=1000, help='number of iterations') # number of time steps
		parser.add_argument('--out', nargs='?',  help='prefix of output filename')  # prefix of the filename to output
		self.args = parser.parse_args()
		
		# Sets default --out to iF(stripped of '.rvc')
		if self.args.out == None:
			self.args.out = self.args.iF[:-4] + '_out'
		else:
			self.args.out = self.args.out + '_out'

	def parse_rvc_file(self, input_file):
		"""
		Parses a RVC (positions/velocities/connectivities) file. The RVC file is read into three global objects:
			- locations: a numpy array with each row representing an atom (row = 0 for atom 1) and three columns store 
			  the atom's position (x,y,z)
			- velocities: a numpy array with each row representing an atom (row = 0 for atom 1) and three columns store 
			  the atom's velocity (Vx,Vy,Vz)
			- connections: a dictionary that maps an atom to a set of bonded atoms (ex. 1:set(2,3,4))

		Parameters: 
		-----------
		input_file : name of input RVC file

		Returns: 
		--------
		first_line : the first line of the RVC file to be used when creating the output RVC file
		"""
		locations = []
		velocities = []
		self.connections = {}
		with open(input_file, 'r') as input_f:
			first_line = input_f.readline().split()
			for line in input_f:
				atom_info = line.split()
				locations.append((float(atom_info[1]),float(atom_info[2]),float(atom_info[3])))
				velocities.append((float(atom_info[4]),float(atom_info[5]),float(atom_info[6])))
				self.connections[int(atom_info[0])] = set()
				for connect in atom_info[7:]:
					self.connections[int(atom_info[0])].add(int(connect))
					
		self.locations = np.array(locations)
		self.velocities = np.array(velocities)
		return first_line

	def compute_reference_lengths(self):
		"""
		Computes the reference distance for bonded and nonbonded atoms. The nonbonded atoms are determined using the 
		nbCutoff parameter. The reference distances are stored in two global dictionaries (bond and nonbond). The dictionaries
		map an atom to another dictionary with all bonded or nonbonded atoms as keys. The second dictionary maps to the length
		of the given bond/nonbond (ex. 1:2:0.4360 where atom : bond/nonbond atom : bond/nonbond length)

		Parameters: 
		-----------
		input_file : name of input RVC file

		Returns: 
		--------
		None
		"""
		self.ref_dist_bond = {}
		self.ref_dist_nonbond = {}
		for atom in self.connections:
			self.ref_dist_bond[atom] = {}
			self.ref_dist_nonbond[atom] = {}
			calc_distances = np.sqrt(((self.locations - self.locations[atom-1])**2).sum(axis = 1))
			for ref_atom in range(len(calc_distances)):
				if (ref_atom+1) in self.connections[atom]:
					self.ref_dist_bond[atom][ref_atom+1] = calc_distances[ref_atom]
				elif calc_distances[ref_atom] < self.args.nbCutoff:
					if atom != (ref_atom+1):
						self.ref_dist_nonbond[atom][ref_atom+1] = calc_distances[ref_atom]

	def do_verlet_iteration(self):
		"""
		Performs one iteration of the Velocity Verlet integrator. 

		Returns: 
		--------
		E_b : the potential energy of all bonded atoms 
		E_nb : the potential energy of all nonbonded atoms
		E_k : the kenetic energy of each atom 
		E_tot : the total energy (E_b+E_nb+E_k)
		"""
		velocities_half_dt = self.velocities + 0.5 * (1/self.args.m) * self.forces * self.args.dt
		self.update_positions(velocities_half_dt)
		
		E_b, E_nb = self.updated_PE_and_forces()

		self.update_velocities(velocities_half_dt)
		E_k = (0.5 * self.args.m * ((self.velocities)**2).sum(axis = 1)).sum()

		return E_b,E_nb,E_k,E_b+E_nb+E_k

	def update_positions(self, velocities_half_dt):
		"""
		Updates atom positions according to the Velocity Verlet integrator.

		Parameters: 
		-----------
		velocities_half_dt : the velocities calculated for dt/2
		"""
		self.locations += velocities_half_dt * self.args.dt

	def updated_PE_and_forces(self):
		"""
		Computes the forces on the atoms according to the Velocity Verlet integrator. Computes the total potential energy for 
		both bonded and nonbonded atom pairs. 

		Returns: 
		--------
		E_b : the potential energy of all bonded atoms 
		E_nb : the potential energy of all nonbonded atoms
		"""
		# Sets energy and forces to zero
		E_b = 0
		E_nb = 0
		self.forces = np.zeros((len(self.connections),3))

		# Tracks PE values already computed (to avoid double counting)
		bond_E_computed = set()
		nonbond_E_computed = set()

		for atom in self.connections:
			distance = np.sqrt(((self.locations - self.locations[atom-1])**2).sum(axis = 1))
			for bond in self.ref_dist_bond[atom]:
				self.forces[atom-1] += (self.args.kB * (distance[bond-1]-self.ref_dist_bond[atom][bond])) * ((self.locations[bond-1] - self.locations[atom-1])/distance[bond-1])
				if (atom,bond) not in bond_E_computed:
					E_b += 0.5 * self.args.kB * (distance[bond-1]-self.ref_dist_bond[atom][bond])**2
					bond_E_computed.add((bond,atom))
			for nonbond in self.ref_dist_nonbond[atom]:
				self.forces[atom-1] += (self.args.kN * (distance[nonbond-1]-self.ref_dist_nonbond[atom][nonbond])) * ((self.locations[nonbond-1] - self.locations[atom-1])/distance[nonbond-1])
				if (atom,nonbond) not in nonbond_E_computed:
					E_nb += 0.5 * self.args.kN * (distance[nonbond-1]-self.ref_dist_nonbond[atom][nonbond])**2
					nonbond_E_computed.add((nonbond,atom))

		return E_b, E_nb

	def update_velocities(self, velocities_half_dt):
		"""
		Updates atom velocities according to the Velocity Verlet integrator.

		Parameters: 
		-----------
		velocities_half_dt : the velocities calculated for dt/2
		"""
		self.velocities = velocities_half_dt + 0.5 * (1/self.args.m) * self.forces * self.args.dt

	def write_rvc_output(self, iter_num, file_handle, first_line, E_tot, write_append):
		"""
		Writes current atom positions and velocities, as well as the iteration
		number, to an output file given by file_handle.

		Parameters: 
		-----------
		iter_num : int, current iteration number.
		file_handle : handle to *.rvc file opened for writing
					  (*NOT* file name)
		first_line : the first line of the input RVC file to be used when creating first line of the file
		E_tot : the total energy at the given iter_num
		write_append : indicates whether the file should be opened in writing mode or append mode 
					   (write is used to overwrite a file at the start, otherwise the information is appended to the file)

		Returns:
		--------
		None        
		"""
		with open(file_handle + '.rvc', write_append) as output_f:
			if(iter_num == 0):
				output_f.write(first_line[0]+' '+first_line[1]+' kB='+'%.1f' %self.args.kB+' kN='+'%.1f' %self.args.kN+
					' nbCutoff='+'%.2f' %self.args.nbCutoff+' dt='+'%.4f' %self.args.dt+' mass='+'%.1f' %self.args.m+' '+
					first_line[-1]+'\n')
			else:
				output_f.write('#At time step ' + str(iter_num) + ',energy = ' + '%.3f' % E_tot +'kJ\n')

			for atom in self.connections:
				output_f.write(str(atom) + '\t' + '\t'.join(map("{:.4f}".format, self.locations[atom-1])) + 
							   '\t' + '\t'.join(map("{:.4f}".format, self.velocities[atom-1])) + 
							  '\t' + '\t'.join(map(str,self.connections[atom])) + '\n')
		
	def write_erg_output(self, iter_num, file_handle, E_k, E_b, E_nb):
		"""
		Writes energy statistics (kinetic energy, potential energy of 
		bonded interactions, potential energy of nonbonded interactions,
		and the sum of the foregoing energies - E_tot) as well as the iteration
		number to an output file given by file_handle.

		Parameters: 
		-----------
		iter_num : int, current iteration number.
		file_handle : handle to *.erg file opened for writing 
					 (*NOT* file name)
		E_b : the potential energy of all bonded atoms 
		E_nb : the potential energy of all nonbonded atoms
		E_k : the kenetic energy of each atom 

		Returns:
		--------
		None
		"""
		# Overwrites file at the beginning, otherwise appends to the file
		if iter_num == 10:
			write_append = 'w'
		else:
			write_append = 'a'

		with open(file_handle + '.erg', write_append) as output_f:
			if iter_num == 10:
				output_f.write('# step'+'\t'+'E_k'+'\t'+'E_b'+'\t'+'E_nB'+'\t'+'E_tot\n')
			output_f.write(str(iter_num) + '\t')
			output_f.write('%.1f' % E_k + '\t')
			output_f.write('%.1f' % E_b + '\t')
			output_f.write('%.1f' % E_nb + '\t')
			output_f.write('%.1f' % (E_k+E_b+E_nb) + '\n')

	def run_md_sim(self):
		"""
		Runs the MD simulation. 
		"""
		first_line = self.parse_rvc_file(self.args.iF)
		self.compute_reference_lengths()
		self.write_rvc_output(0, self.args.out, first_line, 0, 'w')	# writes rvc file for time_step = 0 (copy of input rvc)
		
		intial_E = 0	# variable to track the intial energy
		self.forces = np.zeros((len(self.connections),3))	# numpy array to store forces on each atom (Fx,Fy,Fz)

		for step in range(1,self.args.n+1):
			E_b,E_nb,E_k,E_tot = self.do_verlet_iteration()
			
			if step == 1:
				intial_E = E_tot

			# Throws an error if the energy increases ten fold or decreases ten fold because that is a sign 
			# of a unstable simulation 
			if ((intial_E / E_tot >= 10) | (E_tot / intial_E >= 10)):
				self.write_rvc_output(step, self.args.out, first_line, E_tot,'a')
				self.write_erg_output(step, self.args.out, E_k, E_b, E_nb)
				raise ValueError('An overflow has occurred at iteration number ' + str(step) + '. Check the parameters.')

			# Appends to the rvc and erg files every 10 time steps
			if step % 10 == 0: 
				self.write_rvc_output(step, self.args.out, first_line, E_tot,'a')
				self.write_erg_output(step, self.args.out, E_k, E_b, E_nb)

if __name__ == '__main__':
	MyFirstMDSim = MDSim()
	MyFirstMDSim.run_md_sim()
