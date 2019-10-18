
"""

Usage: python align.py input_file output_file

"""

## Richard Hojel
## Smith-Waterman and Needleman-Wunsch algorithms for sequence alignment.

import sys

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
	"""
	Checks if two floating point numbers are equivalent.
	"""
	epsilon = 10**(-6) 
	return (abs(a - b) < epsilon)
	

#### ------- CLASSES ------- ####

class MatchMatrix(object):
	"""
	Match matrix class stores the scores of matches in a data structure
	"""
	def __init__(self):

		self.match_matrix = {("",""):0}

	def set_score(self, a, b, score):
		"""
		Updates or adds a score for a specified match

		Inputs:
		   a = the character from sequence A
		   b = the character from sequence B
		   score = the score to set it for
		"""
		self.match_matrix[(a,b)] = score

	def get_score(self, a, b):
		"""
		Returns the score for a particular match, where a is the
		character from sequence a and b is from sequence b.

		Inputs:
		   a = the character from sequence A
		   b = the character from sequence B
		Returns:
		   the score of that match
		"""
		return self.match_matrix[(a,b)]

class ScoreMatrix(object):
	"""
	Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
	ScoreEntries that are updated during alignment and used to output the maximum alignment.
	"""

	def __init__(self, name, nrow, ncol):
		self.name = name # identifier for the score matrix - Ix, Iy, or M
		self.nrow = nrow
		self.ncol = ncol
		self.score_matrix = {(0,0):{'score':0., 'pointers':[]}}

		# Sets all end gap penalties to zero in the score matrix and adds a pointer that 
		# points to the next end gap in the matrix
		for i in range(1,nrow+1):
			self.score_matrix[(i,0)] = {'score':0., 'pointers':[(name,i-1,0)]}
		for j in range(1,ncol+1):
			self.score_matrix[(0,j)] = {'score':0., 'pointers':[(name,0,j-1)]}

	def get_score(self, row, col):
		### FILL IN ###
		return self.score_matrix[(row,col)]['score']
		
	def set_score(self, row, col, score):    
		### FILL IN ###
		self.score_matrix[(row,col)] = {'score': score, 'pointers':[]}

	def get_pointers(self, row, col):
		"""
		Returns the indices of the entries that are pointed to
		This should be formatted as a list of tuples:
		 ex. [(1,1), (1,0)]
		"""
		return self.score_matrix[(row,col)]['pointers']

	def set_pointers(self, row, col, pname, prow, pcol): ### FILL IN - this needs additional arguments ###
		### FILL IN ###
		self.score_matrix[(row,col)]['pointers'].append((pname,prow,pcol))

	def print_scores(self):
		"""
		Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

		Example:
		M=
			0.0, 0.0, 0.0, 0.0, 0.0
			0.0, 1.0, 0.0, 0.0, 0.0
			0.0, 1.0, 1.0, 1.0, 1.0
			0.0, 0.0, 1.0, 1.0, 1.0
			0.0, 0.0, 2.0, 2.0, 1.0
			0.0, 0.0, 1.0, 2.0, 3.0

		"""
		print(self.name, "=")
		for i in range(0, self.nrow +1):
			print("\n	", end = '')
			for j in range(0, self.ncol +1):
				print("%.1f" % self.get_score(i,j), end = '')
				if j < self.ncol:
					print(", ", end = '')
		print('\n', end = '')

	def print_pointers(self):
		"""
		Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
		"""
		print(self.name, "=")
		for i in range(0, self.nrow +1):
			print("\n 	", end = '')
			for j in range(0, self.ncol +1):
				print(self.get_pointers(i,j), end = '')
				if j < self.ncol:
					print(", ", end = '')
		print('\n', end = '')

class AlignmentParameters(object):
	
	def __init__(self):
		# default values for variables that are filled in by reading
		# the input alignment file
		self.seq_a = ""
		self.seq_b = ""
		self.global_alignment = False 
		self.dx = 0
		self.ex = 0
		self.dy = 0
		self.ey = 0
		self.alphabet_a = "" 
		self.alphabet_b = ""
		self.len_alphabet_a = 0
		self.len_alphabet_b = 0
		self.match_matrix = MatchMatrix()

	def load_params_from_file(self, input_file): 
		"""
		Reads the parameters from an input file and stores in the object

		Input:
		   input_file = specially formatted alignment input file
		"""
		with open(input_file, 'r') as input_f:

			self.seq_a = input_f.readline().rstrip()
			self.seq_b = input_f.readline().rstrip()

			self.global_alignment = (int(input_f.readline()) == 0)

			penalty = input_f.readline().split()
			self.dx = float(penalty[0])
			self.ex = float(penalty[1])
			self.dy = float(penalty[2])
			self.ey = float(penalty[3])
			
			self.len_alphabet_a = int(input_f.readline())
			self.alphabet_a = input_f.readline().rstrip()

			self.len_alphabet_b = int(input_f.readline())
			self.alphabet_b = input_f.readline().rstrip()

			for matrix_line in input_f:
				matrix_values = matrix_line.split()
				if not matrix_values:
					break
				self.match_matrix.set_score(matrix_values[2],matrix_values[3],float(matrix_values[4])) 	
	
	def no_negatives_in_match_matrix(self):
		## FILL IN
		if not self.global_alignment:
			for a in self.alphabet_a:
				for b in self.alphabet_b:
					if (self.match_matrix.get_score(a,b) < 0):
						return False
			return True
		return False

class AlignmentAlgorithm(object):

	def __init__(self, nrow, ncol):

		### FILL IN - note be careful about how you initialize these! ###
		self.m_matrix = ScoreMatrix("M",nrow,ncol)
		self.ix_matrix = ScoreMatrix("Ix",nrow,ncol)
		self.iy_matrix = ScoreMatrix("Iy",nrow,ncol)

	def populate_score_matrices(self, align_params):
		"""
		Method to populate the score matrices.
		Should call update(i,j) for each entry in the score matrices

		Input:
		   align_params = an AlignmentParameters object with the input
		"""
		for i in range(1,len(align_params.seq_a)+1):
			for j in range(1,len(align_params.seq_b)+1):
				self.update(i, j, align_params)

	def update(self, row, col, align_params):
		"""
		Method to update the matrices at a given row and column index.

		Input:
		   row = the row index to update
		   col = the column index to update
		   align_params = an AlignmentParameters object with the input
		"""
		self.update_m(row, col, align_params)
		self.update_ix(row, col, align_params)
		self.update_iy(row, col, align_params)

	def update_m(self, row, col, align_params):
		### FILL IN ###
		a = align_params.seq_a[row-1]
		b = align_params.seq_b[col-1]

		m = self.m_matrix.get_score(row-1,col-1) + align_params.match_matrix.get_score(a,b)
		ix = self.ix_matrix.get_score(row-1,col-1) + align_params.match_matrix.get_score(a,b)
		iy = self.iy_matrix.get_score(row-1,col-1) + align_params.match_matrix.get_score(a,b)

		if align_params.global_alignment:
			self.m_matrix.set_score(row, col, max(m,ix,iy))
		else:
			self.m_matrix.set_score(row, col, max(m,ix,iy,0)) 
		
		if fuzzy_equals(m,max(m,ix,iy)):
			self.m_matrix.set_pointers(row,col,"M",row-1,col-1)
		if fuzzy_equals(ix,max(m,ix,iy)):
			self.m_matrix.set_pointers(row,col,"Ix",row-1,col-1)
		if fuzzy_equals(iy,max(m,ix,iy)):
			self.m_matrix.set_pointers(row,col,"Iy",row-1,col-1)

	def update_ix(self, row, col, align_params):
		### FILL IN ###
		m = self.m_matrix.get_score(row-1,col) - align_params.dy
		ix = self.ix_matrix.get_score(row-1,col) - align_params.ey
		
		self.ix_matrix.set_score(row, col, max(m,ix))
		if fuzzy_equals(m,max(m,ix)):
			self.ix_matrix.set_pointers(row,col,"M",row-1,col)
		if fuzzy_equals(ix,max(m,ix)):
			self.ix_matrix.set_pointers(row,col,"Ix",row-1,col)

	def update_iy(self, row, col, align_params):
		### FILL IN ###
		m = self.m_matrix.get_score(row,col-1) - align_params.dx
		iy = self.iy_matrix.get_score(row,col-1) - align_params.ex

		self.iy_matrix.set_score(row, col, max(m,iy))
		if fuzzy_equals(m,max(m,iy)): 
			self.iy_matrix.set_pointers(row,col,"M",row,col-1)
		if fuzzy_equals(iy,max(m,iy)):
			self.iy_matrix.set_pointers(row,col,"Iy",row,col-1)	

#### ------- MAIN METHODS ------- ####

def find_traceback_start(align_alg, align_params):
	"""
	Finds the location to start the traceback.
	Think carefully about how to set this up for local 

	Inputs:
		align_alg = an AlignmentAlgorithm object with populated score matrices
		align_params = an AlignmentParameters object with the input

	Returns:
		(max_val, max_loc) where max_val is the best score
		max_loc is a list [] containing tuples with the (i,j) location(s) to start the traceback
		 (ex. [(1,2), (3,4)])
	"""
	max_val = 0
	max_loc = []

	if not align_params.global_alignment:
		for i in range(0, len(align_params.seq_a)+1):
			for j in range(0, len(align_params.seq_b)+1):
				score = align_alg.m_matrix.get_score(i,j)
				if score > max_val:
					max_val = score
					max_loc = []
					max_loc.append(('M',i,j))
				elif fuzzy_equals(score, max_val):
					max_loc.append(('M',i,j))
	else:
		for i in range(0, len(align_params.seq_a)+1):
			score = align_alg.m_matrix.get_score(i,len(align_params.seq_b))
			if score > max_val:
				max_val = score
				max_loc = []
				max_loc.append(('M',i,len(align_params.seq_b)))
		for j in range(0, len(align_params.seq_b)+1):
			score = align_alg.m_matrix.get_score(len(align_params.seq_a),j)
			if score > max_val:
				max_val = score
				max_loc = []
				max_loc.append(('M',len(align_params.seq_a),j))

	return max_val, max_loc

def traceback(align_alg, align_params, output_file): ### FILL IN additional arguments ###
	"""
	Performs a traceback.
	Hint: include a way to printing the traceback path. This will be helpful for debugging!
	   ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)

	Input:
		align_alg = an AlignmentAlgorithm object with populated score matrices
		align_params = an AlignmentParameters object with the input
	"""
	final_score, max_loc = find_traceback_start(align_alg,align_params)

	with open(output_file, 'w') as output_f:
		output_f.write(str(round(final_score,1)))
		output_f.write('\n')

	## Implementation of Depth First Search algorithm to traceback 
	path_stack = []
	path_trace = {}

	for loc in max_loc:
		path_trace[loc] = None
		path_stack.append(loc)

	while not path_stack == []:
		curr_loc = path_stack.pop()

		if curr_loc[0] == 'M':
			next_loc = align_alg.m_matrix.get_pointers(curr_loc[1],curr_loc[2])
		elif curr_loc[0] == 'Ix':
			next_loc = align_alg.ix_matrix.get_pointers(curr_loc[1],curr_loc[2])
		else:
			next_loc = align_alg.iy_matrix.get_pointers(curr_loc[1],curr_loc[2])

		for loc in next_loc:
			if loc[0] == 'M' and ((loc[1] == 0 or loc[2] == 0) or (not align_params.global_alignment 
				and fuzzy_equals(align_alg.m_matrix.get_score(loc[1],loc[2]),0))):
				
				path_trace[loc] = curr_loc
				output_alignment(loc, path_trace, align_params, output_file)

				## Unit Tests to print out the traceback and print out the alignments 
				# print_traceback(loc, path_trace)	# Uncomment to use test
				# print_alignment(loc, path_trace, align_params)	#Uncomment to use tests
			else:
				path_trace[loc] = curr_loc
				path_stack.append(loc)

def print_traceback(loc, path_trace):
	## FILL IN
	if loc is None:
		print()
		return

	print(loc, "<-", end = '')
	print_traceback(path_trace[loc], path_trace) 


def print_alignment(loc_start, path_trace, align_params):
	## FILL IN
	seq_a = True
	print(output_seq(loc_start, path_trace, align_params, seq_a))

	seq_a = False
	print(output_seq(loc_start, path_trace, align_params, seq_a),'\n')

def output_alignment(loc_start, path_trace, align_params, output_file):
	## FILL IN
	with open(output_file, 'a') as output_f:

		output_f.write('\n')
		seq_a = True
		output_f.write(output_seq(loc_start, path_trace, align_params, seq_a))
		
		output_f.write('\n')
		seq_a = False
		output_f.write(output_seq(loc_start, path_trace, align_params, seq_a))
		
		output_f.write('\n')

def output_seq(loc_start, path_trace, align_params, seq_a):
	## FILL IN
	if seq_a:
		seq_letter = 1
		gap_matrix = 'Ix'
		org_seq = align_params.seq_a
	else:
		seq_letter = 2
		gap_matrix = 'Iy'
		org_seq = align_params.seq_b

	loc = path_trace[loc_start]
	seq = ''

	while loc is not None:
		if loc[0] == 'M' or loc[0] == gap_matrix:
			seq += org_seq[loc[seq_letter]-1]
		else:
			seq += '_'
		loc = path_trace[loc]

	return seq

def main():

	# check that the file is being properly used
	if (len(sys.argv) !=3):
		print("Please specify an input file and an output file as args.")
		return
		
	# input variables
	input_file = sys.argv[1]
	output_file = sys.argv[2]

	# read in the alignment data    
	align_params = AlignmentParameters()
	align_params.load_params_from_file(input_file)

	# terminates the program if conducting a local alignment without negative values in the match matrix
	if align_params.no_negatives_in_match_matrix():
		print('Cannot conduct local alignment if the match matrix contains no negatives. Please select global alignment or input a different match matrix.')
	else:
		# populate the score matrices based on the input parameters
		align_object = AlignmentAlgorithm(len(align_params.seq_a), len(align_params.seq_b))
		align_object.populate_score_matrices(align_params)

		# perform a traceback and write the output to an output file
		traceback(align_object,align_params,output_file)

if __name__=="__main__":
	main()
