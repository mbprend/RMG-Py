################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

cimport chem # so we have 

cdef extern from "dictobject.h":
	ctypedef class __builtin__.dict [object PyDictObject]:
		pass

cdef class Vertex:

	def __init__(self):
		# for Extended Connectivity; as introduced by Morgan (1965)
		# http://dx.doi.org/10.1021/c160017a018
		self.connectivity_value_1 = 0
		self.connectivity_value_2 = 0
		self.connectivity_value_3 = 0
		pass

	cpdef bint equivalent(Vertex self, Vertex other):
		return True

cdef class Edge:

	def __init__(self):
		pass

	cpdef bint equivalent(Edge self, Edge other):
		return True

cdef class Graph(dict):
	"""
	A representation of a graph using a dictionary of dictionaries. The keys
	of the outer dictionary are the vertices, while edges are accessed via
	self[vertex1][vertex2].
	"""

	def __init__(self, vertices=None, edges=None):
		self.clear()

	def __reduce__(self):
		return (dict, (), None, None, self.iteritems())

	cpdef list vertices(Graph self):
		"""
		Return a list of the vertices in the graph.
		"""
		return self.keys()
	
	cpdef list edges(Graph self):
		"""
		Return a list of the edges in the graph.
		"""
		cdef list edgelist = list()
		cdef list pairslist = list()
		for v1 in self:
			for v2 in self[v1]:
				if (v1, v2) not in pairslist:
					edgelist.append(self[v1][v2])
					pairslist.append((v1,v2))
					pairslist.append((v2,v1))
		return edgelist
	
	cpdef addVertex(Graph self, vertex):
		"""
		Add a `vertex` to the graph. The vertex is initialized with no edges.
		"""
		self[vertex] = dict()
		return vertex

	cpdef addEdge(Graph self, vertices, edge):
		"""
		Add an `edge` to the graph as an edge connecting the two vertices
		specified in the 2-tuple `vertices`.
		"""
		v1, v2 = vertices
		self[v1][v2] = edge
		self[v2][v1] = edge
		return edge

	cpdef dict getEdges(Graph self, vertex):
		"""
		Return a list of the edges involving the specified `vertex`.
		"""
		return self[vertex]

	cpdef getEdge(Graph self, tuple vertices):
		"""
		Returns the edge connecting vertices in 2-tuple `vertices`, or None if
		no edge exists.
		"""
		v1, v2 = vertices
		return self[v1][v2] if self.hasEdge(vertices) else None

	cpdef bint hasEdge(self, tuple vertices):
		"""
		Returns :data:`True` if vertices in the 2-tuple `vertices` are connected
		by an edge, and :data:`False` if not.
		"""
		v1, v2 = vertices
		if v1 in self:
			return v2 in self[v1]
		return False

	cpdef removeVertex(Graph self, vertex1):
		"""
		Remove `vertex1` and all edges associated with it from the graph. Does
		not remove vertices that no longer have any edges as a result of this
		removal.
		"""
		for vertex2 in self:
			if vertex2 is not vertex1:
				if vertex1 in self[vertex2]:
					del self[vertex2][vertex1]
		del self[vertex1]

	cpdef removeEdge(Graph self, vertices):
		"""
		Remove the edge having vertices as specified in the 2-tuple `vertices
		from the graph. Does not remove vertices that no longer have any edges
		as a result of this removal.
		"""
		v1, v2 = vertices
		del self[v1][v2]
		del self[v2][v1]

	cpdef isIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if two graphs are isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		if len(self) != len(other): return False
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, False, False)
		return ismatch

	cpdef isSubgraphIsomorphic(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		ismatch, map21, map12 = VF2_isomorphism(self, other, map21_0, map12_0, True, False)
		return ismatch

	cpdef findSubgraphIsomorphisms(Graph self, Graph other, dict map12_0, dict map21_0):
		"""
		Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
		otherwise. Uses the VF2 algorithm of Vento and Foggia.
		"""
		return VF2_isomorphism(self, other, map21_0, map12_0, True, True)

	cpdef Graph copy(Graph self):
		"""
		Create a copy of the current graph.
		"""
		cdef Graph other = Graph()

		for vertex in self:
			other.addVertex(vertex)
		for v1 in self:
			for v2 in self[v1]:
				other[v1][v2] = self[v1][v2]
		
		return other

	cpdef Graph merge(Graph self, Graph other):
		"""
		Merge two graphs so as to store them in a single Graph object.
		"""

		# Create output graph
		cdef Graph new = Graph()

		# Add vertices to output graph
		for vertex in self:
			new.addVertex(vertex)
		for vertex in other:
			new.addVertex(vertex)

		# Add edges to output graph
		for v1 in self:
			for v2 in self[v1]:
				new[v1][v2] = self[v1][v2]
		for v1 in other:
			for v2 in other[v1]:
				new[v1][v2] = other[v1][v2]

		return new

	cpdef list split(Graph self):
		"""
		Convert a single Graph object containing two or more unconnected graphs
		into separate graphs.
		"""
		
		# Create potential output graphs
		cdef Graph new1 = self.copy()
		cdef Graph new2 = Graph()
		cdef list verticesToMove
		cdef int index

		if len(self.vertices()) == 0:
			return [new1]

		# Arbitrarily choose last atom as starting point
		verticesToMove = [ self.vertices()[-1] ]

		# Iterate until there are no more atoms to move
		index = 0
		while index < len(verticesToMove):
			for v2 in self.getEdges(verticesToMove[index]):
				if v2 not in verticesToMove:
					verticesToMove.append(v2)
			index += 1

		# If all atoms are to be moved, simply return new1
		if len(new1.vertices()) == len(verticesToMove):
			return [new1]

		# Copy to new graph
		for vertex in verticesToMove:
			new2.addVertex(vertex)
		for v1 in verticesToMove:
			for v2, edge in new1[v1].iteritems():
				new2[v1][v2] = edge

		# Remove from old graph
		for v1 in new2:
			for v2 in new2[v1]:
				if v1 in verticesToMove and v2 in verticesToMove:
					del new1[v1][v2]
		for vertex in verticesToMove:
			new1.removeVertex(vertex)

		new = [new2]
		new.extend(new1.split())
		return new

	cpdef list getSmallestSetOfSmallestRings(Graph self):
		"""
		Return a list of the smallest set of smallest rings in the graph. The
		algorithm implements was adapted from a description by Fan, Panaye,
		Doucet, and Barbu (doi: 10.1021/ci00015a002)

		B. T. Fan, A. Panaye, J. P. Doucet, and A. Barbu. "Ring Perception: A
		New Algorithm for Directly Finding the Smallest Set of Smallest Rings
		from a Connection Table." *J. Chem. Inf. Comput. Sci.* **33**, 
		p. 657-662 (1993).
		"""

		# Make a copy of the graph so we don't modify the original
		cdef Graph graph = self.copy()
		cdef bint done
		cdef list verticesToMove
		cdef list cycleList
		cdef list cycles
		cdef chem.Atom vertex, rootVertex
		cdef bint found
		cdef list cycle, graphs

		# Step 1: Remove all terminal vertices
		done = False
		while not done:
			verticesToRemove = []
			for vertex1, value in graph.iteritems():
				if len(value) == 1: verticesToRemove.append(vertex1)
			done = len(verticesToRemove) == 0
			# Remove identified vertices from graph
			for vertex in verticesToRemove:
				graph.removeVertex(vertex)

		# Step 2: Remove all other vertices that are not part of cycles
		verticesToRemove = []
		for vertex in graph:
			found = graph.isVertexInCycle(vertex)
			if not found:
				verticesToRemove.append(vertex)
		# Remove identified vertices from graph
		for vertex in verticesToRemove:
			graph.removeVertex(vertex)
			
		### also need to remove EDGES that are not in ring

		# Step 3: Split graph into remaining subgraphs
		graphs = graph.split()

		# Step 4: Find ring sets in each subgraph
		cycleList = []
		for graph in graphs:

			while len(graph) > 0:

				# Choose root vertex as vertex with smallest number of edges
				rootVertex = None
				for vertex in graph:
					if rootVertex is None:
						rootVertex = vertex
					elif len(graph[vertex]) < len(graph[rootVertex]):
						rootVertex = vertex

				# Get all cycles involving the root vertex
				cycles = graph.getAllCycles(rootVertex)
				if len(cycles) == 0:
					# this vertex is no longer in a ring.
					# remove all its edges
					neighbours = graph[rootVertex].keys()[:]
					for vertex2 in neighbours:
						graph.removeEdge((rootVertex, vertex2))
					# then remove it
					graph.removeVertex(rootVertex)
					#print("Removed vertex that's no longer in ring")
					continue # (pick a new root Vertex)
#					raise Exception('Did not find expected cycle!')

				# Keep the smallest of the cycles found above
				cycle = cycles[0]
				for c in cycles[1:]:
					if len(c) < len(cycle):
						cycle = c
				cycleList.append(cycle)

				# Remove from the graph all vertices in the cycle that have only two edges
				verticesToRemove = []
				for vertex in cycle:
					if len(graph[vertex]) <= 2:
						verticesToRemove.append(vertex)
				if len(verticesToRemove) == 0:
					# there are no vertices in this cycle that with only two edges
					
					# Remove edge between root vertex and any one vertex it is connected to
					graph.removeEdge((rootVertex, graph[rootVertex].keys()[0]))
				else:
					for vertex in verticesToRemove:
						graph.removeVertex(vertex)
						
		return cycleList

	cpdef bint isVertexInCycle(Graph self, chem.Atom vertex):
		""" 
		Is `vertex` in a cycle?
		Returns :data:`True` if it is in a cycle, else :data:`False`.
		"""
		cdef list chain
		chain = [vertex]
		return self.__isChainInCycle(chain)

	cpdef bint __isChainInCycle(Graph self, list chain):
		""" 
		Is the `chain` in a cycle? 
		Returns True/False.
		Recursively calls itself
		"""
		# Note that this function no longer returns the cycle; just True/False
		cdef chem.Atom vertex2
		cdef chem.Bond edge
		cdef bint found
		
		for vertex2, edge in self[chain[-1]].iteritems():
			if vertex2 is chain[0] and len(chain) > 2:
				return True
			elif vertex2 not in chain:
				# make the chain a little longer and explore again
				chain.append(vertex2)
				found = self.__isChainInCycle(chain)
				if found: return True
				# didn't find a cycle down this path (-vertex2),
				# so remove the vertex from the chain
				chain.remove(vertex2)
		return False

	cpdef list getAllCycles(Graph self, chem.Atom startingVertex):
		"""
		Given a starting vertex, returns a list of all the cycles containing 
		that vertex.
		"""
		cdef list chain, cycleList
		
		cycleList=list()
		chain = [startingVertex]
		
		chainLabels=range(len(self.keys()))
		#print "Starting at %s in graph: %s"%(self.keys().index(startingVertex),chainLabels)
		
		cycleList = self.__exploreCyclesRecursively(chain, cycleList)
		return cycleList
		
	cpdef list __exploreCyclesRecursively(Graph self, list chain, list cycleList):
		"""
		Finds cycles by spidering through a graph.
		Give it a chain of atoms that are connected, `chain`,
		and a list of cycles found so far `cycleList`.
		If `chain` is a cycle, it is appended to `cycleList`.
		Then chain is expanded by one atom (in each available direction)
		and the function is called again. This recursively spiders outwards
		from the starting chain, finding all the cycles.
		"""
		# unless we derive Atom and Bond from Vertex and Edge we can't cdef these:
		#cdef Vertex vertex2
		#cdef Edge edge
		cdef chem.Atom vertex2
		cdef chem.Bond edge
		
		chainLabels=[self.keys().index(v) for v in chain]
		#print "found %d so far. Chain=%s"%(len(cycleList),chainLabels)
		
		for vertex2, edge in self[chain[-1]].iteritems():
			# vertex2 will loop through each of the atoms 
			# that are bonded to the last atom in the chain.
			if vertex2 is chain[0] and len(chain) > 2:
				# it is the first atom in the chain - so the chain IS a cycle!
				cycleList.append(chain[:])
			elif vertex2 not in chain:
				# make the chain a little longer and explore again
				chain.append(vertex2)
				cycleList = self.__exploreCyclesRecursively(chain, cycleList)
				# any cycles down this path (-vertex2) have now been found,
				# so remove the vertex from the chain
				chain.pop(-1)
		return cycleList

	cpdef set_connectivity_values(Graph self):
		"""
		Sets the Extended Connectivity values as introduced by Morgan (1965)
		http://dx.doi.org/10.1021/c160017a018
		
		First CV1 is the number of neighbours  (Morgan proposed non-Hydrogen neighbours)
		CV2 is the sum of neighbouring CV1 values
		CV3 is the sum of neighbouring CV2 values
		"""
		
		cdef int count
		cdef chem.Atom vert1, vert2
		
		for vert1 in self:
			vert1.connectivity_value_1 = len(<dict>self[vert1])
			
		for vert1 in self:
			count = 0
			for vert2 in <dict>self[vert1]:
				count += vert2.connectivity_value_1
			vert1.connectivity_value_2 = count
			
		for vert1 in self:
			count = 0
			for vert2 in <dict>self[vert1]:
				count += vert2.connectivity_value_2
			vert1.connectivity_value_3 = count
		
		
################################################################################

cpdef VF2_isomorphism(Graph graph1, Graph graph2, dict map12, dict map21, \
	bint subgraph=False, bint findAll=False):
	"""
	Returns :data:`True` if two :class:`Graph`s are isomorphic and :data:`False`
	otherwise. Uses the VF2 algorithm of Vento and Foggia. If `subgraph` is
	:data:`True` then graph2 is checked for being a potential subgraph of graph1. 
	`findAll` isused to specify whether all isomorphisms should be returned, 
	or only the first.
	"""

	cdef list map12List = list(), map21List = list()
	cdef bint ismatch
	cdef dict terminals1, terminals2

	graph1.set_connectivity_values() # could probably run this less often elsewhere
	graph2.set_connectivity_values() # as the values don't change often
		
	terminals1 = __VF2_terminals(graph1, map21)
	terminals2 = __VF2_terminals(graph2, map12)

	ismatch = __VF2_match(graph1, graph2, map21, map12, \
		 terminals1, terminals2, subgraph, findAll, map21List, map12List)

	if findAll:
		return len(map21List) > 0, map21List, map12List
	else:
		return ismatch, map12, map21

cdef bint __VF2_feasible(Graph graph1, Graph graph2, chem.Atom vertex1, chem.Atom vertex2, \
	dict map21, dict map12, dict terminals1, dict terminals2, bint subgraph):
	"""
	Returns :data:`True` if two vertices `vertex1` and `vertex2` from graphs
	`graph1` and `graph2`, respectively, are feasible matches. `mapping21` and
	`mapping12` are the current state of the mapping from `graph1` to `graph2`
	and vice versa, respectively. `terminals1` and `terminals2` are lists of
	the vertices that are directly connected to the already-mapped vertices.
	`subgraph` is :data:`True` if graph2 is to be treated as a potential
	subgraph of graph1.

	Uses the VF2 algorithm of Vento and Foggia. The feasibility is assessed
	through a series of semantic and structural checks. Only the combination
	of the semantic checks and the level 0 structural check are both
	necessary and sufficient to ensure feasibility. (This does *not* mean that
	vertex1 and vertex2 are always a match, although the level 1 and level 2
	checks preemptively eliminate a number of false positives.)
	"""

	cdef chem.Bond edge1, edge2
	cdef chem.Atom vert1, vert2

	# Richard's Connectivity Value check
	# not sure where this is best done. Is it more specific or more general?
#	if subgraph:
#		if vertex1.connectivity_value_1 < vertex2.connectivity_value_1: return False
#		if vertex1.connectivity_value_2 < vertex2.connectivity_value_2: return False
#		if vertex1.connectivity_value_3 < vertex2.connectivity_value_3: return False
#	else:
	if not subgraph:
		if vertex1.connectivity_value_1 != vertex2.connectivity_value_1: return False
		if vertex1.connectivity_value_2 != vertex2.connectivity_value_2: return False
		if vertex1.connectivity_value_3 != vertex2.connectivity_value_3: return False
		
	# Semantic check #1: vertex1 and vertex2 must be equivalent
	if not vertex1.equivalent(vertex2):
		return False
	
	
	# Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
	# already mapped should be connected by equivalent edges
	cdef dict edges1 = <dict>graph1[vertex1]
	cdef dict edges2 = <dict>graph2[vertex2]
		
	for vert1 in edges1:
	# for vert1, edge1 in edges1.iteritems(): # if you uncomment this..**
		if vert1 in map21:
			vert2 = map21[vert1]
			if not vert2 in edges2:
				return False
			edge1 = edges1[vert1] # **..then remove this
			edge2 = edges2[vert2]
			if not edge1.equivalent(edge2):
				return False
	
	# Count number of terminals adjacent to vertex1 and vertex2
	cdef int term1Count = 0, term2Count = 0, neither1Count = 0, neither2Count = 0
	
	for vert1 in edges1:
		if vert1 in terminals1:
			term1Count += 1
		elif vert1 not in map21:
			neither1Count += 1
	for vert2 in edges2:
		if vert2 in terminals2:
			term2Count += 1
		elif vert2 not in map12:
			neither2Count += 1

	# Level 2 look-ahead: the number of adjacent vertices of vertex1 and
	# vertex2 that are non-terminals must be equal
	if subgraph:
		if neither1Count < neither2Count:
			return False
	else:
		if neither1Count != neither2Count:
			return False

	# Level 1 look-ahead: the number of adjacent vertices of vertex1 and
	# vertex2 that are terminals must be equal
	if subgraph:
		if term1Count < term2Count:
			return False
	else:
		if term1Count != term2Count:
			return False

	# Level 0 look-ahead: all adjacent vertices of vertex1 already in the
	# mapping must map to adjacent vertices of vertex2
	for vert1 in edges1:
		if vert1 in map21:
			vert2 = map21[vert1]
			if vert2 not in edges2:
				return False
	return True

cdef bint __VF2_match(Graph graph1, Graph graph2, dict map21, dict map12, \
	dict terminals1, dict terminals2, bint subgraph, bint findAll, \
	list map21List, list map12List):
	"""
	A recursive function used to explore two graphs `graph1` and `graph2` for
	isomorphism by attempting to map them to one another. `mapping21` and
	`mapping12` are the current state of the mapping from `graph1` to `graph2`
	and vice versa, respectively. `terminals1` and `terminals2` are lists of
	the vertices that are directly connected to the already-mapped vertices.
	`subgraph` is :data:`True` if graph2 is to be treated as a potential
	subgraph of graph1.

	Uses the VF2 algorithm of Vento and Foggia, which is O(N) in spatial complexity
	and O(N**2) (best-case) to O(N! * N) (worst-case) in temporal complexity.
	"""

	cdef dict new_terminals1, new_terminals2
	cdef chem.Atom vertex1, vertex2
	cdef bint ismatch 
	
	# Done if we have mapped to all vertices in graph2
	if len(map12) >= len(graph2) or len(map21) >= len(graph1):
		return True
	
	#cdef dict terminals1, terminals2
	#terminals1 = __VF2_terminals(graph1, map21)
	#terminals2 = __VF2_terminals(graph2, map12)
	
	# Create list of pairs of candidates for inclusion in mapping
	cdef list pairs = __VF2_pairs(graph1, graph2, terminals1, terminals2)
	
	for vertex1, vertex2 in pairs:
		# propose a pairing
		if __VF2_feasible(graph1, graph2, vertex1, vertex2, map21, map12, \
				terminals1, terminals2, subgraph):
			# Update mapping accordingly
			map21[vertex1] = vertex2
			map12[vertex2] = vertex1
			
			# update terminals
			new_terminals1 = __VF2_new_terminals(graph1, map21, terminals1, vertex1)
			new_terminals2 = __VF2_new_terminals(graph2, map12, terminals2, vertex2)

			# Recurse
			ismatch = __VF2_match(graph1, graph2, \
				map21, map12, new_terminals1, new_terminals2, subgraph, findAll, \
				map21List, map12List)
			if ismatch:
				if findAll:
					map21List.append(map21.copy())
					map12List.append(map12.copy())
				else:
					return True
			# Undo proposed match
			del map21[vertex1]
			del map12[vertex2]
			# changes to 'new_terminals' will be discarded and 'terminals' is unchanged

	return False

cpdef int __atom_sort_value(chem.Atom atom):
	"""
	Used to sort atoms prior to poposing candidate pairs in :method:`__VF2_pairs` 
	The lowest (or most negative) values will be first in the list when you sort, 
	so should be given to the atoms you want to explore first. 
	For now, that is (roughly speaking) the most connected atoms. This definitely helps for large graphs
	but bizarrely the opposite ordering seems to help small graphs. Not sure about subggraphs...
	"""
	return -100 * atom.connectivity_value_1 - 10 * atom.connectivity_value_2  - atom.connectivity_value_3
	
cdef list __VF2_pairs(Graph graph1, Graph graph2, dict terminals1, dict terminals2):
	"""
	Create a list of pairs of candidates for inclusion in the VF2 mapping. If
	there are a nonzero number of terminals in each graph, the candidates are
	selected to be one terminal from the first graph and all terminals from the
	second graph. If there are no terminals, the candidates are	selected to be
	one vertex from the first graph and all vertices from the second graph.
	"""

	cdef list pairs = list()
	cdef chem.Atom vertex1, vertex2, terminal1, terminal2
	cdef list list_to_sort
	# Construct list from terminals if possible
	if len(terminals1) > 0 and len(terminals2) > 0:
		list_to_sort = terminals2.keys()
		list_to_sort.sort(key=__atom_sort_value)
		terminal2 = list_to_sort[0]
		list_to_sort = terminals1.keys()
		list_to_sort.sort(key=__atom_sort_value)
		
		for terminal1 in list_to_sort:
			pairs.append([terminal1, terminal2])
	# Otherwise construct list from all remaining vertices
	else:
		list_to_sort = graph2.keys()
		list_to_sort = graph2.keys()
		list_to_sort.sort(key=__atom_sort_value)
		vertex2 = list_to_sort[0]
		list_to_sort = graph1.keys()
		list_to_sort.sort(key=__atom_sort_value)		
		for vertex1 in list_to_sort:
			pairs.append([vertex1, vertex2])
	
	return pairs

cdef dict __VF2_terminals(Graph graph, dict mapping):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	generate a list of terminals, vertices that are directly connected to
	vertices that have already been mapped.
	"""

	cdef dict terminals = dict() # why won't {} work?
	cdef chem.Atom vertex, vert
	cdef chem.Bond edge
	
	for vertex in mapping:
		for vert in <dict>graph[vertex]:
			if vert not in mapping:
				terminals[vert] = True
	return terminals

cdef dict __VF2_new_terminals(Graph graph, dict mapping, dict old_terminals, new_vertex):
	"""
	For a given graph `graph` and associated partial mapping `mapping`,
	UPDATES a list of terminals, vertices that are directly connected to
	vertices that have already been mapped. You have to pass it the previous 
	list of terminals `old_terminals` and the vertex `vertex` that has been added 
	to the mapping. Returns a new copy of the terminals.
	"""
	
	cdef dict terminals = dict() # why won't {} work?
	
	# copy the old terminals, leaving out the new_vertex
	for vertex in old_terminals:
		if not vertex is new_vertex: 
			terminals[vertex] = True
	
	# add the terminals of new_vertex
	for vertex in <dict>graph[new_vertex]:
		if vertex not in mapping:
			terminals[vertex] = True
	return terminals

################################################################################

