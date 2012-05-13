#!/usr/bin/python
import sys
import math
import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages


"""
The module determines and graphically produces Markov equivalence classes for a given
relational schema.  DAPER models are represented by networkx.DiGraph objects, where 
each vertex is a (entity, variable) tuple.

For a bipartite graph with entity types A, B and variables A.X, A.Z, B.Y, usage would
typically look like this: 

varX = ("A", "X")
varZ = ("A", "Z")
varY = ("B", "Y")
variables = [varX, varZ, varY]
graphGen = PossibleGraphicalGenerator(variables)
graphGen.edgeProhibit(varX, varZ, bidir=True)
graphGen.makeLatent(varZ)

possibleGraphs = graphGen.generate()
equivClasses = getGraphicalEquivClasses(possibleGraphs, True)
drawGraphicalEquivClasses(equivClasses, outfileName, rowsPerPage=20, colsPerPage=4, pageWidth=11.5, pageHeight=64)
"""




class PossibleGraphicalGenerator():
	"""Factory class for generating all possible valid DAGs for a given schema and set of constraints."""

	def __init__(self, entVars, indexVars=None):
		"""Constructor.  Accepts schema information as a list of (entity, variable) tuples."""
		self.entVars = entVars
		self.dagsOnly = True
		self.latents = set()
		self.determines = {}
		self.indexSet = None if indexVars is None else set(indexVars)
		self.edgePossibleSettings = {}
		for i in self.entVars:
			for j in self.entVars:
				if (i < j):
					edge = (i, j)
					self.edgePossibleSettings[edge] = [-1, 0, 1]
		sys.stderr.write("init:\n" + self.dumpEdgePossibleSettings())
		
	def edgeFix(self, s, t):
		"""Add a constraint that a given edge must exist in all DAGs generated."""
		edge = (s, t)
		if (edge in self.edgePossibleSettings):
			self.edgePossibleSettings[edge] = [1]
		else:
			edgeRev = (edge[1], edge[0])
			self.edgePossibleSettings[edgeRev] = [-1]
			
	def edgeProhibit(self, s, t, bidir=False):
		"""Add a constraint that prohibits a given edge in all DAGs generated."""
		edge = (s, t)
		if (edge in self.edgePossibleSettings):
			self.edgePossibleSettings[edge].remove(1)
			if (bidir):
				self.edgePossibleSettings[edge].remove(-1)
		else:
			edgeRev = (t, s)
			self.edgePossibleSettings[edgeRev].remove(-1)
			if (bidir):
				self.edgePossibleSettings[edgeRev].remove(1)

	def edgeProhibitAll(self, s, bidir=False):
		"""Add a constraint that prohibits all edges associated with a vertex."""
		sys.stderr.write("prohibit all before:\n" + self.dumpEdgePossibleSettings())
		edgesPoss = self.edgePossibleSettings
		for t in self.entVars:
			edge = (s, t)
			if (edge in edgesPoss):
				if (1 in edgesPoss[edge]):
					edgesPoss[edge].remove(1)
				if (bidir):
					if (-1 in edgesPoss[edge]):
						edgesPoss[edge].remove(-1)
			edge = (t, s)
			if (edge in edgesPoss):
				if (-1 in edgesPoss[edge]):
					edgesPoss[edge].remove(-1)
				if (bidir):
					if (1 in edgesPoss[edge]):
						edgesPoss[edge].remove(1)

	def makeLatent(self, n):
		"""Designate an edge as latent, which prevents it from being utilized in any
		conditional independence tests used for partitioning DAGs into Markov equvalence
		classes."""
		self.latents.add(n)
		
	def makeDetermines(self, target, sources):
		"""Designate a variable as determinstic, thus affecting the outcome of different
		conditional independence tests."""
		self.determines[target] = set(sources)

	def generate(self):
		"""Workhorse factory method for producing all valid DAGs for this schema and set
		of constraints."""
		graphs = []
		edgeCombos = factorialDict(self.edgePossibleSettings)
		edgesPossible = sorted(self.edgePossibleSettings.keys())
		for edgeCombo in edgeCombos: # edgeCombo is a dict (s, t) -> -1
			graph = nx.DiGraph()
			graphSig = ""
			for i, ev in enumerate(self.entVars):
				lat = True if (ev in self.latents) else False
				det = self.determines.get(ev, None)
				graph.add_node(ev, latent=lat, determines=det, order=i) 
			for edge in edgesPossible: # sorted for indexing
				setting = edgeCombo[edge]
				if (setting == 1):
					graph.add_edge(edge[0], edge[1])
				elif (setting == -1):
					graph.add_edge(edge[1], edge[0])		
				if (self.indexSet is None) or ((edge[0] in self.indexSet) and (edge[1] in self.indexSet)):
					graphSig += str([1, 0, 2][setting + 1]) # translate from [-1, 0, 1] to [1, 0, 2]
			graph.graph['index'] = int(graphSig, 3)
			if (not self.dagsOnly) or (nx.is_directed_acyclic_graph(graph)):
				graphs.append(graph)
		if (len(graphs) < len(edgeCombos)):
			sys.stderr.write("eliminated %d cyclic graphs\n" % (len(edgeCombos) - len(graphs)))
		return sorted(graphs, key=lambda x: x.graph['index'])
	
	def dumpEdgePossibleSettings(self):
		"""Debugging."""
		ret = ""
		for edge, settings in self.edgePossibleSettings.items():
			s, t = edge
			ret += "%s.%s --- %s.%s\t%s\n" % (s[0], s[1], t[0], t[1], str(settings))
		return ret


#
# Convenience methods for producing common cases.
#

def possibleOneMany(latentZ=True, indexId=True):
	varZ = ("A", "Z")
	varX = ("A", "X")
	varY = ("B", "Y")
	index = None if indexId else [varX, varY, varZ]
	graphGen = PossibleGraphicalGenerator([varX, varZ, varY], index)
	graphGen.edgeProhibit(varX, varZ, bidir=True)
	if (latentZ):
		graphGen.makeLatent(varZ)
	return graphGen.generate()

def possibleOneManyId(latentZ=True, latentId=False, indexId=True):
	varId = ("A", "ID")
	varZ = ("A", "Z")
	varX = ("A", "X")
	varY = ("B", "Y")
	index = None if indexId else [varX, varY, varZ]
	graphGen = PossibleGraphicalGenerator([varId, varX, varZ, varY], index)
	graphGen.edgeFix(varId, varX)
	graphGen.edgeFix(varId, varZ)
	graphGen.edgeProhibit(varX, varZ, bidir=True)
	graphGen.edgeProhibit(varId, varY, bidir=True)
	graphGen.makeDetermines(varX, [varId])
	graphGen.makeDetermines(varZ, [varId])
	if (latentZ):
		graphGen.makeLatent(varZ)
	if (latentId):
		graphGen.makeLatent(varId)
	return graphGen.generate()

def possibleOneManyDegDisp(undirDegFy=False, latentDeg=False):	
	X = ("A", "X")
	fY = ("A", "fY")
	deg = ("A", "deg")
	Y = ("B", "Y")
	H = ("A", "H")
	entvars = [X, deg, H, fY, Y] if undirDegFy else [X, deg, fY, Y]
	graphGen = PossibleGraphicalGenerator(entvars)
	graphGen.edgeFix(Y, fY)
	graphGen.edgeProhibit(X, fY, bidir=True)
	graphGen.edgeProhibit(Y, deg, bidir=True)
	graphGen.makeLatent(Y)
	if (latentDeg):
		graphGen.makeLatent(deg)		
	if (undirDegFy):
		graphGen.edgeProhibit(deg, H, bidir=False)
		graphGen.edgeProhibit(fY, H, bidir=False)
		graphGen.edgeProhibit(H, X, bidir=True)
		graphGen.edgeProhibit(H, Y, bidir=True)
		graphGen.edgeProhibit(deg, fY, bidir=True)
		graphGen.makeLatent(H)
	return graphGen.generate()

def possibleBlocking(linksId=True, latentZ=True, indexId=True):
	varId = ("A", "id")
	varZ = ("A", "Z")
	varX = ("B", "X")
	varY = ("B", "Y")
	index = None if indexId else [varX, varY, varZ]
	graphGen = PossibleGraphicalGenerator([varId, varZ, varX, varY], index)
	if (latentZ):	
		graphGen.makeLatent(varZ)
	if (linksId):
		graphGen.edgeFix(varId, varZ)
	else:
		graphGen.edgeProhibit(varId, varZ, bidir=True)
	graphGen.edgeProhibit(varId, varX, bidir=True)
	graphGen.edgeProhibit(varId, varY, bidir=True)
	graphGen.makeDetermines(varZ, [varId])
	return graphGen.generate()

def possibleBlockingZH(linksId=True, indexId=True):
	varId = ("A", "id")
	varZ = ("A", "Z")
	varH = ("A", "H")
	varX = ("B", "X")
	varY = ("B", "Y")
	index = None if indexId else [varX, varY, varZ, varH]
	graphGen = PossibleGraphicalGenerator([varId, varZ, varH, varX, varY], index)
	graphGen.makeLatent(varH)
	if (linksId):
		graphGen.edgeFix(varId, varZ)
		graphGen.edgeFix(varId, varH)
	else:
		graphGen.edgeProhibit(varId, varZ, bidir=True)
		graphGen.edgeProhibit(varId, varH, bidir=True)
	graphGen.edgeProhibit(varId, varX, bidir=True)
	graphGen.edgeProhibit(varId, varY, bidir=True)
	graphGen.makeDetermines(varZ, [varId])
	graphGen.makeDetermines(varH, [varId])
	return graphGen.generate()

def possibleManyMany():
	varIdA = ("A", "idA")
	varX = ("A", "X")
	varZ = ("A", "Z")
	varIdB = ("B", "idB")
	varY = ("B", "Y")
	varW = ("B", "W")
	graphGen = PossibleGraphicalGenerator([varIdA, varZ, varX, varIdB, varW, varY])
	graphGen.makeLatent(varZ)
	graphGen.makeLatent(varW)
	graphGen.edgeFix(varIdA, varX)
	graphGen.edgeFix(varIdA, varZ)	
	graphGen.edgeFix(varIdB, varY)
	graphGen.edgeFix(varIdB, varW)	
	graphGen.edgeProhibit(varIdA, varY, bidir=True)
	graphGen.edgeProhibit(varIdA, varW, bidir=True)
	graphGen.edgeProhibit(varIdA, varIdB, bidir=True)
	graphGen.edgeProhibit(varIdB, varX, bidir=True)
	graphGen.edgeProhibit(varIdB, varZ, bidir=True)
	graphGen.makeDetermines(varX, [varIdA])
	graphGen.makeDetermines(varZ, [varIdA])
	graphGen.makeDetermines(varY, [varIdB])
	graphGen.makeDetermines(varW, [varIdB])	
	graphGen.edgeProhibit(varX, varZ, bidir=True)
	graphGen.edgeProhibit(varY, varW, bidir=True)
	graphGen.edgeProhibit(varZ, varW, bidir=True)
	return graphGen.generate()




def condIndys(gDir):
	"""Return a list of independence relations for all measured variables in 
	the form of tuples: (s, t, condSet, isIndy)"""
	sys.stderr.write("\nfinding cond indys for graph %s\n" % (gDir.edges()))
	for e in sorted(gDir.edges()):
		sys.stderr.write("%s -> %s\n" % (e[0][1], e[1][1]))
	g = gDir.to_undirected()
	nodes = g.nodes()

	determines = dict([ [n[0], n[1]["determines"]] for n in g.nodes(data=True) if (n[1]["determines"] is not None) ])
	latents = set([ n for n in g.nodes() if g.node[n]["latent"] ])
	condSets = [ set(x) for x in powerSet(nodes) if (len(latents.intersection(x)) == 0) ]
	measureds = [ n for n in nodes if (n not in latents) ]

	indys = {} # (s1, t1, cons) -> True
	for i, s in enumerate(measureds):
		for t in measureds[i+1:]:
			for condSet in condSets:
				if (s not in condSet) and (t not in condSet):
					detCondS = (s in determines) and (determines[s] <= condSet)
					detCondT = (t in determines) and (determines[t] <= condSet)
					if not (detCondS or detCondT):
						tup = (s, t, tuple(condSet)) if (s < t) else (t, s, tuple(condSet))
						indys[tup] = True # assume conditional independence

	paths = findAllPaths(g)
	if (latents):
		paths = [ p for p in paths if ((p[0] not in latents) and (p[-1] not in latents)) ]
	for path in paths:
		s = path[0]
		t = path[-1]
		for condSet in condSets:
			tup = (s, t, tuple(condSet)) if (s < t) else (t, s, tuple(condSet))
			if (tup in indys) and (dconn(gDir, path, condSet)):
				indys[tup] = False

	indyTups = []
	for key in sorted(indys.keys()):
		s, t, condSet = key 
		tup = (s, t, condSet, indys[key])
		indyTups.append(tup)
	return sorted(indyTups)  
	
def condIndy(g, x, y, conds=[]):
	"""Check the conditional independence of a given pair of variables."""
	for path in findAllPathsSingle(g, x, y):
		if (dconn(g, path, conds)):
			return True
	return False
		
def dconn(gDir, path, conds):
	"""Check whether the path d-connects start and end nodes."""
	for idx in range(len(path) - 2):
		p1 = path[idx]
		p2 = path[idx + 1]
		p3 = path[idx + 2]
		if (gDir.has_edge(p1, p2) and gDir.has_edge(p3, p2)): # p1 -> p2 <- p3
			if (p2 not in conds):
				if not (len(set(nx.dfs_successors(gDir, p2)).intersection(conds)) > 0):					
					return False
		else:
			det = gDir.node[p2]["determines"]
			if (p2 in conds) or ((det is not None) and (det <= conds)):
				return False	
	return True
		
def findAllPaths(gOrig, excludes=[]):
	"""Exhaustively identify all paths from every node to every node within an undirected graph
	and return as a list."""
	g = gOrig.copy()	
	for ex in excludes:
		g.remove_node(ex)
	nodes = g.nodes()	
	pathsAll = []
	for source in nodes:
		path0 = [source]
		pathsCurrent = [path0]
		for radius in range(1, len(nodes)):
			pathsExtended = extendPaths(g, pathsCurrent)
			pathsAll.extend(pathsExtended)
			pathsCurrent = pathsExtended
	return pathsAll

def findAllPathsSingle(gOrig, source, dest=None, excludes=[]):
	"""Get all paths from a single point of origin."""
	if (excludes):
		g = gOrig.copy()	
		for ex in excludes:
			g.remove_node(ex)		
	else:
		g = gOrig	
	pathsAll = []
	path0 = [source]
	pathsCurrent = [path0]
	for radius in range(1, len(g)):
		pathsExtended = extendPaths(g, pathsCurrent)
		if (pathsExtended):
			pathsAll.extend(pathsExtended)
			pathsCurrent = pathsExtended
		else:
			break	
	if (dest is None):
		return pathsAll
	else:
		return [ path for path in pathsAll if (path[-1] == dest) ]

def extendPaths(g, paths):
	"""Helper method for findAllPathsSingle()"""
	rets = []
	for path in paths:
		extendeds = extendPath(g, path)
		rets += extendeds	
	return rets

def extendPath(g, path):
	# sys.stderr.write("\textendPath() for %s\n" % str(path))
	"""Helper helper method for extendPaths(), findAllPathsSingle()"""
	rets = []
	used = set(path)
	for next in g.neighbors(path[-1]):
		if next not in used:
			rets.append(path + [next])				
	return rets


def groundGraph(graphical, dataOrig, typeAttrName="type"):
	"""From a directed graphical model and a(n) (un)directed data graph, produce a ground graph."""
	class GroundNode(object):
		def __init__(self, dataNode, dataNodeType, attrName):
			self.dataNode = dataNode
			self.dataNodeType = dataNodeType
			self.attrName = attrName
		def getEntVar(self):
			return ("%s (%s)" % (self.dataNode, self.dataNodeType), self.attrName)
		def getDataNode(self):
			return self.dataNode

	data = dataOrig.to_undirected()
	ground = nx.DiGraph()

	attrsByType = {} # type -> attr
	for ent, var in graphical.nodes():
		attrsByType.setdefault(ent, []).append(var)
	
	# create a ground node for every attr of every data node
	groundNodesByEntVar = {} # (A, x) -> [ (joe, A, x), (clyde, A, x), ... ]
	for dataNode in data.nodes():
		dataNodeType = data.node[dataNode][typeAttrName]
		for attrName in attrsByType[dataNodeType]:
			groundNode = GroundNode(dataNode, dataNodeType, attrName)
			ground.add_node(groundNode.getEntVar(), type=dataNodeType)
			groundNodesByEntVar.setdefault((dataNodeType, attrName), []).append(groundNode)
	for entVarS, entVarT in graphical.edges():
		entS, varS = entVarS
		entT, varT = entVarT
		if (entS == entT): # intra-node influence
			for groundNodeS in groundNodesByEntVar[(entS, varS)]:
				groundNodeS.getDataNode()
				groundNodeT = GroundNode(dataNode, entT, varT)
				ground.add_edge(groundNodeS.getEntVar(), groundNodeT.getEntVar())
		else: # inter-node influence
			for groundNodeS in groundNodesByEntVar[(entS, varS)]:
				dataNodeS = groundNodeS.getDataNode()
				for dataNodeT in data.neighbors(dataNodeS):
					groundNodeT = GroundNode(dataNodeT, entT, varT)
					ground.add_edge(groundNodeS.getEntVar(), groundNodeT.getEntVar())					
	return ground		




def getGraphicalEquivClasses(graphs, verbose=False):
	"""Given a list of DAGs represented by networkx.DiGraph objects, calculate Markov
	equivalence classes based on conditional independence tests."""
	equivClasses = {} # indys -> [ graph1, graph2, ... ]
	for g, graph in enumerate(graphs):
		independences = tuple(condIndys(graph))		
		if (verbose):
			sys.stderr.write("graph %d\n" % g)
			sys.stderr.write("\t%s\n" % str(graph))
			for indy in independences:
				sys.stderr.write("\t" + str(indy) + "\n")
		equivClasses.setdefault(independences, []).append(graph)
	sys.stderr.write("found %d equivalence classes for %d graphs\n" % (len(equivClasses), len(graphs)))
	if (verbose):
		equivs = equivClasses.keys()
		for i, equiv in enumerate(equivs):
			for condIndy in equiv:
				for j, equivOther in enumerate(equivs):
					if condIndy not in equivOther:
						sys.stderr.write("\tcond indy '%s' is in equiv %d, not in %d\n" % (str(condIndy), i, j))
	return equivClasses

def filterGraphicalEquivClasses(equivClasses, s, t, condSet, isIndy):
	"""Filter graphical equiv classes by the existence of a given relationship."""
	tup = (s, t, tuple(condSet), isIndy) if (s < t) else (t, s, tuple(condSet), isIndy)
	sys.stderr.write("filtering equiv classes with tup: %s\n" % str(tup))
	rets = {}
	for i, (independences, graphicals) in enumerate(equivClasses.items()):
		sys.stderr.write("class %d independences: %s\n" % (i, str(independences)))
		if (tup in independences):
			rets[independences] = graphicals
			sys.stderr.write("keeping %d graphicals in equiv class %d\n" % (len(graphicals), i))
		else:
			sys.stderr.write("discarding %d graphicals in equiv class %d\n" % (len(graphicals), i))
	return rets

def getCommonGraph(graphicals): 
	"""Calculate a "summary graph" for a set of DAGs.  In the summary graph, each
	edge is annotated to indicate whether it exists in all, some, or none of the
	input DAGs."""
	graphCommon = nx.DiGraph()

	# assume that all graphicals have same vars
	vertices = graphicals[0].nodes()
	for node in vertices:
		graphCommon.add_node(node, order=graphicals[0].node[node]['order'])

	num = len(graphicals)
	for i in range(len(vertices)):
		s = vertices[i]
		for j in range(i+1, len(vertices)):
			t = vertices[j]
			countForward = 0
			countBackward = 0
			countMissing = 0
			for graphical in graphicals:
				if (graphical.has_edge(s, t)):
					countForward += 1
				elif (graphical.has_edge(t, s)):
					countBackward += 1
				else:
					countMissing += 1
			if (countForward == num):
				graphCommon.add_edge(s, t, style="directed")
			elif (countBackward == num):
				graphCommon.add_edge(t, s, style="directed")
			elif (countMissing == num):
				pass
			elif (countForward+countBackward == num):
				graphCommon.add_edge(s, t, style="undirected")
			elif (countForward+countMissing == num):
				graphCommon.add_edge(s, t, style="dirmiss")
			elif (countBackward+countMissing == num):
				graphCommon.add_edge(t, s, style="dirmiss")
			else:
				graphCommon.add_edge(s, t, style="undirmiss")
	return graphCommon




# 
# Drawing functions
#

def drawEntity(ax, x, y, rad, variables, latents=[], determined=[]):
	"""Draw a single entity. Returns a dict containing x,y coordinates for each
	variable (for use causal edge drawing)."""
	boxWidth = 3*rad*len(variables) + rad
	boxHeight = 4*rad
	rect = mpatches.Rectangle((x, y), boxWidth, boxHeight, ec="gray", fc="#FFFFFF")
	ax.add_patch(rect)
	variableCoords = {}
	for i, v in enumerate(variables):
		circX = x + rad*(2 + i*3) 
		circY = y + rad*2
		line = "dotted" if (v in latents) else "solid"
		circ = mpatches.Circle(xy=(circX, circY), radius=rad,ec="black", fc="none", linestyle=line)
		variableCoords[v] = (circX, circY)
		ax.add_patch(circ)
		ax.text(circX, circY, v, ha="center", va="center", fontsize=9)
		if (v in determined):
			circ2 = mpatches.Circle(xy=(circX, circY), radius=rad*0.75,ec="black", fc="none", linestyle=line)
			variableCoords[v] = (circX, circY)
			ax.add_patch(circ2)
	return variableCoords

def drawArrow(ax, coordsS, coordsT, head, shrink=0, style="directed", width=1):
	"""Given the coordinates of to variables, draw a causal edge between them."""
	# need to do a little trig to make the arrows look pretty
	distX = coordsT[0] - coordsS[0]
	distY = coordsT[1] - coordsS[1]
	theta = math.atan2(distY, distX)
	offX = shrink*math.cos(theta)
	offY = shrink*math.sin(theta)			
	dx = coordsT[0] - coordsS[0] - 2*offX
	dy = coordsT[1] - coordsS[1] - 2*offY	
	if (coordsS[1] == coordsT[1]):
		rad = -0.75
	else:
		rad = -0.2
	connStyle = mpatches.ConnectionStyle.Arc3(rad=rad)
	if (style == 'directed') or (style == 'dirmiss'):
		arrStyle = "-|>"
	elif (style == 'undirected') or (style == 'undirmiss'):
		arrStyle = "-"
	if (style == "directed") or (style == "undirected"):
		linStyle = 'solid'
	else:
		linStyle = 'dashed'
	mutScale = 10 + 5*(width - 1)
	arr = mpatches.FancyArrowPatch((coordsS[0] + offX, coordsS[1] + offY), (coordsT[0] - offX, coordsT[1] - offY), connectionstyle=connStyle, arrowstyle=arrStyle, mutation_scale=mutScale, ec="black", fc="black", linestyle=linStyle, linewidth=width)
	ax.add_patch(arr)
	
def drawGraphical(ax, x, y, graphical, shrink=1.0, isDirected=True, arrowWidth=1):
	"""Draw a DAG."""
	entToVars = {}
	entToLatents = {}
	entToDetermined = {}
	for n in graphical.nodes(True):
		ent, variable = n[0]
		entToVars.setdefault(ent, []).append(variable)
		if (n[1].get("latent", False)):
			entToLatents.setdefault(ent, set()).add(variable)
		if (n[1].get("determines", False)):
			entToDetermined.setdefault(ent, set()).add(variable)

	radius = shrink / (len(entToVars)*6 - 2)
	coords = {}
	maxVars = max([ len(v) for v in entToVars.values() ])
	for i, ent in enumerate(sorted(entToVars.keys(), reverse=True)):
		gap = (maxVars - len(entToVars[ent]))/2.0
		entX = x + (gap*radius*3)
		entY = y + radius*6*i
		variables = sorted(entToVars[ent], key=lambda n: graphical.node[(ent, n)]['order'])
		coordsEnt = drawEntity(ax, entX, entY, radius, variables, entToLatents.get(ent, []), entToDetermined.get(ent, []))
		for var, xy in coordsEnt.items():
			coords[(ent, var)] = xy		
		# relational connector
		if (i < (len(entToVars) -1)):
			y1 = entY + radius*4
			y2 = y1 + radius*2
			x1 = entX + radius*2
			x2 = x1
			ax.add_line(Line2D([x1, x2], [y1, y2], color="gray", linestyle=":", linewidth=2))
	for s, t, attrs in graphical.edges(data=True):
		drawArrow(ax, coords[s], coords[t], radius/2, 0.8*radius, attrs.get("style", "directed"), arrowWidth)

# equivClasses: independences -> [ graphical1, graphical2, ... ]
def drawGraphicalEquivClasses(equivClasses, outfileName, rowsPerPage=4, colsPerPage=4, pageWidth=8.5, pageHeight=11.0):
	"""Draw the set of Markov equivalence classes.  Each class is listed with a list of conditional
	independence relationships as well as a summary graph for all DAGs in the class."""

	pdf = PdfPages(outfileName)
	def newFig(oldFig=None):
		if (oldFig is not None):
			pdf.savefig(oldFig)
		newFig = plt.figure(figsize=(pageWidth, pageHeight), frameon=False)
		newFig.subplots_adjust(left=0.05, right=0.95, top=0.975, bottom=0.05, wspace=0.035, hspace=0.035)
		return newFig
	
	# translate the base-3 model indices
	modelIdxs = []
	for graphicals in equivClasses.values():
		modelIdxs.extend([ g.graph['index'] for g in graphicals ])
	modelIdx_modelNum = dict([ (idx, num) for num, idx in enumerate(sorted(modelIdxs)) ])
		
	rowCurr = rowsPerPage
	fig = None
	for classNum, independences in enumerate(sorted(equivClasses.keys(), key=lambda x: min([ g.graph['index'] for g in equivClasses[x] ]))):
		graphicals = equivClasses[independences]
		sys.stderr.write("class %d, %d graphicals\n" % (classNum, len(graphicals)))		
		if (rowCurr >= rowsPerPage):
			fig = newFig(fig)
			rowCurr = 0
		rowStart = rowCurr
		axSubIdx = 1 + rowCurr*(colsPerPage + 2)
		sys.stderr.write("class %d add_subplot(%d, %d, %d)\n" % (classNum, rowsPerPage, colsPerPage+2, axSubIdx))
		axSub = fig.add_subplot(rowsPerPage, colsPerPage + 2, axSubIdx, frameon=False, xticks=[], yticks=[]) 
		
		# this is a lot of work to figure out how far down the indy statements should start
		inchesPerRow = pageHeight / rowsPerPage
		skipSize = 0.5 # inches
		#  1.0		  s
		# ----		----
		# 4.3in		0.5in
		# s*inchesPerRow = 1.0*skipSize
		s = skipSize / inchesPerRow 

		axSub.text(0, 1.0, "Class %s (%d)" % (romanize(classNum + 1), len(graphicals)), va="top", fontweight='bold', fontsize=24)
		axSub.text(0.05, 1.0-s, indyTups2String(independences, indyOnly=False), fontsize=18, va="top")
	
		colCurr = 0
		for graphical in sorted(graphicals, key=lambda x: x.graph['index']):
			colCurr += 1
			if (colCurr > colsPerPage):
				colCurr = 1
				rowCurr += 1
				if (rowCurr >= rowsPerPage):
					fig = newFig(fig)
					rowCurr = 0
			axSubIdx = rowCurr*(colsPerPage + 2) + 1 + colCurr
			axSub = fig.add_subplot(rowsPerPage, colsPerPage + 2, axSubIdx, frameon=False, xticks=[], yticks=[], aspect="equal", anchor="NW")				
			#axSub.text(0, 1.0, "model %d" % graphical.graph['index'], va="top")
			axSub.text(0, 1.0, "model %d" % (modelIdx_modelNum[graphical.graph['index']]), va="top")
			drawGraphical(axSub, 0, 0.125, graphical, shrink=0.75)			
			
		axSubIdx = (rowStart + 1)*(colsPerPage + 2)
		axSub = fig.add_subplot(rowsPerPage, colsPerPage + 2, axSubIdx, frameon=True, xticks=[], yticks=[], aspect="equal", anchor="NW", axis_bgcolor='#EEEEEE')				
		#axSub.text(0, 1.0, "common edges", va="top")
		graphCommon = getCommonGraph(graphicals)
		drawGraphical(axSub, 0.125, 0.125, graphCommon, shrink=0.75, arrowWidth=2)			
		
		rowCurr += 1

	pdf.close()
	fig.savefig(outfileName, format="pdf", bbox_inches='tight', pad_inches=0.5)
	#os.system("open " + outfileName)


def indyTups2String(indyTups, indyOnly=False):
	"""Represent a set of conditional independence relationships as a formatted string."""
	ret = ""
	for i, indyTup in enumerate(indyTups):
		s, t, condSet, indy = indyTup
		if (indy):
			ret +=  " %s $ \perp $ %s " % (s[1], t[1])
			if (condSet):	
				ret += " | " + ",".join([ c[1] for c in condSet])
			ret += "\n"
		else:
			if not (indyOnly):
				ret +=  " %s $ \leftrightarrow $ %s " % (s[1], t[1])
				if (condSet):	
					ret += " | " + ",".join([ c[1] for c in condSet])
				ret += "\n"
	return ret[:-1]


# 
# Misc utility functions
#

def factorialList(*paramLists):	
	"""Returns all possible combinations of the args given"""
	combos = [ [] ]
	paramListsRev = [x for x in paramLists ]
	paramListsRev.reverse()
	for paramList in paramListsRev:
		combosNew = []
		for paramVal in paramList:
			for combo in combos:
				comboNew = [paramVal] + combo  
				combosNew.append(comboNew)
		combos = combosNew
	return [ tuple(x) for x in combos ]

def factorialDict(paramListDict):	
	"""Returns list of dicts param -> setting for all combinations of args"""
	params = sorted(paramListDict.keys())
	paramLists = [ paramListDict[x] for x in params ]
	comboTups = factorialList(*paramLists)
	comboDicts = []
	for comboTup in comboTups:
		comboDicts.append(dict(zip(params, comboTup)))
	return comboDicts

def powerSet(lstguy):
	"""Returns the powerset of list.  Would be shorter in Scheme."""
	if (len(lstguy) == 0):
		return [[]]
	else:
		car = lstguy[0]
		cdr = lstguy[1:]
		smaller = powerSet(cdr)
		return smaller + [ [car] + s for s in smaller ] 	

def romanize(n):
	"""Turns an integer into a Roman numeral.  At what price beauty?"""
	s = str(n)
	roman = ""
	if (n >= 1000):
		roman += "m"*int(s[:-3])
	for i in [3, 2, 1]:
		if (len(s) >= i):
			charOne, charFive, charTen = [("i", "v", "x"), ("x", "l", "c"), ("c", "d", "m")][i-1]
			digit = int(s[-i])
			if (digit < 4):
				roman += charOne*digit
			elif (digit == 4):
				roman += charOne + charFive
			elif (digit == 9):
				roman += charOne + charTen
			else: # (digit > 4) and (digit < 9)
				roman += charFive + charOne*(digit - 5)
	return roman



		
###################################################

if (__name__ == '__main__'):

	outfileName = sys.argv[1]
	possibleGraphs = possibleOneManyId(latentZ=True, latentId=False, indexId=False)
	
	equivClasses = getGraphicalEquivClasses(possibleGraphs, True)
	drawGraphicalEquivClasses(equivClasses, outfileName, rowsPerPage=20, colsPerPage=4, pageWidth=16, pageHeight=50)
	




