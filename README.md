CausalRelational
================

Set of tools for reasoning with causal relational models.

Requires: 
matplotlib: http://matplotlib.sourceforge.net/
networkx: http://networkx.lanl.gov/



MarkovEquivalence.py
=====================================================================================
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
