import sys
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import bernoulli
from random import choice,randrange,randint
from collections import Counter
from statistics import mean

## Initial Parameters
k = 32 #k quadratic (2D) lattice
p = 0.65 #bernoulli probability for bond percolation
q = 0.95 #bernoulli probability for terminal arrival
RepeaterMark = "s" #display repaters as squares, i.e., 's'
TerminalMark = "o" #display terminals as discs, i.e., 'o'


## Processing Parameters
DrawGrid=False
ShowLabels=False
BondPercolation=True
ComputePaths=True
PathSearchAlgorithm=1 #1=shortestPaths 2=peelingPaths
AdditionalRing=True
CSVFormat=True

#Types of nodes, either terminals or repeaters
nodeShapes = [RepeaterMark,TerminalMark]


#compute the intersection of two lists
def intersection(A, B):
    C = [c for c in A if c in B]
    return C

def process_graph(G,k,p,q):
    """
    Given an initial grid_2d_graph, return a copy of the graph with
    bond percolation, peripheral edges as terminals, ...
    :param G: the original graph
    :param k: k quadratic (2D) lattice
    :param p: bernoulli probability for bond percolation
    :param q: bernoulli probability for terminal arrival
    :return: the processed copy of the original graph
    """
    H = G.copy()

    Edges=[]
    for edge in H.edges():
        Edges.append(edge)

    rmEdges=[]

    periphery=[]

    ###### Periphery Removal ######
    if AdditionalRing:
        #remove all peripheral nodes
        peripheryExtra=[]
        for i in range(k):
            periphery.append(i)#1st row [e.g., if k=10, nodes 0 to 9]
            periphery.append(((k-1) * k)+i)#last row [e.g., nodes 90 to 99]

        for j in range(1,k-1):
            periphery.append(j*k)#1st column [e.g., nodes 10 to 80]
            periphery.append((j*k)+(k-1))#1st column [e.g., nodes 19 to 89]
            peripheryExtra.append((j*k)+1)#2nd column [e.g., 11 to 81]
            peripheryExtra.append((j*k)+(k-2))#prev-last column [e.g., 18 to 88]
            peripheryExtra.append(((k-2) * k)+j)#prev-last row [e.g., 81 to 88]

        # get edges to be removed
        rmEdges = [(u,v)
                   for u,v in H.edges()
                   if (u or v) in list(periphery+peripheryExtra) ]

        H.remove_nodes_from(periphery)
        for edge in rmEdges:
            Edges.remove(edge)
    ###### Periphery Removal ######

    ###### Bond Percolation ######
    if BondPercolation:
        # Bond percolation with bernoulli probability p
        for i in range(len(Edges)):
            coinBondPercolation = bernoulli.rvs(p, size=1)
            if (not coinBondPercolation[0]):
                H.remove_edge(Edges[i][0],Edges[i][1]) #edge dissapears
    ###### Bond Percolation ######

    #Select the largest clique structure as the main subgraph
    subH=max(list(nx.connected_component_subgraphs(H)),key=len)

    repeaterNodes=[]
    restoredEdges=[]

    #All nodes in subH (all of them repeaters) are connected
    for node in subH.nodes():
        repeaterNodes.append(node)

    terminalNodes = []

    ###### Periphery Restauration ######
    if AdditionalRing:
        # Periphery restauration with bernoulli probability 1-q
        for i in range(len(rmEdges)):
            coinBondRestauration = bernoulli.rvs(1-q, size=1)
            if (not coinBondRestauration[0]):
                if rmEdges[i][0] not in repeaterNodes:
                    restoredEdges.append(rmEdges[i])
                else:
                    if rmEdges[i][1] not in repeaterNodes:
                        restoredEdges.append(rmEdges[i])

        subH.add_edges_from(restoredEdges) #peripheral edges restored
        #those restored nodes are now terminals

        #remove isolated nodes, if any
        subH = max(list(nx.connected_component_subgraphs(subH)),key=len)

        #Remove adjacent terminals
        for node in list(subH.nodes):
            if node in periphery:
                H1 = nx.ego_graph(G, node, center=False,distance=1)
                adjacentNodes=intersection(list(H1),list(periphery))
                if (len(adjacentNodes) > 0):
                    for adjacentNode in adjacentNodes:
                        ebunch=[(node,adjacentNode)]
                        try:
                            subH.remove_edges_from(ebunch)
                            #print "removing edge between",node,"and",adjacentNode
                        except nx.exception:
                            print "Exception"

        #remove isolated nodes, if any
        subH = max(list(nx.connected_component_subgraphs(subH)),key=len)

        terminalNodes = [u for u in subH.nodes()
                         if u in periphery]

    ###### Periphery Restauration ######

    #Add, as well, those nodes with just one edge, not in periphery
    for node in subH.nodes():
        if node not in periphery:
            if subH.degree(node)==1:
                coinBondRestauration = bernoulli.rvs(1-q, size=1)
                if (not coinBondRestauration[0]):
                    terminalNodes.append(node)


    repeaterNodes=[]
    for node in subH.nodes():
        if node not in terminalNodes:
            repeaterNodes.append(node)
            nx.set_node_attributes(subH, {node:{'shape':RepeaterMark}})
        else:
            nx.set_node_attributes(subH, {node:{'shape':TerminalMark}})

    numRepeaters=len(repeaterNodes)
    numTerminals=len(terminalNodes)

    if ShowLabels:
        print "The graph contains",numRepeaters,"repeaters [",repeaterNodes,"] and",numTerminals,"terminals [",terminalNodes,"]."
        print ""
        print "Paths:"
    elif not CSVFormat:
        print "The graph contains",numRepeaters,"repeaters and",numTerminals,"terminals."
    else:
        print numRepeaters,
        print "\t",
        print numTerminals,
        print "\t",

    ###### Search of Paths ######
    if ComputePaths and len(terminalNodes)<2:
        print "Not enough terminals. Try again."

    elif ComputePaths and PathSearchAlgorithm==1:
        visitedRepeaters=[]
        pathLength=[]
        for i in range(numTerminals):
            for j in range(i+1,numTerminals):
                Found=False
                if ShowLabels:
                    print terminalNodes[i],"->",terminalNodes[j],":",
                try:
                    path=nx.shortest_path(subH,terminalNodes[i],terminalNodes[j]) #This can fail
                    if not intersection(list(path[1:-1]),list(terminalNodes)):#discard paths containing terminals
                        if ShowLabels:
                            print path,
                            print "[shortest path between both terminals]"
                        visitedRepeaters.append(path)
                        pathLength.append(len(path))
                        Found=True
                    else:
                        #print " (2nd try)",
                        allPaths=nx.all_shortest_paths(subH,terminalNodes[i],terminalNodes[j]) ##This can still fail
                        for newpath in allPaths:
                            if not intersection(list(newpath[1:-1]),list(terminalNodes)):#discard paths containing terminals"
                                if ShowLabels:
                                    print path,
                                    print "[found around",len(list(allPaths)),"other shortest paths]"
                                visitedRepeaters.append(newpath)
                                pathLength.append(len(newpath))
                                Found=True
                                break
                            if not Found:
                                #print " (3rd try needed)",
                                allPaths=nx.all_simple_paths(subH,terminalNodes[i],terminalNodes[j])
                                for newpath in allPaths:
                                    if not intersection(list(newpath[1:-1]),list(terminalNodes)):#discard paths containing terminals"
                                        if ShowLabels:
                                            print path,
                                            print "[found around",len(list(allPaths)),"other paths]"
                                        visitedRepeaters.append(newpath)
                                        pathLength.append(len(newpath))
                                        Found=True
                                        break
                                    if not Found:
                                        print "I couldn't find any path from node",terminalNodes[i],"to node",terminalNodes[j],"!"
                                        sys.exit()
                except nx.NetworkXNoPath:
                    print "No path between them!"

    elif ComputePaths and PathSearchAlgorithm==2:
        visitedRepeaters=[]
        pathLength=[]
        REDO=True
        for i in range(numTerminals):
            for j in range(i+1,numTerminals):
                if REDO:
                    H = subH.copy()
                    REDO=False

                if ShowLabels:
                    print terminalNodes[i],"->",terminalNodes[j],":",
                try:
                    path=nx.shortest_path(H,terminalNodes[i],terminalNodes[j])
                    visitedRepeaters.append(path)
                    pathLength.append(len(path))
                    #remove edges in the path
                    path_edges = zip(path,path[1:])
                    path_edges = path_edges[1:-1]#omit first and last pairs
                    H.remove_edges_from(path_edges)
                    if ShowLabels:
                        #print path
                        print "Path=",path,"Edges=",path_edges
                except nx.NetworkXNoPath:
                    #print "Enabling redo, case 1"
                    REDO=True
                except nx.NetworkXException:
                    #print "Enabling redo, case 2"
                    REDO=True

    if ComputePaths:
        summary=[]
        if not CSVFormat:
            print ""
        for repeater in visitedRepeaters:
            summary+=repeater
            #print summary
        congestion=Counter(summary).most_common(3) #show the (three) most congested repeaters
        if ShowLabels:
            print "Congestion = ",congestion[0][-1]," (Repeater",congestion[0][0],"appears in",congestion[0][-1],"paths,",
            print "repeater",congestion[1][0],"appears in",congestion[1][-1],"paths,",
            print "repeater",congestion[2][0],"appears in",congestion[2][-1],"paths, etc.) "
            #print "path lenght is",pathLength
            print "Average length is",int(round(mean(pathLength),3))
            print ""
        elif not CSVFormat:
            print "Congestion = ",congestion[0][-1]
            print "Average length = ",int(round(mean(pathLength),3))
            print ""
            print ""
        else:
            print congestion[0][-1],"\t",int(round(mean(pathLength),3))
    ###### Search of Paths ######


    return subH


############## MAIN ############

G = nx.grid_2d_graph(k,k) #Create a networkx k^2 2D lattice

#Instead of labelling nodes as (i,j), label them as 0,1,..k^2
labels = dict( ((i, j), i + (k-1-j) * k ) for i, j in G.nodes() )
#print "labels=",labels
nx.relabel_nodes(G,labels,False)#relabel the node identifiers
pos = {y:x for x,y in labels.iteritems()}#prepare the springbox for matplotlib

subGraph=process_graph(G,k,p,q)

if DrawGrid:


    if ShowLabels:
        nx.draw(subGraph, pos=pos, node_shape=RepeaterMark, node_color='#000000', with_labels = True, node_size = 1)
        for aShape in nodeShapes:
            nx.draw_networkx_nodes(subGraph, pos, node_shape = aShape, nodelist = [sNode[0] for sNode in filter(lambda x: x[1]["shape"]==aShape,subGraph.nodes(data = True))], node_color='#000000', with_labels = False, node_size = 1)
    else:
        nx.draw(subGraph, pos=pos, node_shape=RepeaterMark, node_color='#000000', with_labels = False, node_size = 1)
        for aShape in nodeShapes:
            nx.draw_networkx_nodes(subGraph, pos, node_shape = aShape, nodelist = [sNode[0] for sNode in filter(lambda x: x[1]["shape"]==aShape,subGraph.nodes(data = True))], node_color='#000000', with_labels = False, node_size = 30)

    plt.axis('off')
    plt.show()
