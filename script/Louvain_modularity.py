#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Prototype of the Louvain algorithm method. Implemented based on
Fast unfolding of communities in large networks, Vincent D et al., Journal of Statistical Mechanics: Theory and Experiment(2008)
"""


import sys
import numpy as np
import networkx as nx
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import community as community_louvain
import itertools
from numpy.linalg import eigh
xrange=range

def partition(G, refinement=True):
    partition = community_louvain.best_partition(G)

    # draw the graph
    pos = nx.spring_layout(G)
    # color the nodes according to their partition
    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40, cmap=cmap, node_color=list(partition.values()))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.savefig("network_louvain.png")
    plt.gca().clear()



def main(arg):
    """ Main entry
    """
    fp = open(arg)
    G = nx.Graph()
    for row in fp:
        if row[0] == '%': # Comments
            continue
        if len(row.split()) == 2:
            a, b = row.split()
            G.add_edge(a, b)
        else:
            a, b, w = row.split()
            w = int(w)
            G.add_edge(a, b, weight=w)

    print (G.edges())
    s = partition(G, refinement=True)
#    for node, part in zip(itertools.repeat(G.nodes(), s)):
#        print (node, part)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_color=range(G.number_of_nodes()),
            cmap=plt.cm.Blues)
    plt.savefig("graph.png")
    plt.gca().clear()
    nx.draw(G, pos, node_color=s, cmap=plt.cm.jet)
    plt.savefig("graph.partitioned.png")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ("Usage:", file=sys.stderr)
        print ("\t{} karate/out.ucidata-zachary".format(sys.argv[0]), file=sys.stderr)
        print ("\t{} rice/out".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
