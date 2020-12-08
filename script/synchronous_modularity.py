#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Prototype of the Label propagation algorithm for community detection based on belows
[1] Raghavan, Usha Nandini, Réka Albert, and Soundar Kumara. "Near linear time algorithm to detect community structures in large-scale networks." Physical review E 76.3 (2007): 036106.
[2] Cordasco, Gennaro, and Luisa Gargano. "Community detection via semi-synchronous label propagation algorithms." 2010 IEEE International Workshop on: Business Applications of Social Network Analysis (BASNA). IEEE, 2010.
"""


import sys
import numpy as np
import networkx as nx
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import itertools
xrange=range

def partition(G):
    #c_iter = community.asyn_lpa_communities(G=G, weight=None)   # Asynchronous
    c_iter = nx.algorithms.community.label_propagation_communities(G=G)   # Semi-synchronous
    max_k_w = []
    for c in c_iter:
        max_k_w += [c]
    community_num_group = len(max_k_w)
    color_list_community = [[] for i in range(len(G.nodes()))]
    for i in range(len(G.nodes())):
        for j in range(community_num_group):
            if i in max_k_w[j]:
                color_list_community[i] = j
    fig = plt.figure()
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]

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
    s = partition(G)
    node_size = 5
    pos = nx.spring_layout(G)
    s = partition(G)
    nx.draw(G, pos, node_color=s, cmap=plt.cm.Blues)
    #nx.draw_networkx_nodes(G, pos, node_size=node_size,  cmap='jet', vmin=0, vmax=max(color_list_community))
    #nx.draw_networkx_edges(G, pos, width=1)
    #nx.draw_networkx_labels(G, pos, font_size=1, font_color="black")
    #plt.xticks([])
    #plt.yticks([])
    plt.savefig("Semi-synchronous.partitioned.png")
    plt.gca().clear()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ("Usage:", file=sys.stderr)
        print ("\t{} karate/out.ucidata-zachary".format(sys.argv[0]), file=sys.stderr)
        print ("\t{} rice/out".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
