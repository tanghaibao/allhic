#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Prototype of the Newman modularity inference method. Implemented based on
Newman. (2006) Modularity and community structure in networks.
"""


import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#from scipy.linalg import eigh
from numpy.linalg import eigh


def partition(G):
    """ Newman modularity inference method

    Q = 1/4m Sum_ij(A_ij - k_i * k_j / 2m) s_i * s_j

    m: total number of edges
    A_ij: adjancency matrix
    k_i, k_j: degree of node i, j
    s_i, s_j: partition of node i, j, either +1 or -1

    We can conveniently wrtie Q in matrix form

    Q = 1/4m s.T * B * s

    Where B_ij = A_ij - k_i * k_j / 2m
    """
    A = nx.adjacency_matrix(G)
    A = np.array(A.todense())
    print A
    k = A.sum(axis=0)
    m = k.sum() / 2
    n = len(k)
    k = np.resize(k, (n, 1)).astype(np.float)
    B = A - np.dot(k, k.T) / (2 * m)
    #print B.sum(axis=0)  # Make sure this is almost all zero
    print B
    # Only interested in the largest eigenvector
    #evals_large, evecs_large = eigh(B, eigvals=(N - 1, N - 1))
    w, v = eigh(B)
    print w[-1]
    print v[-1]


def main():
    """ Main entry
    """
    fp = open("karate/out.ucidata-zachary")
    G = nx.Graph()
    for row in fp:
        if row[0] == '%': # Comments
            continue
        a, b = row.split()
        a, b = int(a), int(b)
        G.add_edge(a, b)

    print G.edges()
    partition(G)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_color=range(G.number_of_nodes()),
            cmap=plt.cm.Blues)
    plt.savefig("graph.png")


if __name__ == '__main__':
    main()
