#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Prototype of the Newman modularity inference method. Implemented based on
Newman. (2006) Modularity and community structure in networks.
"""


import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

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
    print "A total of {} nodes and {} edges".format(n, m)
    k = np.resize(k, (n, 1)).astype(np.float64)
    B = A - np.dot(k, k.T) / (2 * m)
    #print B.sum(axis=0)  # Make sure this is almost all zero
    print B
    # Only interested in the largest eigenvector
    #evals_large, evecs_large = eigh(B, eigvals=(N - 1, N - 1))
    w, v = eigh(B)
    z = np.argmax(w)
    print "Largest eigenvalue:", w[z]
    s = np.sign(v[:, z])
    orig_score = evaluate(s, n, m, B)
    print "Q:", orig_score

    # Fine-tuning - iteratively modify the partition of s_i
    seen = set()
    while True:
        ans = []
        for i in xrange(n):
            if i in seen:
                continue
            s[i] = -s[i]
            new_score = evaluate(s, n, m, B)
            ans.append((new_score, i))
            s[i] = -s[i]

        best_score, best_i = max(ans)
        if best_score > orig_score:
            print "ACCEPTED: Q = {}, Q' = {}".format(orig_score, best_score)
            s[best_i] = -s[best_i]
            seen.add(best_i)
            orig_score = best_score
        else:
            break

    return s


def evaluate(s, n, m, B):
    """ Calculate Q
    """
    ps = np.resize(s, (n, 1))
    Q = np.dot(np.dot(ps.T, B), ps) / (4 * m)
    return Q[0][0]


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

    print G.edges()
    s = partition(G)
    for node, part in zip(G.nodes(), s):
        print node, part
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_color=range(G.number_of_nodes()),
            cmap=plt.cm.Blues)
    plt.savefig("graph.png")
    plt.gca().clear()
    nx.draw(G, pos, node_color=s, cmap=plt.cm.jet)
    plt.savefig("graph.partitioned.png")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >> sys.stderr, "Usage:"
        print >> sys.stderr, "\t{} karate/out.ucidata-zachary".format(sys.argv[0])
        print >> sys.stderr, "\t{} rice/out".format(sys.argv[0])
        sys.exit(1)

    main(sys.argv[1])
