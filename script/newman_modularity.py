#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Prototype of the Newman modularity inference method. Implemented based on
Newman. (2006) Modularity and community structure in networks.
"""


import networkx as nx


def main():
    """ Main method
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
    A = nx.adjacency_matrix(G)
    print A.todense()


if __name__ == '__main__':
    main()
