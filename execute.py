# Python program for Kruskal's algorithm to find
# Minimum Spanning Tree of a given connected,
# undirected and weighted graph
import sys
from collections import defaultdict
import numpy as np


# Class to represent a graph
class Graph:

    def __init__(self, vertices):
        self.V = vertices  # No. of vertices
        self.graph = []  # default dictionary
        self.primGraph = [[0 for column in range(vertices)]
                          for row in range(vertices)]
        # to store graph

    # function to add an edge to graph
    def add_edge(self, u, v, w):
        self.graph.append([u, v, w])

        # A utility function to find set of an element i

    def degree_inv(self, vertex):
        degrees = 0
        for i in self.graph:
            if i[1] == vertex:
                degrees += 1
        return degrees

    def degree_all(self):
        for i in range(self.V):
            print("the degree of", i, "is", self.degree_inv(i))

    def check_euler(self):
        for i in range(self.V):
            if self.degree_inv(i):
                print("This graph is NOT an euler circuit")
                return False
        print("This graph is an euler circuit")
        return True

    # (uses path compression technique)
    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])

        # A function that does union of two sets of x and y

    # (uses union by rank)
    def union(self, parent, rank, x, y):
        xroot = self.find(parent, x)
        yroot = self.find(parent, y)

        # Attach smaller rank tree under root of
        # high rank tree (Union by Rank)
        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot

            # If ranks are same, then make one as root
        # and increment its rank by one
        else:
            parent[yroot] = xroot
            rank[xroot] += 1

    # The main function to construct MST using Kruskal's
    # algorithm
    def KruskalMST(self):

        result = []  # This will store the resultant MST

        i = 0  # An index variable, used for sorted edges
        e = 0  # An index variable, used for result[]

        # Step 1:  Sort all the edges in non-decreasing
        # order of their
        # weight.  If we are not allowed to change the
        # given graph, we can create a copy of graph
        self.graph = sorted(self.graph, key=lambda item: item[2])

        parent = []
        rank = []

        # Create V subsets with single elements
        for node in range(self.V):
            parent.append(node)
            rank.append(0)

            # Number of edges to be taken is equal to V-1
        while e < self.V - 1:

            # Step 2: Pick the smallest edge and increment
            # the index for next iteration
            u, v, w = self.graph[i]
            i = i + 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            # If including this edge does't cause cycle,
            # include it in result and increment the index
            # of result for next edge
            if x != y:
                e = e + 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)
                # Else discard the edge

        # print the contents of result[] to display the built MST
        print("Following are the edges in the constructed MST by Kruskal")
        sum_weight = 0
        for u, v, weight in result:
            # print str(u) + " -- " + str(v) + " == " + str(weight)
            sum_weight += weight
            print("%d --> %d weight %d" % (u, v, weight))
        # Driver code
        print("Sum of weight", sum_weight)

    def minKey(self, key, mstSet):

        # Initilaize min value
        min = sys.maxsize

        for v in range(self.V):
            if key[v] < min and mstSet[v] == False:
                min = key[v]
                min_index = v

        return min_index

        # Function to construct and print MST for a graph

    # represented using adjacency matrix representation
    def primMST(self):

        # Key values used to pick minimum weight edge in cut
        key = [sys.maxsize] * self.V
        parent = [None] * self.V  # Array to store constructed MST
        # Make key 0 so that this vertex is picked as first vertex
        key[0] = 0
        mstSet = [False] * self.V

        parent[0] = -1  # First node is always the root of

        for cout in range(self.V):

            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is always equal to src in first iteration
            u = self.minKey(key, mstSet)

            # Put the minimum distance vertex in
            # the shortest path tree
            mstSet[u] = True

            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shotest path tree
            for v in range(self.V):
                # graph[u][v] is non zero only for adjacent vertices of m
                # mstSet[v] is false for vertices not yet included in MST
                # Update the key only if graph[u][v] is smaller than key[v]
                if 0 < self.primGraph[u][v] < key[v] and mstSet[v] == False:
                    key[v] = self.primGraph[u][v]
                    parent[v] = u
        print("Following are the edges in the constructed MST by Prim")
        self.printMST(parent)

    # A utility function to print the constructed MST stored in parent[]
    def printMST(self, parent):
        sum_weight = 0
        print("Edge \tWeight")
        for i in range(1, self.V):
            sum_weight += self.primGraph[i][parent[i]]
            print(parent[i], "-", i, "\t", self.primGraph[i][parent[i]])
        print("Sum of weight", sum_weight)


def main(filepath):
    data = np.loadtxt(filepath)
    print("The input data is")
    print(data)
    g = create_graph(data)
    g.primGraph = data
    g.degree_all()
    print()
    g.check_euler()
    print()
    g.KruskalMST()
    print()
    g.primMST()


def create_graph(data):
    graph = Graph(len(data))
    for i in range(len(data)):
        for j in range(len(data[i])):
            if data[i][j] > 0:
                graph.add_edge(i, j, data[i][j])
    return graph


main('Matrix_Group10_A.txt')
print()
main('Matrix_Group10_B.txt')

# g.addEdge(0, 1, 10)
# g.addEdge(0, 2, 6)
# g.addEdge(0, 3, 5)
# g.addEdge(1, 3, 15)
# g.addEdge(2, 3, 4)
#
# g.KruskalMST()
