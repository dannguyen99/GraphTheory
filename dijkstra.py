import csv
import sys

import numpy as np

np.set_printoptions(threshold=np.inf)

cities = []


def form_matrix(file_path, id_path):
    with open(id_path, 'r') as file:
        reader = csv.reader(file, skipinitialspace=True)
        for row in reader:
            cities.append(row[0])
    matrix = np.zeros((len(cities), len(cities)))
    with open(file_path, 'r') as file:
        reader = csv.reader(file, skipinitialspace=True)
        for row in reader:
            index_row = cities.index(row[0])
            index_column = cities.index(row[1])
            matrix[index_row][index_column] = int(row[2])
            matrix[index_column][index_row] = int(row[2])
    return matrix


class Graph:

    def __init__(self, vertices):
        self.graph = None
        self.V = vertices

    def printSolution(self, dist):
        print("Vertex\tDist from Src\tName")
        for node in range(self.V):
            print(node + 1, "\t\t", int(dist[node]), "\t\t", cities[node])

    # A utility function to find the vertex with
    # minimum distance value, from the set of vertices
    # not yet included in shortest path tree
    def minDistance(self, dist, sptSet):

        # Initialize minimum distance for next node
        min_ = sys.maxsize
        min_index = 0
        # Search not nearest vertex not in the
        # shortest path tree
        for v in range(self.V):
            if dist[v] < min_ and not sptSet[v]:
                min_ = dist[v]
                min_index = v

        return min_index

        # Funtion that implements Dijkstra's single source

    # shortest path algorithm for a graph represented
    # using adjacency matrix representation
    def dijkstra(self, src):
        src -= 1
        dist = [sys.maxsize] * self.V
        dist[src] = 0
        sptSet = [False] * self.V

        for i in range(self.V):

            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is al   ways equal to src in first iteration
            u = self.minDistance(dist, sptSet)

            # Put the minimum distance vertex in the
            # shortest path tree
            sptSet[u] = True

            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shortest path tree
            for v in range(self.V):
                if not (not (self.graph[u][v] > 0) or sptSet[v] or not (dist[v] > dist[u] + self.graph[u][v])):
                    dist[v] = dist[u] + self.graph[u][v]

        print("The city as source is", cities[src])
        self.printSolution(dist)


def main(src=1, matrix_path='Vietnam_Distances.txt', id_path='Placed_ID.txt'):
    data = form_matrix(matrix_path, id_path)
    np.savetxt('vietnam_matrix.txt', form_matrix(matrix_path, id_path), fmt='%d')
    g = Graph(np.size(data, 0))
    g.graph = data
    g.dijkstra(src)


main(1)
