import sys
import string
import numpy as np
from pyspark import SparkConf, SparkContext

conf = SparkConf()
sc = SparkContext(conf=conf)

lines = sc.textFile(sys.argv[1])
data = lines.map(lambda line: np.fromstring(line,dtype=float, sep=' '))

inputCentroids = sc.textFile(sys.argv[2])
centroids = inputCentroids.map(lambda centroid: np.fromstring(centroid,dtype=float, sep=' '))
centroids = centroids.zipWithIndex()
centroids = centroids.map(lambda c: (c[1],c[0]))

def findCentroid(A):
    distances = list(A)
    closest = distances[0]

    for x in distances:
        if x[1] < closest[1]:
            closest = x

    return closest

def calcEuclideanDist(A,B):
    return np.sqrt(np.sum((A-B)**2))

def calcManhattanDist(A,B):
    return np.sum(np.abs(A-B))


max_iter = 5
cost_measure = []

for i in range(max_iter):

    allPairs = data.cartesian(centroids)
    # Change to calcEuclideanDist for euclidean distance
    pairDist = allPairs.map(lambda pair: (tuple(pair[0]),(pair[1][0],calcEuclideanDist(pair[0],pair[1][1]))))
    pairDist = pairDist.groupByKey()
    pairDist = pairDist.map(lambda point: (findCentroid(point[1]),np.asarray(point[0])))

    allPairs.unpersist()

    errorComp = pairDist.map(lambda point: point[0][1]**2)  # for euclidean add **2
    phi = int(errorComp.reduce(lambda n1,n2: n1+n2))
    cost_measure.append(phi)

    errorComp.unpersist()
    centroids.unpersist()

    centroidGroups = pairDist.map(lambda point: (point[0][0],(point[1], 1)))
    centroidGroups = centroidGroups.reduceByKey(lambda point1, point2: (point1[0]+point2[0],point1[1]+point2[1]))
    centroids = centroidGroups.map(lambda centroid: (centroid[0], centroid[1][0]/centroid[1][1]))
    centroids = centroids.coalesce(1)

    pairDist.unpersist()
    centroidGroups.unpersist()

print('\n\n')
print(cost_measure)
print('\n')
sc.stop()


# Plotting done in jupyter notebooks (code included commented below)

# # Euclidean
# c1 = [616361381, 507045289, 479900973, 458007246, 454990956, 454581137, 454317576, 453943413, 453428248, 452832704,
#       452252105, 451862530, 451576906, 451314212, 450941190, 450551113, 450015083, 448981223, 446698251, 437730979]
# c2 = [401956195, 232605752, 178615711, 153916767, 139802438, 131548078, 123701453, 113252578, 104708368, 98386791, 
#  93050652, 89876836, 87387562, 85989207, 84951189, 84294878, 83801034, 83398486, 83162751, 83068090]
# 
# plt.plot(c1)
# plt.title('Cost Function over Number of Iterations \n Using Euclidean Distance (C1 Intialized)')
# plt.ylabel('Cost Function')
# plt.xlabel('Iteration #')
# plt.show
# 
# plt.plot(c2)
# plt.title('Cost Function over Number of Iterations \n Using Euclidean Distance (C2 Intialized)')
# plt.ylabel('Cost Function')
# plt.xlabel('Iteration #')
# plt.show

# # Manhattan
# c1_M = [510222, 437040, 443295, 455765, 459120, 456269, 452159, 445790, 430613, 415630, 
#         420504, 421559, 421371, 421246, 420447, 421522, 419462, 416532, 415383, 415963]
# c2_M = [1296785, 1005192, 899195, 819434, 785836, 761825, 734131, 699048, 665722, 643036, 
#         629870, 621002, 615305, 610618, 607733, 606184, 603270, 600538, 598977, 597877]

# plt.plot(c1_M)
# plt.title('Cost Function over Number of Iterations \n Using Manhattan Distance (C1 Intialized)')
# plt.ylabel('Cost Function')
# plt.xlabel('Iteration #')
# plt.show

# plt.plot(c2_M)
# plt.title('Cost Function over Number of Iterations \n Using Manhattan Distance (C2 Intialized)')
# plt.ylabel('Cost Function')
# plt.xlabel('Iteration #')
# plt.show

