## Code for pageRank.py
## Calculates the page rank of every node through 40 power iterations 
import sys
import numpy as np
from pyspark import SparkConf, SparkContext
from pyspark.accumulators import AccumulatorParam

conf = SparkConf()
sc = SparkContext(conf=conf)

lines = sc.textFile(sys.argv[1])
data = lines.map(lambda line: line.split())
edgeList = data.map(lambda l: (int(l[0]),int(l[1])))
edgeList = edgeList.distinct()

degree = data.map(lambda l: (int(l[0]),1))
degree = degree.reduceByKey(lambda n1,n2: n1+n2)

edgeList = edgeList.join(degree)

n = degree.count()
beta = 0.8
l = np.ones(n)
r_old = 1/n * l

def update_r(edge, r_old, accum):
    update = np.zeros(len(r_old))
    update[edge[1][0]-1] = beta * r_old[edge[0]-1]/edge[1][1]
    accum.add(update)

class VectorAccumulatorParam(AccumulatorParam):
    def zero(self, value):
        return np.zeros(len(value))

    def addInPlace(self, v1, v2):
        v1 += v2
        return v1

for i in range(40):

    r_new = (1-beta)/n * l
    accum = sc.accumulator(r_new, VectorAccumulatorParam())

    edgeList.foreach(lambda edge: update_r(edge, r_old, accum)) 
    r_old = accum.value

sort = np.argsort(r_old)
# print(list(r_old))
print("Worst: ", sort[:5]+1)
print("Best: ", sort[990:]+1)