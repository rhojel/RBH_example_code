## Code for HITS.py
## Calculates the hubbiness and authority of all nodes (uses 40 iterations)
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

n = edgeList.map(lambda k: k[0]).distinct().count()
h = np.ones(n)
a = np.zeros(n)

def update_a(edge, h, a_accum):
    update = np.zeros(len(h))
    update[edge[1]-1] = h[edge[0]-1]
    a_accum.add(update)

def update_h(edge, a, h_accum):
    update = np.zeros(len(h))
    update[edge[0]-1] = a[edge[1]-1]
    h_accum.add(update)

class VectorAccumulatorParam(AccumulatorParam):
    def zero(self, value):
        return np.zeros(len(value))

    def addInPlace(self, v1, v2):
        v1 += v2
        return v1

for i in range(40):

    a_accum = sc.accumulator(a, VectorAccumulatorParam())
    edgeList.foreach(lambda edge: update_a(edge, h, a_accum)) 
    a = a_accum.value

    a = a / np.amax(a)

    h_accum = sc.accumulator(h, VectorAccumulatorParam())
    edgeList.foreach(lambda edge: update_h(edge, a, h_accum)) 
    h = h_accum.value

    h = h / np.amax(h)

asort = np.argsort(a)
# print(list(a))
print("a Worst: ", asort[:5]+1)
print("a Best: ", asort[990:]+1)
hsort = np.argsort(h)
# print(list(h))
print("h Worst: ", hsort[:5]+1)
print("h Best: ", hsort[990:]+1)