import sys
import string
from pyspark import SparkConf, SparkContext

conf = SparkConf()
sc = SparkContext(conf=conf)

lines = sc.textFile(sys.argv[1])
friendList = lines.map(lambda line: line.split('\t'))

# Find potential pairs
friendPairs = friendList.flatMapValues(lambda friends: friends.split(','))
allPossiblePairs = friendPairs.join(friendPairs)
allPossiblePairs = allPossiblePairs.map(lambda pairs: pairs[1]).filter(lambda pairs: pairs[0] != pairs[1])

# Filter for friend suggestions and number of mututal friends
suggestedPairs = allPossiblePairs.subtract(friendPairs)
suggestedPairs = suggestedPairs.map(lambda pairs: (pairs,1))
suggestedPairs = suggestedPairs.reduceByKey(lambda n1, n2: n1+n2)
suggestedPairs = suggestedPairs.map(lambda pairs: (int(pairs[0][0]), (int(pairs[0][1]),pairs[1])))
suggestedPairs = suggestedPairs.groupByKey()
suggestedFriends = suggestedPairs.map(lambda pairs: (pairs[0],topTen(pairs[1])))

# Function that returns a list of the top ten friends in sorted order
def topTen(iterFriendSuggestions):
	friendSuggest = list(iterFriendSuggestions)
	topTen = sorted(friendSuggest, key=lambda k: (k[1],-k[0]), reverse=True)[:10]
	for x in range(len(topTen)):
		topTen[x] = topTen[x][0]
	return topTen

# Creates a final mapping for file to be exported
outputSuggestions = suggestedFriends.map(lambda suggest: str(suggest[0]) + '\t' + ','.join(str(friend) for friend in suggest[1]))

outputSuggestions.saveAsTextFile(sys.argv[2])
sc.stop()