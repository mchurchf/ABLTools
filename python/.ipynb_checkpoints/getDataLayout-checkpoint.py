# getDataLayout.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 May 2017
#
# This is a module that contains methods to get the layout of time dependent data
# sampled from atmospheric LES.

 





# If the data is stored with each time sample in its own directory, this determines
# how many time directories there are, what they are called, and puts them in numerical 
#order.

def getOutputTimes(dir):
  # Import necessary modules
  import os


  data = os.listdir(dir)
  outputTimesI = []
  outputTimes = []
  nTimes = len(data)
  ii = 0
  for i in range(nTimes):
     if (data[i][0]).isnumeric():
        outputTimesI.append(data[i])
        ii = ii + 1

  nTimes = len(outputTimesI)

  outputTimesIndex = 0
  outputTimesSort = []
  for i in range(nTimes):
     outputTimesSort.append([i,float(outputTimesI[i])])

  outputTimesSort = sorted(outputTimesSort, key=lambda index: index[1])

  for i in range(nTimes):
     outputTimes.append(outputTimesI[outputTimesSort[i][0]])


  return nTimes, outputTimes
