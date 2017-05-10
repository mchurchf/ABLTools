# readData.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 May 2017
#
# This is a module that contains methods to read in data sampled from atmospheric LES.







# Read in structured VTK data.  The data is returned as a list called fields that
# contains numpy arrays.  For example, the VTK file may contain both pressure
# and temperature, or velocity and vorticity, so each quantity will be contained
# in a separate array in the list.  Then for each array, the first index is
# the number of components, so for scalars, there is only one component, but for
# vectors there will be three.  The remaining dimensions are the number of 
# sample points in x, y, and z.

def structuredVTK(fileName):
  # Import necessary modules
  import numpy as np
  
  
  # Open the file
  f = open(fileName,'r')
  
  
  # Get data set name.
  f.readline()
  dataSetName = f.readline()
  
  
  # Get data dimensions.
  f.readline()
  f.readline()
  d = f.readline()
  d = d[11:]
  c = d.split()
  dims = []
  dims.append(int(c[0]))
  dims.append(int(c[1]))
  dims.append(int(c[2]))
  dims = np.asarray(dims)
  
  
  # Get data origin.
  d = f.readline()
  d = d[7:]
  c = d.split()
  origin = []
  origin.append(float(c[0]))
  origin.append(float(c[1]))
  origin.append(float(c[2]))
  origin = np.asarray(origin)
  
  
  # Get data spacing in each direction.
  d = f.readline()
  d = d[8:]
  c = d.split()
  spacing = []
  spacing.append(float(c[0]))
  spacing.append(float(c[1]))
  spacing.append(float(c[2]))
  spacing = np.asarray(spacing)
  
  
  # Form data point structured grid.
  if (dims[0] > 1):
      x = np.linspace(origin[0],origin[0]+spacing[0]*(dims[0]-1),dims[0])
  else:
      x = np.array([1])
      x[0] = origin[0]
  if (dims[1] > 1):
      y = np.linspace(origin[1],origin[1]+spacing[1]*(dims[1]-1),dims[1])
  else:
      y = np.array([1])
      y[0] = origin[1]
  if (dims[2] > 1):
      z = np.linspace(origin[2],origin[2]+spacing[2]*(dims[2]-1),dims[2])
  else:
      z = np.array([1])
      z[0] = origin[2]
  
  
  # Read header for field data
  f.readline()
  d = f.readline()
  d = d[18:]
  nFields = int(d)
  
  
  # Read field data.
  field = []
  fieldName = []
  fieldDim = []
  for m in range(nFields):
      if (m > 0):
          f.readline()
          
      d = f.readline()
      c = d.split()
      fieldName.append(c[0])
      fieldDim.append(int(c[1]))
      dataArray = np.zeros((fieldDim[m], dims[0], dims[1], dims[2]))
      for i in range(dims[0]):
          for j in range(dims[1]):
              for k in range(dims[2]):
                  l = f.readline()
                  l = l.split()
                  for n in range(fieldDim[m]):
                      dataArray[n][i][j][k] = float(l[n])
                      
      field.append(dataArray)
  
  
  # Close file.
  f.close()
  
  
  # Return the data.
  return dataSetName, dims, origin, spacing, x, y, z, nFields, fieldName, fieldDim, field
