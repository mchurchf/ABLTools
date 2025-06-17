# writeData.py
#
# Paula Doubrawa
# National Renewable Energy Laboratory
# 9 May 2017
#
# This is a module that contains methods to write out postprocessed data.







# Write out Turbsim format.
def turbSim(fileName):
    
  # Open the file
  f = open(fileName,'w')
  
  
  # To do here...
  
  
  # Close file.
  f.close()   
    
   
  
  
  
  
    
    
# Write out structured VTK data.  
def structuredVTK(fileName,dataName,dims,origin,spacing,fieldDims,fields,fieldNames):
  # Import necessary modules.
  import numpy as np
  
  
  # Open the file.
  f = open(fileName,'w')
  

  # Get the total number of points.
  nPoints = dims[0] * dims[1] * dims[2]
  nFields = len(fields)
    
  # Write the header.
  f.write('# vtk DataFile Version 3.0' + '\n')
  f.write(dataName + '\n')
  f.write('ASCII' + '\n')
  f.write('DATASET STRUCTURED_POINTS' + '\n')
  f.write('DIMENSIONS ' + str(dims[0]) + ' ' + str(dims[1]) + ' ' + str(dims[2]) + '\n')
  f.write('ORIGIN ' + str(origin[0]) + ' ' + str(origin[1]) + ' ' + str(origin[2]) + '\n')
  f.write('SPACING ' + str(spacing[0]) + ' ' + str(spacing[1]) + ' ' + str(spacing[2]) + '\n')
  f.write('POINT_DATA ' + str(nPoints) + '\n')
  f.write(' FIELD attributes ' + str(nFields) + '\n')


  # Write the data.
  for m in range(nFields):
      f.write(fieldNames[m] + ' ' + str(fieldDims[m]) + ' ' + str(nPoints) + ' ' + 'float' + '\n')
      
      for i in range(nPoints):
          fieldString = str(fields[m][i,0])
          for j in range(1,fieldDims[m]):
              fieldString = fieldString + ' ' + str(fields[m][i,j])
          fieldString = fieldString + '\n'
          f.write(fieldString)
      
      f.write('\n')
        
  
  # Close file.
  f.close()  
    
   
  
  
  
  
    
    
# Write out Ensight data.  
def ensight(fileNameMesh,fileNameField,writeFieldOnly,dims,fieldDim,x,y,z,field):
  # Import necessary modules
  import numpy as np
  
  
  # Write the mesh file.
  if (writeFieldOnly == 0):
      f = open(fileNameMesh,'w')
      
      f.write('Ensight Geometry File' + '\n')
      f.write('written by OpenFOAM-2.4.x' + '\n')
      f.write('node id assign' + '\n')
      f.write('element id assign ' + '\n')
      f.write('part' + '\n')
      f.write('         1' + '\n')
      f.write('internalMesh' + '\n')
      f.write('coordinates' + '\n')
      f.write(str(dims) + '\n')
      
      for i in range(dims):
          f.write(str(x[i]) + '\n')
      
      for i in range(dims):
          f.write(str(y[i]) + '\n')
      
      for i in range(dims):
          f.write(str(z[i]) + '\n')
          
      f.write('point' + '\n')
      
      f.write(str(dims) + '\n')
      
      for i in range(dims):
          f.write(str(i) + '\n')
      
      
      f.close()
      
      
      


  


  # Write the field file.
  f = open(fileNameField,'w')
  
  if (fieldDim == 1):
      f.write('scalar' + '\n')
  elif (fieldDim == 3):
      f.write('vector' + '\n')
  elif (fieldDim == 9):
      f.write('tensor' + '\n')
      
  f.write('part' + '\n')
  f.write('         1' + '\n')
  f.write('coordinates' + '\n')
  
  for i in range(fieldDim):
      for j in range(dims):
          f.write(str(field[j,i]) + '\n')
  
    
  f.close()