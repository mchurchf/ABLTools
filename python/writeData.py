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

def structuredVTK(fileName):
  # Import necessary modules
  import numpy as np
  
  
  # Open the file
  f = open(fileName,'w')
  
  
  # To do here...
  
  
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