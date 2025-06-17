# vortexTools.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# Initially created: 23 July 2018
# Last updated: 23 July 2018
#
#
# This is a module that contains methods to identify and visualize vortices
# in fluid flows.







# This follows the method of Sujudi and Haimes (give reference here) to find
# a series of line segments that follows the centerline of vortices within
# a three-dimensional flowfield.  The algorithm assumes it will be supplied
# with a structured set of data.

def findVortexCenterline(x,y,z,u):
  # Import necessary modules
  import numpy as np
  #from numpy.linalg import inv
  from numpy import linalg as la
  import matplotlib.pyplot as plt
  
  
  # Get the i, j, k length of the supplied points of data.
  ni = len(x)
  nj = len(y)
  nk = len(z)
  
  print 'Number of hexahedral elements =',(ni-1)*(nj-1)*(nk-1)
  
  
  
  # This code assumes it will be given structured data in which the points can
  # be arranged into hexahedral elements.  Those elements are then broken into
  # tetrahedral elements as outlined by Sujudi and Haimes.  The hexahedral
  # element vertices are numbered 0 to 7 for each element.  The following
  # list gives the four vertices for each of the six resultant tetrahedral
  # elements.
  tetrahedronVertices = []
  tetrahedronVertices.append((0,1,3,5))
  tetrahedronVertices.append((0,2,3,5))
  tetrahedronVertices.append((2,3,5,7))
  tetrahedronVertices.append((0,4,5,6))
  tetrahedronVertices.append((2,5,6,7))
  tetrahedronVertices.append((0,2,5,6))
  
  
  
  # Each tetrahedral element is made up of four triangular faces.  Here we list 
  # the vertices that make up each triangular face.
  faceVertices = []
  for i in range(6):
      tmp = []
      tmp.append((tetrahedronVertices[i][0],tetrahedronVertices[i][1],tetrahedronVertices[i][2]))
      tmp.append((tetrahedronVertices[i][0],tetrahedronVertices[i][1],tetrahedronVertices[i][3]))
      tmp.append((tetrahedronVertices[i][1],tetrahedronVertices[i][2],tetrahedronVertices[i][3]))
      tmp.append((tetrahedronVertices[i][0],tetrahedronVertices[i][2],tetrahedronVertices[i][3]))
      faceVertices.append(tmp)

  
    
  vortexElements = 0
  totalElements = 0
  vortexPointsX = []
  vortexPointsY = []
  vortexPointsZ = []

  # Loop element by element
  #for i in range(ni-1):
  for i in range(100,101):
      for j in range(nj-1):
          for k in range(nk-1):
              
              totalElements = totalElements + 1
              if (totalElements % 1000 == 0):
                  print 'Total elements processed =',totalElements
                  print '    number of elements with vortices =',vortexElements
  #for i in range(1):
  #    for j in range(1):
  #        for k in range(1):
              
              # Define the hexahedral element vertex coordinates and velocities.
              vertexCoords = np.zeros((8,3))
              vertexCoords[0,:] = (x[i]  ,y[j+1],z[k]  )
              vertexCoords[1,:] = (x[i]  ,y[j]  ,z[k]  )
              vertexCoords[2,:] = (x[i+1],y[j+1],z[k]  )
              vertexCoords[3,:] = (x[i+1],y[j]  ,z[k]  )
              vertexCoords[4,:] = (x[i]  ,y[j+1],z[k+1])
              vertexCoords[5,:] = (x[i]  ,y[j]  ,z[k+1])
              vertexCoords[6,:] = (x[i+1],y[j+1],z[k+1])
              vertexCoords[7,:] = (x[i+1],y[j]  ,z[k+1])
              
              vertexVel = np.zeros((8,3))
              vertexVel[0,:] = (u[:,i  ,j+1,k  ])
              vertexVel[1,:] = (u[:,i  ,j  ,k  ])
              vertexVel[2,:] = (u[:,i+1,j+1,k  ])
              vertexVel[3,:] = (u[:,i+1,j  ,k  ])
              vertexVel[4,:] = (u[:,i  ,j+1,k+1])
              vertexVel[5,:] = (u[:,i  ,j  ,k+1])
              vertexVel[6,:] = (u[:,i+1,j+1,k+1])
              vertexVel[7,:] = (u[:,i+1,j  ,k+1])
              
              if ((i == 100) and (j == 0) and (k==0)):
                  print vertexCoords
                  print ' '
              
              
              # Split the hexahedral element into six tetrahedral elements and
              # apply the algorithm.
              realSumHex = 0
              for m in range(6):
                  
                  # Compute the coefficients of the linear interpolant of the
                  # velocity within the element.  M is the right-hand side
                  # matrix, r is the left-hand vector, a is the set of
                  # solution vectors, and A is the velocity-gradient tensor.
                  M = np.zeros((4,4))
                  r = np.zeros((4,3))
                  a = np.zeros((3,4))
                  for n in range(4):
                      M[n,0:3] = vertexCoords[tetrahedronVertices[m][n]]
                      M[n,3] = 1.0
                      
                      r[n] = vertexVel[tetrahedronVertices[m][n]]
                      
                      
                      
                  if ((i == 100) and (j == 0) and (k==0)):
                     print M
                     print ' '
                          
                  Minv = la.inv(M)
                  a[0,:] = np.matmul(Minv,r[:,0])
                  a[1,:] = np.matmul(Minv,r[:,1])
                  a[2,:] = np.matmul(Minv,r[:,2])
                  A = a[0:3,0:3]
                  
                  # Compute the eigenvalues and eigenvectors of the velocity-
                  # gradient tensor.
                  w,v = la.eig(A)
                  realCheck = np.isreal(w)
                  realSum = 0
                  for p in range(3):
                      if (realCheck[p]):
                          realSum = realSum + 1
                          
                  if (realSum == 3):
                      realSumHex = realSumHex + 1
                      
              if (realSumHex < 6):
                  vortexElements = vortexElements + 1
                  vortexPointsX.append(vertexCoords[0,0])
                  vortexPointsY.append(vertexCoords[0,1])
                  vortexPointsZ.append(vertexCoords[0,2])
                      

  print vortexElements            
              
  return vortexPointsX,vortexPointsY,vortexPointsZ
              
              