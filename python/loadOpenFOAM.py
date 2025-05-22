# loadOpenFOAM.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 7 May 2025
#
# This is a class to load OpenFOAM native data.  In writing this,
# I borrowed heavily from https://github.com/fluiddyn/fluidfoam/.
# If you want more capability, consider that package.  The purpose
# of this class is to be able to load basic OpenFOAM data so that
# on can manipulate it in python. 


import numpy as np
import string
import struct




class loadOpenFOAM:
    
    # Initialize the class.
    def __init__(self,intType='int64',floatType='double'):
        
        self.intType = intType
        self.floatType = floatType   
        
        self.polyMeshDir = []
        
        self.nPoints = []
        self.points = []
        
        self.nOwners = []
        self.owners = []

        self.nNeighbours = []
        self.neighbours = []
        
        self.nFaces = []
        self.faces = []
        self.faceCenters = []
        self.faceAreas = []
        self.faceNormals = []

        self.nCells = []
        self.cellnFaces = []
        self.cellCenters = []
        self.cellVolumes = []

        self.percentCompleteInterval = 20
        
        if (self.intType == 'int64'):
            self.structFormatInt = "l"
        else:
            self.structFormatInt = "i"


        if (floatType == 'double'):
            self.structFormatFloat = "d"
        else:
            self.structFormatFloat = "f"



    

    
    # Completion message
    def completionMessage(self,string,ii,total,interval):
        if (ii % int(total/interval) == 0):
                percentComplete = round(100.0*(ii/total),1)
                
                if ((percentComplete > 99.0) and (percentComplete <= 101.0)):
                    endCharacter = '\n'
                else:
                    endCharacter = '\r'
                 
                print(string + ': ' + str(percentComplete) + '% complete...', end=endCharacter)





    
    # Compute the cross product between two vectors.
    def cross(self,v1,v2):
        #                |    i    j    k |
        # v3 = v1 x v2 = | v1_0 v1_1 v1_2 |
        #                | v2_0 v2_1 v2_2 |
        #
        # v3 = (v1_1*v2_2 - v2_1*v1_2)*i + (-v1_0*v2_2+v2_0*v1_2)*j + (v1_0*v2_1-v2_0*v1_1)*k

        v3 = np.zeros((3,))
        v3[0] =  v1[1]*v2[2] - v2[1]*v1[2]
        v3[1] = -v1[0]*v2[2] + v2[0]*v1[2]
        v3[2] =  v1[0]*v2[1] - v2[0]*v1[1]
        
        return v3



    

    
    # Compute the dot product between two vectors.
    def dot(self,v1,v2):

        # v3 = v1*v2 = v1_0*v2_0 + v1_1*v2_1 + v1_2*v2_2

        v3 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
        
        return v3



    

    
    # Compute the magnitude of a general vector.
    def mag(self,v):
        nComps = len(v)

        sumSquare = 0.0
        if (nComps > 1):
            for i in range(nComps):
                sumSquare += v[i]*v[i]
        else:
            sumSquare = v*v

        return np.sqrt(sumSquare)



    

    
    # Read and load the mesh information.
    def readMesh(self,polyMeshDir,
                 pointsFile=[],
                 facesFile=[],
                 ownersFile=[],
                 neighboursFile=[],
                 makeStructured=False):
        
        print('Reading the mesh from ' + polyMeshDir)

        # Set the polyMesh directory
        self.polyMeshDir = polyMeshDir

        # Set the filenames to default locations.
        if not pointsFile:
            pointsFile = polyMeshDir + '/points'
        if not facesFile:
            facesFile = polyMeshDir + '/faces'
        if not ownersFile:
            ownersFile = polyMeshDir + '/owner'
        if not neighboursFile:
            neighboursFile = polyMeshDir + '/neighbour'
        
        # Call functions to read in the mesh and compute cell centers/volume and
        # face centers/areas/normals
        self.readPoints(pointsFile)
        self.readFaces(facesFile)
        self.readOwners(ownersFile)
        self.readNeighbours(neighboursFile)
        self.computeFaceCentersAndAreas()
        self.computeCellCentersAndVolumes()

        # Write summary of grid statistics to screen
        print('\n')
        print('Summary:')
        print('    -number of cells:', self.nCells)
        print('    -number of faces:', self.nFaces)
        print('    -number of owners:', self.nOwners)
        print('    -number of neighbours:', self.nNeighbours)
        print('    -number of points:', self.nPoints)

        # If this is a structured mesh, reorganize into structured format.
        if (makeStructured == True):
            print("make structured...")



    

    
    # Read and load the mesh information.
    def findMeshStructure(self,cellCenters,precision=8):

        # The data comes out with differing digits out around the 10th to 12th
        # decimal place, so round the data before sorting
        cellCentersRounded = np.round(cellCenters,decimals=precision)


        # Find the unique coordinates to get the number of points in x, y, and z
        nx = len(np.unique(cellCentersRounded[:,0]))
        ny = len(np.unique(cellCentersRounded[:,1]))
        nz = len(np.unique(cellCentersRounded[:,2]))

        
        # Sort by x, y, and then z, and return the sorting index.
        iSort = np.lexsort((cellCentersRounded[:,2], 
                            cellCentersRounded[:,1],
                            cellCentersRounded[:,0]))
        

        return nx,ny,nz,iSort


    

    

    
    # Read the points file and load points information.
    def readPoints(self,pointsFile):
        print('    -Reading the points file...')
        

        # Open the points file (currently only coded to read binary)
        fileID = open(pointsFile,'rb')
        content = fileID.read()
        
        
        # Split the data into lines and then find the number of points.
        lines = content.split(b"\n")
        self.nPoints = int(lines[18])
        
        
        # Split the data where the actual points begin using .split, and
        # then keep only the non-header data (the [1:] does this), and 
        # then join the data back into a single entity with .join.
        data = b"\n(".join(content.split(b"\n(")[1:])
        
        
        # Use struct unpack to take the binar data and turn it into a list of
        # of doubles, and then put those into a numpy array.
        structFormatFloatBrace = "{}" + self.structFormatFloat
        values = np.array(struct.unpack(structFormatFloatBrace.format(3*self.nPoints),
                                        data[: 3*self.nPoints * struct.calcsize(self.structFormatFloat)]))
        
        
        # Take the points data and organize into x, y, and z.
        self.points = np.zeros((self.nPoints,3))
        self.points[:,0] = values[0::3]
        self.points[:,1] = values[1::3]
        self.points[:,2] = values[2::3]

        #print(self.points[0:10,:])
        
        
        # Close the binary points file and delete unnecessary variables.
        fileID.close()
        del content, data, lines, values



    


    # Read in the faces file data.
    def readFaces(self,facesFile):
        print('    -Reading the faces file...')

        # Open the faces file (currently only coded to read binary)
        fileID = open(facesFile,'rb')
        content = fileID.read()

        
        
        # Split the data into lines and then find the number of faces.
        lines = content.split(b"\n")
        self.nFaces = int(lines[18])-1
        
        
        # Split the data where the actual faces begin using .split, and 
        # then keep only the non-header data (the [1:] does this), and 
        # then join the data back into a single entity with .join.
        #data = content.split(lines[18], 1)[1]
        #data = b"\n(".join(data.split(b"\n(")[1:])
        data = b"\n(".join(content.split(b"\n(")[1:])
        
        
        
        # The OpenFOAM binary faces is stored as a faceCompactList, which
        # begins with a list of indices that mark the start of a new face
        # in a second list that contains point numbers that make up each
        # face.  For example, the beginning of the first list might be something
        # like 0, 4, 8, and so on.  That means that the points that make up the 
        # first face start at index 0 and end at index 3 in the following
        # point number list.  The second face starts at index 4 and ends at 7.
        
        nb_numbers = self.nFaces + 1
        structFormatIntBrace = "{}" + self.structFormatInt
        pointsbyface = struct.unpack(structFormatIntBrace.format(nb_numbers),
                                     data[0: nb_numbers * struct.calcsize(self.structFormatInt)],)
        
        data = content.split(str.encode(str(pointsbyface[-1])))[1]
        data = b"\n(".join(data.split(b"\n(")[1:])
        
        self.faces = {}
        for i in range(self.nFaces):
            self.faces[i] = {}
            self.faces[i]["nPts"] = pointsbyface[i + 1] - pointsbyface[i]
            self.faces[i]["idPts"] = np.array(
                struct.unpack(
                    structFormatIntBrace.format(self.faces[i]["nPts"]),
                    data[pointsbyface[i]*struct.calcsize(self.structFormatInt):
                         pointsbyface[i + 1]*struct.calcsize(self.structFormatInt)],
                )
            )
        
        
        
        # Close the binary points file and delete unnecessary data.
        fileID.close()
        del content, data, lines, pointsbyface





    
    # Read the owner file and load owner cell information.
    def readOwners(self,ownersFile):
        print('    -Reading the owners file...')
        

        # Open the owner file (currently only coded to read binary)
        fileID = open(ownersFile,'rb')
        content = fileID.read()
        

        
        # Split the data into lines and then find the number of owners.
        lines = content.split(b"\n")
        self.nOwners = int(lines[19])


        
        # Read the indices of the cells that own each face
        data = b"\n(".join(content.split(b"\n(")[1:])
        structFormatIntBrace = "{}" + self.structFormatInt
        self.owners = np.array(
            struct.unpack(
                structFormatIntBrace.format(self.nOwners),
                data[: self.nOwners * struct.calcsize(self.structFormatInt)],
            )
        )


        # Set the number of grid cells which is equivalent to unique face owners.
        # (i.e., every grid cell must own a face, but not every grid cell must
        # be a neighbor to a face).
        self.nCells = len(np.unique(self.owners))


        # Close the binary points file and delete unnecessary data.
        fileID.close()
        del content, data, lines



    

    
    # Read the owner file and load neighbour cell information.
    def readNeighbours(self,neighboursFile):
        print('    -Reading the neighbours file...')
        

        # Open the neighbour file (currently only coded to read binary)
        fileID = open(neighboursFile,'rb')
        content = fileID.read()
        

        
        # Split the data into lines and then find the number of neighbours.
        lines = content.split(b"\n")
        self.nNeighbours = int(lines[19])


        
        # Read the indices of the cells that own each face
        data = b"\n(".join(content.split(b"\n(")[1:])
        structFormatIntBrace = "{}" + self.structFormatInt
        self.neighbours = np.array(
            struct.unpack(
                structFormatIntBrace.format(self.nNeighbours),
                data[: self.nNeighbours * struct.calcsize(self.structFormatInt)],
            )
        )



        # Close the binary points file and delete unnecessary data.
        fileID.close()
        del content, data, lines



    

    
    # Compute face centers and areas.
    # This is done following how OpenFOAM does it in: Foam::primitiveMesh::makeFaceCentresAndAreas
    def computeFaceCentersAndAreas(self):

        self.faceCenters = np.zeros((self.nFaces,3))
        self.faceAreas = np.zeros((self.nFaces,))
        self.faceNormals = np.zeros((self.nFaces,3))

        ii = 0
        for fi in range(self.nFaces):
            ii += 1

            # Write out a completion message.
            self.completionMessage('Computing cell face centers and areas',ii,self.nFaces,self.percentCompleteInterval)

            # Get number of vertices of this face.
            nPoints = self.faces[fi]["nPts"]

            # For efficiency, if the face is a triangle, directly compute the face center and area.
            # The center is just the average of the three points.
            # The area normal vector is 1/2 the cross product of vectors formed by two of the triangle edges.
            if (nPoints == 3):
                self.faceCenters[fi,:] = (1.0/3.0) * (self.points[self.faces[fi]["idPts"][0]] 
                                                    + self.points[self.faces[fi]["idPts"][1]]
                                                    + self.points[self.faces[fi]["idPts"][2]]);
                self.faceAreaNormals[fi,:] = 0.5 * (self.cross((self.points[self.faces[fi]["idPts"][1]]
                                                              - self.points[self.faces[fi]["idPts"][0]]),
                                                               (self.points[self.faces[fi]["idPts"][2]]
                                                              - self.points[self.faces[fi]["idPts"][0]])))
                self.faceAreas[fi] = self.mag(self.faceAreaNormals[fi])


            # If the face is not a triangle, then deal with it as follows:
            else:
                # Initialize the variables we need to zero and make proper size.
                approxCenter = np.zeros((3,))
                trueCenter = np.zeros((3,))
                
                area = 0.0
                center = np.zeros((3,))
                areaNormal = np.zeros((3,))
                areaCenter = np.zeros((3,))

                areaSum = 0.0
                centerSum = np.zeros((3,))
                areaNormalSum = np.zeros((3,))
                areaCenterSum = np.zeros((3,))

                # First find the approximate center of the face by averaging the vertex points.
                for pi in range(nPoints):
                    approxCenter += self.points[self.faces[fi]["idPts"][pi]]
                approxCenter /= nPoints


                # Now break the face into triangles where one vertex is the approximate center, and
                # the other two points are face vertex points.  Compute the area and center of each
                # of these triangles.  The center is the average of the three triangle vertices.  The
                # triangle face area normal vector is the cross product of two vectors formed from two
                # edges of the triangle (but choose them consistently for all triangles of the face).
                # The center of the face is the area-weighted average of each triangle's center.
                for pi in range(nPoints):
                    piNext = self.faces[fi]["idPts"][(pi+1) % nPoints]
                    thisPoint = self.points[self.faces[fi]["idPts"][pi]]
                    nextPoint = self.points[piNext]
                    center = (1.0/3.0)*(thisPoint + nextPoint + approxCenter)

                    areaNormal = 0.5*self.cross(-(nextPoint-thisPoint),(approxCenter-nextPoint))
                    area = self.mag(areaNormal)
                    areaCenter = area*center

                    areaSum += area
                    centerSum += center
                    areaNormalSum += areaNormal
                    areaCenterSum += areaCenter

                trueCenter = areaCenterSum/areaSum

                self.faceCenters[fi] = trueCenter
                self.faceAreas[fi] = areaSum
                self.faceNormals[fi] = areaNormalSum/areaSum
                    
                #if (fi == 0):
                #    print("face center: ", self.faceCenters[fi])
                #    print("face area: ", self.faceAreas[fi])
                #    print("face normal: ", self.faceNormals[fi])



    

    
    # Compute cell centers and volumes.
    # This is done following how OpenFOAM does it in: Foam::primitiveMesh::makeCellCentresAndVols
    def computeCellCentersAndVolumes(self):
        

        self.cellCenters = np.zeros((self.nCells,3))
        self.cellVolumes = np.zeros((self.nCells,))
        self.cellnFaces = np.zeros((self.nCells,))

        
        # Compute the approximate cell center for all cells as the average of face centers.  We do this by 
        # first looping over all "owned" faced and summing into the cell owned by that face.  Then we loop over
        # all "neighboured" faces and summing into the cell neighboured by that face.  Finally we loop over
        # all cells and divide the sum of face centers by number of face centers.
        approxCenter = np.zeros((self.nFaces,3))

        totalSums = 2.0*(self.nOwners + self.nNeighbours + self.nCells)
        ii = 0
        for fi in range(self.nOwners):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            ci = self.owners[fi]
            self.cellnFaces[ci] += 1
            approxCenter[ci] += self.faceCenters[fi]

        for fi in range(self.nNeighbours):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            ci = self.neighbours[fi]
            self.cellnFaces[ci] += 1
            approxCenter[ci] += self.faceCenters[fi]

        for ci in range(self.nCells):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            approxCenter[ci] /= self.cellnFaces[ci]


        # Compute the true cell center for all the cells.  Do this by breaking each volume into a pyramid.  Each
        # face is the base of a pyramid, and each pyramid's top is the approximate cell center.  Find the centroid
        # and volume of each pyramid and then the cell's center is the volume weighted average of the pyramid
        # centroids.  Pyramid volume is 1/3 base area times height.  We do this by dotting the base area normal 
        # vector with the vector between the base centroid and the pyramid top and then multiplying the result by
        # 1/3.  The pyramid centroid lies 1/4 the distance along the vector between the base centroid and the 
        # pyramid top.  Again, we have to loop by owned faces, then neighboured faces, and finally over all cells
        # to divide by volume to complete the volume-weighted averaging.
        trueCenter = np.zeros((self.nFaces,3))
        volumeSum = np.zeros((self.nFaces,))
        volumeCenterSum = np.zeros((self.nFaces,3))

        for fi in range(self.nOwners):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            ci = self.owners[fi]
            pyramidVolume = (1.0/3.0)*self.dot((self.faceAreas[fi] * self.faceNormals[fi]), 
                                                -(approxCenter[ci] - self.faceCenters[fi]))
            pyramidCenter = (3.0/4.0)*(self.faceCenters[fi]) + (1.0/4.0)*(approxCenter[ci])

            volumeSum[ci] += pyramidVolume
            volumeCenterSum[ci] += pyramidCenter * pyramidVolume
            

        for fi in range(self.nNeighbours):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            ci = self.neighbours[fi]
            pyramidVolume = (1.0/3.0)*self.dot((self.faceAreas[fi] * self.faceNormals[fi]), 
                                               (approxCenter[ci] - self.faceCenters[fi]))
            pyramidCenter = (3.0/4.0)*(self.faceCenters[fi]) + (1.0/4.0)*(approxCenter[ci])

            volumeSum[ci] += pyramidVolume
            volumeCenterSum[ci] += pyramidCenter * pyramidVolume
            

            
        for ci in range(self.nCells):
            ii += 1
            
            # Write out a completion message.
            self.completionMessage('Computing cell centers and volumes',ii,totalSums,self.percentCompleteInterval)
                
            self.cellVolumes[ci] = volumeSum[ci]
            self.cellCenters[ci] = volumeCenterSum[ci] / volumeSum[ci]
            

        #ci = 0
        #print("cell center: ", self.cellCenters[ci])
        #print("cell volume: ", self.cellVolumes[ci])



    

    
    # Read a solution field.
    def readField(self,timeDir,varName):
        

        # Open the data file (currently only coded to read binary)
        fileID = open(timeDir + '/' + varName,'rb')
        content = fileID.read()


        # Get the variable name stored in the file
        lines = content.split(b"\n")
        varNameRead = lines[13]
        varNameRead = varNameRead.split(b"object")[1]
        varNameRead = varNameRead.split(b";")[0]
        varNameRead = varNameRead.split(b" ")[-1]
        varNameRead = varNameRead.decode('ascii')
        print('    -Reading the',varNameRead,'file...')


        # Determine type of data: scalar, vector, symmTensor, tensor
        dataType = lines[19]
        dataType = dataType.split(b"<")[1]
        dataType = dataType.split(b">")[0]
        dataType = dataType.decode('ascii')
        print('    -Data type is',dataType,'\n')
        
        if (dataType == 'scalar'):
            nComp = 1
        elif (dataType == 'vector'):
            nComp = 3
        elif (dataType == 'symmTensor'):
            nComp = 6
        elif (dataType == 'tensor'):
            nComp = 9
        else:
            nComp = 1

        
        # Get the number of data points.
        n = int(lines[20])
        
        
        # Split the data where the actual vector data begins using .split, and
        # then keep only the non-header data (the [1:] does this), and 
        # then join the data back into a single entity with .join.
        data = b"\n(".join(content.split(b"\n(")[1:])
        
        
        # Use struct unpack to take the binary data and turn it into a list of
        # of doubles, and then put those into a numpy array.
        structFormatFloatBrace = "{}" + self.structFormatFloat
        values = np.array(struct.unpack(structFormatFloatBrace.format(nComp*n),
                                        data[: nComp*n * struct.calcsize(self.structFormatFloat)]))
        
        
        # Take the vector data and organize into x, y, and z components.
        if (nComp > 1):
            data = np.zeros((n,nComp))
            for i in range(nComp):
                data[:,i] = values[i::nComp]
        else:
            data = values
        
        
        # Close the binary points file and delete unnecessary variables.
        fileID.close()
        del content, lines, values

        
        return data, nComp




