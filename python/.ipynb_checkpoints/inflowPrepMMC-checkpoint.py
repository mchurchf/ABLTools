# inflowPrepMMC.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 25 August 2017
#
# This is a module to deal with preparing inflow for mesoscale-microscale 
# coupled cases over complex terrain.


import numpy as np










# Class to deal with the mesoscale-microscale inflow preparation.
class inflowPrepMMC:
    
    # Initialize the class.
    def __init__(self):
        self.name = []
        self.nPoints = []
        self.xyz = []
        
        
        
        
    # Read OpenFOAM boundaryData format points file.
    def readBoundaryDataPoints(self,pointsFileName):
        
        # Figure out header length.
        f = open(pointsFileName,'r')
        lineText = f.readline()
        i = 0
        while (lineText[0] != '('):
            i = i + 1
            lineText = f.readline()
            
        headerLength = i - 1
        f.close()
        
        
        # Read in data.
        f = open(pointsFileName,'r')
        for i in range(headerLength):
            f.readline()
            
        self.nPoints = int(f.readline())
        
        f.readline()
        self.xyz = np.zeros((self.nPoints,3))
        for i in range(self.nPoints):
            data = f.readline().split()
            x = float(data[0][1:])
            y = float(data[1]) 
            z = data[2]
            z = float(data[2][0:len(z)-1])
            self.xyz[i,0] = x
            self.xyz[i,1] = y
            self.xyz[i,2] = z

        f.close()
    



    
    # Write OpenFOAM boundaryData format points file.
    def writeBoundaryDataPoints(self,pointsFileName):
        
        f = open(pointsFileName,'w')
        
        f.write('\n')
        nPoints = len(self.xyz)
        nDim = 3
        
        f.write(str(nPoints) + '\n')
        f.write('(\n')
        
        for i in range(nPoints):
            f.write('(')
                
            for j in range(nDim):
                if (j < nDim-1):
                    f.write(str(self.xyz[i,j]) + ' ')
                else:
                    f.write(str(self.xyz[i,j]))
                
            f.write(')\n')
                
        
        f.write(')')
        
        f.close()
        
        
        
        

    # Write WRFtoFOAM format data file.
    def writeDataWRFtoFOAM(self,pointsFileName,UTMoffsetX,UTMoffsetY,UTMzone):
        
        f = open(pointsFileName,'w')
        
        nPoints = len(self.xyz)
        
        for i in range(nPoints):
            [lat,lon] = self.UTMtoLatLon(self.xyz[i,0]+UTMoffsetX,self.xyz[i,1]+UTMoffsetY,UTMzone)
            
            f.write(str(lat) + ' ' + str(lon) + ' ' + str(self.xyz[i,2]) + '\n')
            
        f.close()
        
        
        
        
    
    
    # Read OpenFOAM boundaryData format field file.
    def readBoundaryDataField(self,fieldFileName):
        
        # Figure out the header length.
        f = open(fieldFileName,'r')
        lineText = f.readline()
        if (len(lineText.strip()) > 0):
            lineText = lineText.strip()
        lineTextPrevious = lineText
        i = 0
        while (((lineText[0] != '(')) or
               ((lineText[0] == '(') and (lineTextPrevious[0] == '/'))):
            i = i + 1
            lineTextPrevious = lineText
            lineText = f.readline()
            if (len(lineText.strip()) > 0):
                lineText = lineText.strip()
                
            
        headerLength = i - 1
        f.close()
        

        # Get data sizes.
        f = open(fieldFileName,'r')
        for i in range(headerLength):
            f.readline()
        
        self.nPoints = int(f.readline())
        f.readline()

        
        data = f.readline().split()
        
        # Remove parentheses
        if (data[0] == '('):
            del(data[0])
        if (data[-1] == ')'):
            del(data[-1])
        if (data[0][0] == '('):
            data[0] = data[0][1:]
        if (data[-1][-1] == ')'):
            data[-1] = data[-1][:-1]
        
        nDim = len(data)

        field = np.zeros((self.nPoints,nDim))
        
        f.close()
        
        
        # Read through the header.
        f = open(fieldFileName,'r')
        for i in range(headerLength+2):
            f.readline()
        
        # Read in the data stripping out parentheses, if they exist, and
        # newline characters.
        for i in range(self.nPoints):
            data = f.readline().strip().strip('(').strip(')').split()
            for j in range(nDim):
                field[i,j] = float(data[j])
                    
        f.close()
                    
        return field     






    # Write OpenFOAM boundaryData format field file.
    def writeBoundaryDataField(self,fieldFileName,field):
        
        print(field)
        f = open(fieldFileName,'w')
        
        f.write('\n')
        nPoints = len(field)
        nDim = np.size(field[0])
        
        f.write(str(nPoints) + '\n')
        f.write('(\n')
        
        for i in range(nPoints):
            if (nDim == 1):
                f.write(str(field[i]) + '\n')
            
            else:
                f.write('(')
                
                for j in range(nDim):
                    if (j < nDim-1):
                        f.write(str(field[i,j]) + ' ')
                    else:
                        f.write(str(field[i,j]))
                
                f.write(')\n')
                
        
        f.write(')')
        
        f.close()
        
        
        
        
        
        

    # Convert foamToVTK VTK files to xyz numpy array.
    def foamToVTKtoXYZ(self,VTKfileName):
        
        # Open the foamToVTK generated boundary file.
        f = open(VTKfileName,'r')
        
        
        # Read through the header.
        f.readline()
        self.name = f.readline()
        f.readline()
        f.readline()
        sizing = f.readline().split()
        self.nPoints = int(sizing[1])
        
        
        # Make an array for the x, y, z data.
        data = np.zeros((3*self.nPoints))
        
        
        # Read in the data.
        totalReadLength = 0
        while totalReadLength < 3*self.nPoints:   
            d = f.readline().split()
            nd = len(d)
            for i in range(nd):
                data[totalReadLength + i] = float(d[i])
            
            totalReadLength = totalReadLength + nd
            
                    
        # Close the VTK file.
        f.close()
            
        
        # Make an array for the ordered x, y, z data and assign the data to it.
        self.xyz = np.zeros((self.nPoints,3))
        for j in range(3):
            self.xyz[:,j] = data[j::3]
            
            
            
            
            
            
    # Sort the xyz data.
    def sortXYZ(self,planeType,sortOrder=None,sortIndex=None):
        
        # Based on the plane type, set the sort order to follow OpenFoam standard
        # boundary sorting order.
        if (sortIndex == None):
            if (sortOrder == None):
                if planeType == 'xy':
                    sortOrder = [0,1]
                elif planeType == 'xz':
                    sortOrder = [0,2]
                elif planeType == 'yz':
                    sortOrder = [2,1]
                else:
                    print('Invalid specification of planeType:')
                    print('Need xy, xz, or yz.')
                    
            sortInd = np.lexsort((self.xyz[:,sortOrder[0]],self.xyz[:,sortOrder[1]]))
            
        else:
            sortInd = sortIndex
                
        
        
        self.xyz = self.xyz[sortInd]
        
        return sortInd
    
    
    
    # Given a point, find the surrounding point index
    def findXYZ(self,p):
        
        xyzRel = p - self.xyz
        disRel = np.sqrt(np.square(xyzRel[:,0])+np.square(xyzRel[:,1])+np.square(xyzRel[:,2]))
        sortInd = np.argsort(disRel)
        
        #print(sortInd)
        #print(self.xyz)
        #print(xyzRel)
        #print("disRel")
        #print(disRel)
        #print("disRel[sortInd]")
        #print(disRel[sortInd])
        
        surroundingPoints = self.xyz[sortInd[0:4]]
        sortInd4 = np.lexsort((surroundingPoints[:,0],surroundingPoints[:,2]))
        surroundingPoints = surroundingPoints[sortInd4]
        sortInd4 = sortInd[0:4][sortInd4]
        
        #print("sortInd4")
        #print(sortInd4)
    
        #print("surroundingPoints")
        #print(surroundingPoints)
        
        #print("xyz[sortInd[0:4]]")
        #print(self.xyz[sortInd4])
        
        w = np.zeros((4,))
        dx = surroundingPoints[1,0] - surroundingPoints[0,0]
        dy = surroundingPoints[3,2] - surroundingPoints[0,2]
        dxl = p[0] - surroundingPoints[0,0]
        dxr = surroundingPoints[1,0] - p[0]
        dyu = surroundingPoints[2,2] - p[2]
        dyb = p[2] - surroundingPoints[0,2]
        
        w[0] = dyu*dxr/(dx*dy)
        w[1] = dyu*dxl/(dx*dy)
        w[2] = dyb*dxr/(dx*dy)
        w[3] = dyb*dxl/(dx*dy)
        
        
        #print(dx,dy,dxl,dxr,dyu,dyb)
        #print(w)
        
        return surroundingPoints,sortInd4,w
        
        
      
        
            
    
    
    
    
    
    # Translate the xyz data.
    def translateXYZ(self,T):
        
        self.xyz = self.xyz + T
        
        
        
    
    # Rotate the xyz data.
    def rotateXYZ(self,R):
        
        for i in range(self.nPoints):
            self.xyz[i] = np.matmul(self.xyz[i],R)
        
    
    
    
    
    
    
    # Given x,y,z data, find the bottom surface.
    def findSurface(self,planeType,boundingBox,searchRadius):
        
        # Apply bounding box.
        ind = []
        for i in range(self.nPoints):
            x = self.xyz[i,0]
            y = self.xyz[i,1]
            z = self.xyz[i,2]
            if ((x > boundingBox[0,0]) and (x < boundingBox[0,1]) and 
                (y > boundingBox[1,0]) and (y < boundingBox[1,1]) and
                (z > boundingBox[2,0]) and (z < boundingBox[2,1])):
                ind.append(i)
                        
        keepInd = np.asarray(ind)
                
        xyzBound = self.xyz[keepInd]
        
        
        # Now find lower left point.
        if (planeType == 'yz'):
            ind = self.xyz[:,1].argmin()
            
        print(ind)
            
        return xyzBound
    
    
    
    
    
    
    # Driver for conversion from latitude-longitude to UTM coordinates.
    def LatLonToUTM(self,latitude,longitude,UTMzone=None):
        
        if isinstance(latitude,np.ndarray):
            
            # Get the dimensions of the coordinates passed into here.
            dims = latitude.shape
            ndims = len(dims)
            ni = dims[0]
            nj = 0
            nk = 0
            if (ndims > 1):
                nj = dims[1]
            if (ndims > 2):
                nk = dims[2]
                
            UTMx = np.zeros(dims)
            UTMy = np.zeros(dims)
            
            if UTMzone is None:
                UTMzoneArray = np.chararray(dims,itemsize=4)
                
            
            for i in range(ni):
                for j in range(nj):
                    for k in range(nk):
                        if UTMZone is None:
                            [UTMx[i,j,k],UTMy[i,j,k],UTMzone[i,j,k]] = self.LatLonToUTM_elem(latitude[i,j,k],longitude[i,j,k])
                        else:
                            [UTMx[i,j,k],UTMy[i,j,k],UTMzone[i,j,k]] = self.LatLonToUTM_elem(latitude[i,j,k],longitude[i,j,k],UTMzone[i,j,k])
         
        else:
            if UTMzone is None:
                [UTMx,UTMy,UTMzone] = self.LatLonToUTM_elem(latitude,longitude)
            else:
                [UTMx,UTMy,UTMzone] = self.LatLonToUTM_elem(latitude,longitude,UTMzone)
                
                
                
                
                
            
    # Convert from latitude-longitude to UTM coordinates.
    def LatLonToUTM_elem(self,latitude,longitude,UTMzone=None):
        
        # Compute the UTM zone if not given
        if UTMzone is None:
            z = int((longitude/6.0) + 31.0)
                        
            if latitude < -72.0:
                l = 'C'
            elif latitude < -64.0:
                l = 'D'
            elif latitude < -56.0:
                l = 'E'
            elif latitude < -48.0:
                l = 'F'
            elif latitude < -40.0:
                l = 'G'
            elif latitude < -32.0:
                l = 'H'
            elif latitude < -24.0:
                l = 'J'
            elif latitude < -16.0:
                l = 'K'
            elif latitude < -8.0:
                l = 'L'
            elif latitude < 0.0:
                l = 'M'
            elif latitude < 8.0:
                l = 'N'
            elif latitude < 16.0:
                l = 'P'
            elif latitude < 24.0:
                l = 'Q'
            elif latitude < 32.0:
                l = 'R'
            elif latitude < 40.0:
                l = 'S'
            elif latitude < 48.0:
                l = 'T'
            elif latitude < 56.0:
                l = 'U'
            elif latitude < 64.0:
                l = 'V'
            elif latitude < 72.0:
                l = 'W'
            else:
                l = 'X'
                            
            if z < 10:
                zStr = '0' + str(z)
            else:
                zStr = str(z)
            
            UTMzone = zStr + ' ' + l
                        
                        
        # Compute the UTM coordinates given the zone.               
        sa = 6378137.0
        sb = 6356752.314245
        
        e2 = (((sa**2) - (sb**2))**0.5) / sb
        e2sqr = e2**2
        c = (sa**2) / sb
        
        lat = latitude * (np.pi/180.0)
        lon = longitude * (np.pi/180.0)
                    
        z = float(UTMzone[0:2])
        S = ((6.0*z) - 183.0)
        deltaS = lon - (S * (np.pi/180.0))
                    
        a = np.cos(lat) * np.sin(deltaS);
        epsilon = 0.5 * np.log((1.0 +  a)/(1.0 - a))
        nu = np.arctan(np.tan(lat)/np.cos(deltaS)) - lat
        v = (c / ((1.0 + (e2sqr * (np.cos(lat))**2)))**0.5) * 0.9996
        ta = (0.5*e2sqr) * epsilon**2 * (np.cos(lat))**2
        a1 = np.sin(2.0*lat)
        a2 = a1 * (np.cos(lat))**2
        j2 = lat + (0.5*a1)
        j4 = ((3.0*j2) + a2) / 4.0
        j6 = ((5.0 * j4) + (a2 * (np.cos(lat))**2))/3.0
        alpha = (3.0/4.0) * e2sqr
        beta = (5.0/3.0) * alpha**2
        gamma = (35.0/27.0) * alpha**3
        Bm = 0.9996 * c * (lat - alpha * j2 + beta * j4 - gamma * j6)
        UTMx = epsilon * v * (1.0 + (ta/3.0)) + 500000.0
        UTMy = nu * v * (1.0 + ta) + Bm
        if (UTMy < 0.0):
            UTMy = UTMy + 9999999.0

        return (UTMx,UTMy,UTMzone)
    
    
    
    
    
    
    # Driver for conversion from UTM coordinates to latitude-longitude.
    def UTMtoLatLon(self,UTMx,UTMy,UTMzone):
    
        if isinstance(UTMx,np.ndarray):
            
            # Get the dimensions of the coordinates passed into here.
            dims = UTMx.shape
            ndims = len(dims)
            ni = dims[0]
            nj = 0
            nk = 0
            if (ndims > 1):
                nj = dims[1] 
            if (ndims > 2):
                nk = dims[2]
            
            latitude = np.zeros(dims)
            longitude = np.zeros(dims)
            
            for i in range(ni):
                for j in range(nj):
                    for k in range(nk):
                        [latitude[i,j,k],longitude[i,j,k]] = self.UTMtoLatLon_elem(self,UTMx[i,j,k],UTMy[i,j,k],UTMzone[i,j,k])
                        
        else:
            [latitude,longitude] = self.UTMtoLatLon_elem(UTMx,UTMy,UTMzone)
           
           
        return (latitude,longitude)
 

       
   
    
    
    # Convert from UTM coordinates to latitude-longitude.
    def UTMtoLatLon_elem(self,UTMx,UTMy,UTMzone):
        
        # Perform the conversion from UTM to latitude-longitude
        sa = 6378137.0
        sb = 6356752.314245
        
        e2 = (((sa**2) - (sb**2))**0.5) / sb
        e2sqr = e2**2
        c = (sa**2) / sb
                    
        # Deal with the UTM zone.
        zoneNumber = float(UTMzone[0:2])
                    
        if UTMzone[3] > 'M':
            y = UTMy
        else:
            y = UTMy - 10000000.0
                    
        sa = 6378137.0
        sb = 6356752.314245
        e2 = (((sa**2) - (sb**2))**0.5) / sb
        c = (sa**2) / sb
        x = UTMx - 500000.0
                    
        s = ((6.0*zoneNumber) - 183.0)
                    
        lat = y / (6366197.724*0.9996)
                    
        v = (c / ((1.0 + (e2sqr*(np.cos(lat))**2)))**0.5) * 0.9996
        a = x / v
        a1 = np.sin(2.0*lat)
        a2 = a1 * (np.cos(lat))**2
        j2 = lat + (0.5*a1)
        j4 = ((3.0*j2) + a2)/4.0
        j6 = ((5.0*j4) + (a2*(np.cos(lat))**2))/3.0
        alpha = (3.0/4.0)*e2sqr
        beta = (5.0/3.0)*(alpha**2)
        gamma = (35.0/27.0)*alpha**3
        Bm = 0.9996*c*(lat - alpha*j2 + beta*j4 - gamma*j6)
        b = (y - Bm) / v
        Epsi = ((e2sqr * a**2) / 2.0) * (np.cos(lat))**2
        Eps = a * (1.0 - (Epsi/3.0))
        nab = (b * (1.0 - Epsi)) + lat
        senoheps = (np.exp(Eps) - np.exp(-Eps)) / 2.0
        Delt = np.arctan(senoheps / (np.cos(nab) ) )
        TaO = np.arctan(np.cos(Delt) * np.tan(nab))
        longitude = (Delt *(180.0/np.pi)) + s
        latitude = (lat + (1.0 + e2sqr*(np.cos(lat)**2) - (3.0/2.0) * e2sqr * np.sin(lat) * np.cos(lat) * (TaO - lat)) * (TaO - lat)) * (180.0/np.pi)                    
        
        return (latitude,longitude)
        
        
        
        
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
