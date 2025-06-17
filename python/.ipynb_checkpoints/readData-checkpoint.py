# readData.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
#
# Created: 8 May 2017
# Last Updated: 17 June 2025
#
#
#
# This is a module that contains methods to read in data sampled from atmospheric LES.

# isStringInList
# getOutputTimes
# findNearest
# nekSlice
# amrwind_netcdf
# unstructuredVTK
# structuredVTK
# ensight
# readFluxesSOWFA
# uStarSOWFA
# qWallSOWFA
# planarAveragesSOWFA
# planarAveragesAMRWind
# planarAveragesAMRWind_ascii
# planarAveragesAMRWind_netcdf
# probeOpenFOAM
# ncdump





import numpy as np






       
    
# Test to see if a string is contained in a list of strings.
def isStringInList(str,strList):
    
    inList = False
    
    for i in range(len(strList)):
        if (strList[i] == str):
            inList = True
            
    return inList
    


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




def findNearest(array,value):
    array = np.asarray(array)
    idx = int((np.abs(array - value)).argmin())
    
    return idx, array[idx]



def findBoundingIndices(array,value):
    




def nekSlice(fileName):
    
    # Read the file.
    data = np.loadtxt(fileName,skiprows=3)
    print(data.shape)
    
    xyzList = data[:,0:1]
    
    sort_y = data[:,1].argsort(kind='stable')
    data = data[sort_y]
    sort_x = data[:,0].argsort(kind='stable')
    data = data[sort_x]
    
    
    nx = len(np.unique(data[:,0]))
    ny = len(np.unique(data[:,1]))
    print(nx,ny)
    
    xy = np.zeros((nx,ny,2))
    uvw = np.zeros((nx,ny,3))
    
    
    for j in range(ny):
        xy[:,j,:] = data[j*ny:(j+1)*ny,0:2]
        uvw[:,j,:] = data[j*ny:(j+1)*ny,2:]
    

    
    return xy,uvw
    
    
    
    
def amrwind_netcdf(fileName,samplerName):
    # Import necessary modules
    import netCDF4 as ncdf
    
    
    # Read the netcdf file.
    netcdfData = ncdf.Dataset(fileName, 'r')

    
    # Extract time.
    netcdfTime = netcdfData.variables['time']
    print('netcdfTime',netcdfTime)
    netcdfMaxTime = np.max(netcdfTime)
    print('netcdfMaxTime',netcdfMaxTime)
    print(netcdfTime[-10:])

    # Get the shape of the plane.
    planeShape = (netcdfData[samplerName].ijk_dims[2],
                  netcdfData[samplerName].ijk_dims[0],
                  netcdfData[samplerName].ijk_dims[1])

    nx = planeShape[1]
    ny = planeShape[2]
    nz = planeShape[0]
    nt = len(netcdfTime)
    

    x = netcdfData[samplerName]['coordinates'][:,0]
    y = netcdfData[samplerName]['coordinates'][:,1]
    z = netcdfData[samplerName]['coordinates'][:,2]
    
    x = x.reshape(planeShape)
    y = y.reshape(planeShape)
    z = z.reshape(planeShape)
    
    
    xyz = np.zeros((nx,ny,nz,3))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xyz[i,j,k,0] = x[k,i,j]
                xyz[i,j,k,1] = y[k,i,j]
                xyz[i,j,k,2] = z[k,i,j]


    field = []
    fieldName = []
    fieldSize = []
    
    # Get latest u,v,w
    aaa = netcdfData[samplerName]['velocityx'].shape
    print(aaa)
    u = netcdfData[samplerName]['velocityx'][-1] # Take the most recent value
    v = netcdfData[samplerName]['velocityy'][-1]
    w = netcdfData[samplerName]['velocityz'][-1] # not currently using
    theta = netcdfData[samplerName]['temperature'][-1] # not currently using

    u = u.reshape(planeShape)
    v = v.reshape(planeShape)
    w = w.reshape(planeShape) # not currently using
    theta = theta.reshape(planeShape) # not currently using
    
    fieldSize.append(3)
    fieldName.append('U')
    field_ = np.zeros((nx,ny,nz,3))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                field_[i,j,k,0] = u[k,i,j]
                field_[i,j,k,1] = v[k,i,j]
                field_[i,j,k,2] = w[k,i,j]
    field.append(field_)
                      
    fieldSize.append(1)
    fieldName.append('theta')
    field_ = np.zeros((nx,ny,nz,1))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                field_[i,j,k,0] = theta[k,i,j]
    field.append(field_)
                      
                      
    return xyz,field,fieldName,fieldSize,nx,ny,nz,nt



def unstructuredVTK(fileName):
    # Read the number of points.
    nHeader = 4
    dd = np.loadtxt(fileName,skiprows=nHeader,max_rows=1,dtype='str')
    nPoints = int(dd[1])
    
    
    # Read the points.
    nHeader = nHeader + 1
    points = np.loadtxt(fileName,skiprows=nHeader,max_rows=nPoints)

    
    # Read the number of cells.
    nHeader = nHeader + nPoints + 1
    dd = np.loadtxt(fileName,skiprows=nHeader,max_rows=1,dtype='str')
    nCells = int(dd[1])
    
    
    # Read the cells.
    nHeader = nHeader + 1
    cells = np.loadtxt(fileName,skiprows=nHeader,max_rows=nCells,dtype='int')
    
    # Compute cell centroids.
    xyzList = np.zeros((nCells,3))
    for i in range(nCells):
        for n in range(cells[i][0]):
            xyzList[i] = xyzList[i] + points[cells[i][n+1],:]
        xyzList[i] = (1.0/float(cells[i][0]))*xyzList[i]

    
    x = np.unique(xyzList[:,0])
    y = np.unique(xyzList[:,1])
    z = np.unique(xyzList[:,2])
    
    nx = len(x)
    ny = len(y)
    nz = len(z)
    
    
    sort_z = xyzList[:,2].argsort()
    xyzList = xyzList[sort_z]
    sort_y = xyzList[:,1].argsort(kind='stable')
    xyzList = xyzList[sort_y]
    sort_x = xyzList[:,0].argsort(kind='stable')
    xyzList = xyzList[sort_x]
    
    ii = 0
    xyz = np.zeros((nx,ny,nz,3))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xyz[i,j,k,:] = xyzList[ii,:]
                ii = ii + 1    
        
    
    # Read the number of fields.
    nHeader = nHeader + nCells + 1
    dd = np.loadtxt(fileName,skiprows=nHeader,max_rows=1,dtype='str')
    nFields = int(dd[2])
    
    
    # Read the fields.
    fieldName = []
    field = []
    fieldSize = []
    for m in range(nFields):
        nHeader = nHeader + 1
        dd = np.loadtxt(fileName,skiprows=nHeader,max_rows=1,dtype='str')
        fieldName.append(dd[0])
        fieldSize.append(int(dd[1]))
        
        nHeader = nHeader + 1
        fieldList = np.loadtxt(fileName,skiprows=nHeader,max_rows=nCells)
        
        fieldList = fieldList[sort_z]
        fieldList = fieldList[sort_y]
        fieldList = fieldList[sort_x]
        
        ii = 0
        field_ = np.zeros((nx,ny,nz,fieldSize[m]))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    field_[i,j,k,:] = fieldList[ii,:]
                    ii = ii + 1
        
        field.append(field_)
        
    
    
    return xyz,field,fieldName,fieldSize,nx,ny,nz
    
    
        





# Read in structured VTK data.  The data is returned as a list called fields that
# contains numpy arrays.  For example, the VTK file may contain both pressure
# and temperature, or velocity and vorticity, so each quantity will be contained
# in a separate array in the list.  Then for each array, the first index is
# the number of components, so for scalars, there is only one component, but for
# vectors there will be three.  The remaining dimensions are the number of 
# sample points in x, y, and z.

def structuredVTK(fileName):
    print('    -reading structured VTK data...')
  
  
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
    
    print('      -data dimensions: ',dims)
    print('      -data origin: ',origin)
    print('      -data spacing: ',spacing)
  
  
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
        print('      -reading field:',m+1,'of',nFields)
        if (m > 0):
            f.readline()
          
        d = f.readline()
        c = d.split()
        fieldName.append(c[0])
        fieldDim.append(int(c[1]))
        print('        -field name:',fieldName[m],'of dimensions:',fieldDim[m])
        dataArray = np.zeros((fieldDim[m], dims[0], dims[1], dims[2]))
        for k in range(dims[2]):
            progressBar = '';
            progress = int(20.0*(k/dims[2]))
            percentComplete = np.round(100.0*(k/dims[2]),0)
            for kk in range(20):
                if (kk < progress):
                    progressBar = progressBar + '='
                elif (kk == progress):
                    progressBar = progressBar + '>'
                else:
                    progressBar = progressBar + ' '
            progressBar = progressBar + '|'
                
            print('          -',percentComplete,'% complete ',progressBar, end='\r')
            for j in range(dims[1]):
                for i in range(dims[0]):
                    l = f.readline()
                    l = l.split()
                    for n in range(fieldDim[m]):
                        dataArray[n][i][j][k] = float(l[n])
                      
        field.append(dataArray)
        print(' ')
  
  
    # Close file.
    f.close()
  
  
    # Return the data.
    return dataSetName, dims, origin, spacing, x, y, z, nFields, fieldName, fieldDim, field







# Read in Ensight format data.  The data is returned as a list called fields that
# contains numpy arrays.  For example, the Ensight file may contain both pressure
# and temperature, or velocity and vorticity, so each quantity will be contained
# in a separate array in the list.  Then for each array, the first index is
# the number of components, so for scalars, there is only one component, but for
# vectors there will be three.  The remaining dimensions are the number of 
# sample points in x, y, and z.

def ensight(fileNameMesh,fileNameField,readFieldOnly):  
    x = []
    y = []
    z = []
    # Read in the geometry.
    if (readFieldOnly == 0):
        # Open the mesh file
        f = open(fileNameMesh,'r')
  
  
        # Get the data length.
        for i in range(8):
            f.readline()
        dims = int(f.readline())
        print(dims)
        #f.close()
  
  
        # Read in the x,y,z data.
        print('np.genfromtxt x,y,z data...')
        mesh = np.genfromtxt(fileNameMesh,skip_header=9,max_rows=dims*3)
        print(mesh.shape)
        #df = pd.read_csv(f,header=None,nrows=3*dims).values
        #mesh = pd.DataFrame(data=df.reshape((dims,3),order='F'), columns=['x','y','z']).to_numpy()
        x = mesh[0:dims]
        y = mesh[dims:(2*dims)]
        z = mesh[(2*dims):(3*dims)]
        
        print(x.shape)
      
      
      
      
  
    
    # Read in the field data.
    f = open(fileNameField,'r')
  
  
    # Get the data type
    fieldDim = 0
    dataType = f.readline().strip('\n')
  
    if (dataType == 'scalar'):
        fieldDim = 1
    elif (dataType == 'vector'):
        fieldDim = 3
    elif (dataType == 'tensor'):
        fieldDim = 9
        
    print(fieldDim)

    f.close()
    f = open(fileNameField,'r')
    for i in range(4):
        f.readline()
  
      
    # Read the field
    field = np.zeros((dims,fieldDim))
    print('np.genfromtxt field data...')
    data = np.genfromtxt(fileNameField,skip_header=4,max_rows=dims*fieldDim)
    #vals = pd.read_csv(f,header=None,nrows=fieldDim*dims).values
    #print(vals.shape)
    #field = pd.DataFrame(data=vals.reshape((dims,fieldDim),order='F')).to_numpy()
  
    for i in range(fieldDim):
        field[0:dims,i] = data[i*dims:(i+1)*dims]

  
  
    # Return the data.
    return dims, x, y, z, fieldDim, field




# Read in uStar data written from the SOWFA output file.
def readFluxesSOWFA(uStarDirectory):
    # Find the time directories.
    [nTimes,outputTimes] = getOutputTimes(uStarDirectory)
    
    print('   -reading from averaging directories ', outputTimes)
    
    # For each time directory, read the data and then concatenate.

    t = []
    dt = []
    data = []
    uStar = []
    qWall = []
    TWall = []
    for i in range(nTimes):
        fileName = uStarDirectory + '/' + outputTimes[i] + '/' + 'surfaceFluxHistory'
        
        dataInt = np.genfromtxt(fileName,skip_header=1)
        tInt = dataInt[:,0]
        dtInt = dataInt[:,1]
        dataInt = dataInt[:,2:]
        
        if (i < nTimes-1):
            dtLast = tInt[-1]-tInt[-2]
            tNext = float(outputTimes[i+1])
            index = np.argmax(np.abs(tInt - tNext) <= 1.5*dtLast)
        else:
            index = len(tInt)
            

        if (i == 0):
            t = tInt[0:index]
            dt = dtInt[0:index]
            data = dataInt[0:index,:]
        else:
            t = np.concatenate((t,tInt[0:index]),axis=0)
            dt = np.concatenate((dt,dtInt[0:index]),axis=0)
            data = np.concatenate((data,dataInt[0:index,:]),axis=0)
       
    # Look for time repeats.  Keep the second one
    keep = []
    lastTime = 1.0E15
    for i in range(len(t)-1,0,-1):
        if (t[i] < lastTime):
            lastTime = t[i]
            keep.append(i)
            
    keep = np.flipud(keep)
            
    t = t[keep]
    dt = dt[keep]
    data = data[keep]
    uStar = data[:,0]
   #Twall = data[:,1]
   #qWall = data[:,2]
        
        
    return t, dt, uStar, qWall, TWall





# Read in uStar data grep'd from the SOWFA output file.
def uStarSOWFA(uStarDirectory):
    # Find the time directories.
    [nTimes,outputTimes] = getOutputTimes(uStarDirectory)
    
    print('   -reading from averaging directories ', outputTimes)
    
    # For each time directory, read the data and then concatenate.

    uStar = []
    for i in range(nTimes):
        dataInt = []
        fileName = uStarDirectory + '/' + outputTimes[i] + '/' + 'uStar'
        f = open(fileName,'r')
        
        datatxt_ = f.readlines()
        print(len(datatxt_))
        for i in range(4):
            print(datatxt_[i])
        datatxt = datatxt_[1::3]
        print(len(datatxt))
        for i in range(4):
            print(datatxt[i])
        
        for line in datatxt:
            lineOfInterest = line.find('Avg iterations')
            if (lineOfInterest != -1):
                indStart = line.find('=') + 1
                indEnd = line.find('L') - 1
                dataInt.append(float(line[indStart:indEnd]))
            
        dataInt = np.array(dataInt)
        
        if (i == 0):
            uStar = dataInt
        else:
            uStar = np.concatenate((uStar,dataInt),axis=0)
    
    return uStar






# Read in surface temperature flux data grep'd from the SOWFA output file.
def qWallSOWFA(uStarDirectory):
    # Find the time directories.
    [nTimes,outputTimes] = getOutputTimes(uStarDirectory)
    
    print('   -reading from averaging directories ', outputTimes)
    
    # For each time directory, read the data and then concatenate.

    qWall = []
    for i in range(nTimes):
        dataInt = []
        fileName = uStarDirectory + '/' + outputTimes[i] + '/' + 'uStar'
        f = open(fileName,'r')
        
        datatxt_ = f.readlines()
        print(len(datatxt_))
        for i in range(4):
            print(datatxt_[i])
        datatxt = datatxt_[0::3]
        print(len(datatxt))
        for i in range(4):
            print(datatxt[i])
        
        for line in datatxt:
            lineOfInterest = line.find('qwMean')
            if (lineOfInterest != -1):
                indStart = line.find('qwMean') + 8
                indEnd = line.find('LMean') - 1
                dataInt.append(float(line[indStart:indEnd]))
            
        dataInt = np.array(dataInt)
        
        if (i == 0):
            qWall = dataInt
        else:
            qWall = np.concatenate((qWall,dataInt),axis=0)
    
    return qWall
        

    
    
    
    
    
    

# Read in the planar averaged data output by SOWFA.

def planarAveragesSOWFA(averagingDirectory,varName,version):
    # Find the time directories.
    [nTimes,outputTimes] = getOutputTimes(averagingDirectory)
    
    print('   -reading from averaging directories ', outputTimes)

    
    # Set the data file name and header length.
    fileName = averagingDirectory + '/' + outputTimes[0] + '/' + varName
    nHeader = 1

    
    if ((version == 'new') or (version == '6')):
        print("   -reading data from SOWFA-6")
        
        # Open the data file to gather metadata.
        f = open(fileName,'r')
        
        
        # The new version has a header to the file.
        nHeader = 3


        # Read the heights of data points.
        data = (f.readline()[12:].split())
        nz = len(data)
        z = np.zeros((nz,))
        for i in range(nz):
            z[i] = float(data[i])
        
        # Close the file.
        f.close()
        
    else:
        print("   -reading data from SOWFA-2.4.x")
        
        # Read the heights file.
        z = np.genfromtxt(averagingDirectory + '/' + outputTimes[0] + '/hLevelsCell')    
        
  
    # Open the data file to gather metadata.
    f = open(fileName,'r')
    
    
    # Find out the dimension of the data (i.e., scalar, vector, tensor)
    for i in range(nHeader):
        data = f.readline()
        
    iStart = data.find('(')
    iEnd = data.find(')')
    dim = 1
    if (iStart != -1):
        data = data[iStart+1:iEnd].split()
        dim = len(data)
   

    # Close the file.
    f.close()

    
    # For each time directory, read the data and then concatenate.
    tInt = []
    dtInt = []
    dataInt = []
    t = []
    dt = []
    data = []
    for i in range(nTimes):
        fileName = averagingDirectory + '/' + outputTimes[i] + '/' + varName
        
        if (version == 'new'):
            f = open(fileName,'r')
            datatxt = f.readlines()[2:]
            f.close()
        else:
            datatxt = np.genfromtxt(averagingDirectory + '/' + outputTimes[i] + '/' + varName)

        datafloat = []
        for line in datatxt:
            line = [ float(val) for val in
                    line.replace('(','').replace(')','').split() ]
            datafloat.append(line)
        datafloat = np.array(datafloat)

        tInt = datafloat[:,0]
        dtInt = datafloat[:,1]
        dataInt = np.zeros((len(tInt),nz,dim))
        for j in range(dim):
            dataInt[:,:,j] = datafloat[:,2+j::dim]

        if (i < nTimes-1):
            dtLast = tInt[-1]-tInt[-2]
            tNext = float(outputTimes[i+1])
            index = np.argmax(np.abs(tInt - tNext) <= 1.5*dtLast)
        else:
            index = len(tInt)

        if (i == 0):
            t = tInt[0:index]
            dt = dtInt[0:index]
            data = dataInt[0:index,:]
        else:
            t = np.concatenate((t,tInt[0:index]),axis=0)
            dt = np.concatenate((dt,dtInt[0:index]),axis=0)
            data = np.concatenate((data,dataInt[0:index,:]),axis=0)


    return z, t, dt, data, nz, dim






# Read in the planar averaged data output by AMR-Wind.  (This just calls
# the specific function to read netcdf or ascii data, but assumes netcdf
# if not otherwise specified).
def planarAveragesAMRWind(averagingDirectory,format='netcdf'):
    z = []
    t = []
    data = []
    profileData = []
    if ((format == 'netcdf') or (format == 'NETCDF') or (format == 'netCDF') or (format == 'NetCDF')):
        z, t, data, profileData = planarAveragesAMRWind_netcdf(averagingDirectory)
    elif ((format == 'ascii') or (format == 'ascii') or (format == 'ascii') or (format == 'ascii')):
        z, t, data, profileData = planarAveragesAMRWind_ascii(averagingDirectory)
        
    return z, t, data, profileData    


        

    
# Read in the planar averaged data output by AMR-Wind in ASCII format.     
def planarAveragesAMRWind_ascii(averagingDirectory):
    import os
    import time
    #import matplotlib.pyplot as plt
    
    
    print('Using the AMRWind ABLStats ASCII reader...')
    
    # Currently because of precision, AMR-Wind's output files lose precision for times greater than 100,000 s,
    # so we have to do some special guessing about precision in that case.  This should be fixed in AMR-Wind,
    # and then this reader has to be updated to reflect that.
    precisionLossThreshold = 100000.0

    z = []
    profileData = []
    
    tic = time.perf_counter()
    
    # Figure out how many abl statistics scalar data (not profiles) files there are.
    files = []
    files_ = (os.listdir(averagingDirectory))
    for i in range(len(files_)):
        if ('abl_statistics' in files_[i]):
            files.append(files_[i])
            
    nFiles = len(files)
    print('   -Files read: ',files)
    
    
    # Get the starting time of each scalar data file.
    tStart = []
    for i in range(nFiles):
        asciidata = np.loadtxt(averagingDirectory + '/' + files[i],delimiter=',',skiprows=1)
        tStart.append(asciidata[0,0])
        del(asciidata)
        
    sortedInd = np.argsort(np.asarray(tStart))    
    
    
    # Read in the scalar data files.
    for i in range(len(sortedInd)):
        j = sortedInd[i]
        asciidata = np.loadtxt(averagingDirectory + '/' + files[j],delimiter=',',skiprows=1)

        t_ = np.array(asciidata[:,0])
        
        t0 = t_[0]
        if (i == 0):
            dt = t_[1]-t_[0]
        #print(t0,t_[-1],dt,t_)
        
        # This is a fix for the lack of precision for large numbers in the output files. 
        for m in range(len(t_)):
            if (t0 >= precisionLossThreshold):
                t_[m] = t0 + (m+1)*dt
            else:
                t_[m] = t0 + m*dt
        #print(t0,t_[-1],dt,t_)
        
        data_ = np.zeros((len(t_),8))
        
        for m in range(8):
            data_[:,m] = asciidata[:,m+1]
        
        if (i < nFiles-1):
            dtLast = t_[-1]-t_[-2]
            
            tNext = tStart[sortedInd[i+1]]
            if (tNext >= precisionLossThreshold):
                tNext = tNext + dtLast
                
            #print('tNext = ',tNext)
            index = np.argmax(np.abs(t_ - tNext) <= 1.0*dtLast)+1
        else:
            index = len(t_)
            
        #print(index,len(t_))
            

        if (i == 0):
            t = t_[0:index]
            data = data_[0:index,:]
            #profileData = profileData_[0:index,:,:]
        else:
            t = np.concatenate((t,t_[0:index]),axis=0)
            data = np.concatenate((data,data_[0:index,:]),axis=0)
            #profileData = np.concatenate((profileData,profileData_[0:index,:,:]),axis=0)
    
        del(asciidata)
        
    del(t_,data_)
    
    
    # Get the profile data heights (both the level 0 and finest level).
    #  level 0
    asciidata = np.loadtxt(averagingDirectory + '/' + 'plane_average_temperature.txt',delimiter=',')
    dataDims = asciidata.shape
    i = 0
    while (asciidata[i+1,0] == asciidata[i,0]):
        i = i + 1
    nzLevel0 = i + 1
    zLevel0 = asciidata[0:nzLevel0,2]
    del(asciidata)
    
    #  finest level
    asciidata = np.loadtxt(averagingDirectory + '/' + 'plane_average_temperature_fine.txt',delimiter=',')
    dataDims = asciidata.shape
    i = 0
    while (asciidata[i+1,0] == asciidata[i,0]):
        i = i + 1
    nzLevelFine = i + 1
    zLevelFine = asciidata[0:nzLevelFine,2]
    
    
    # Get the time information in the profiles.  Note, that since the profile files concatenate, data
    # can get repeated, so look for places this happens and eliminate them making time completely monotonic.
    # Do this from the fine data because it has higher precision (this really needs to be fixed in AMR-Wind).
    nt_ = int(dataDims[0]/nzLevelFine)
    t_ = asciidata[::nzLevelFine,1]
    
    #  for some reason, the profile time stamps are offset from the scalar data by one time step, so try to figure
    #  the time step size and offset (again this needs to be fixed in AMR-Wind).
    t_offset = t[0] - t_[0]
    t_ = t_ + t_offset
    #print(len(t_),len(t))
    #plt.figure(0)
    #plt.plot(t_,'k-')
    
    # This is a fix for the lack of precision for large numbers in the output files. 
    #for m in range(len(t_)):
    #    t_[m] = t0 + m*dt
    #print(t0,dt,t_)
    
    #  because of file concatenation, data can get repeated.  Identify places where time is not monotonic.
    tLastOkay = t_[0]
    indEliminate = []
    for i in range(nt_-1):
        if (t_[i+1] > tLastOkay):
            tLastOkay = t_[i+1]
        else:
            indEliminate.append(i+1)
                
    t_ = np.delete(t_,indEliminate)

    nt_ = len(t_)
    
    del(asciidata)
    
    #plt.figure(1)
    #plt.plot(t,'k-')
    #plt.plot(t_,'r--')
    
    
    
    
    profileDataLevelFine_ = np.zeros((nt_,nzLevelFine,24))
    
    
    
    # Get the actual profile data.
    dataFiles = ['plane_average_velocity_fine.txt',
                 'plane_average_temperature_fine.txt',
                 'plane_average_velocity_mueff.txt',
                 'second_moment_temperature_velocity.txt',
                 'second_moment_velocity_velocity.txt',
                 'third_moment_velocity_velocity_velocity.txt']
    
    componentsToLoad = [[0, 1, 2],
                        [0],
                        [0],
                        [0, 1, 2],
                        [0, 4, 8, 1, 2, 5],
                        [0, 1, 2]]
    
    derivedComponents = [1,
                         0,
                         0,
                         0,
                         0,
                         0]
   
    j = 0
    for i in range(len(dataFiles)):
        
        # Load the text file data.
        asciidata = np.loadtxt(averagingDirectory + '/' + dataFiles[i], delimiter=',')
        
        # Get the raw shape of the data.
        nComponents = asciidata.shape[1] - 3
        nTimes = int(asciidata.shape[0]/nzLevelFine)
        
        # Reshape the data to easily access time, height, and component.
        #print('shape = ',asciidata.shape)
        asciidata = np.reshape(asciidata,(nTimes,nzLevelFine,nComponents+3))
        #print('shape = ',asciidata.shape)
        
        # Delete the previously found repeated data due to restarts and file concatenation.
        asciidata = np.delete(asciidata,indEliminate,axis=0)
        #print('shape = ',asciidata.shape)
        
        # Load the data into the variable profileData
        for m in componentsToLoad[i]:
            #print('Loading profile data index ',j)
            # This is for all non-derived data straight from the data files.
            profileDataLevelFine_[:,:,j] = asciidata[:,:,m+3]
            j = j + 1
            
        for m in range(derivedComponents[i]):
            #print('Loading profile data index ',j)
            if (i == 0 and m == 0):
                # This is the derived quantity horizontal wind speed.
                profileDataLevelFine_[:,:,j] = np.sqrt(np.square(profileDataLevelFine_[:,:,0]) +
                                                       np.square(profileDataLevelFine_[:,:,1]))
                j = j + 1
            
        
        
        del(asciidata)
    profileData = profileDataLevelFine_
        

        
     
    toc = time.perf_counter()
    print(f"   -Read data in {toc - tic:0.4f} s...")
    
    
    # For now, set the z to the fine level.  One day we want to blend the two.
    z = zLevelFine
    
    
    # Return the data.
    return z, t, data, profileData       

    

    
# Read in the planar averaged data output by AMR-Wind in NetCDF format.
def planarAveragesAMRWind_netcdf(averagingDirectory):
    # Import the necessary modules.
    import netCDF4 as ncdf
    import os
    import time
    
    print('Using the AMRWind ABLStats NetCDF reader...')
  
    tic = time.perf_counter()
    
    
    # Figure out how many abl statistics files there are.
    files = []
    files_ = (os.listdir(averagingDirectory))
    for i in range(len(files_)):
        if ('abl_statistics' in files_[i]):
            files.append(files_[i])
            
    nFiles = len(files)
    print('   -Files read: ',files)
    
    
    # Get the starting time of each file's data.
    tStart = []
    for i in range(nFiles):
        ncData = ncdf.Dataset(averagingDirectory + '/' + files[i])
        tStart.append(np.array(ncData["time"][:])[0])
        del(ncData)
        
    sortedInd = np.argsort(np.asarray(tStart))    
    
    
    # Read in the data file
    for i in range(len(sortedInd)):

        
        # Lists of variable names to be read.
        bulkDataVarNames = ["Q","Tsurf","ustar","wstar","L","zi",
                            "abl_forcing_x","abl_forcing_y"]
        
        profileDataVarNames = ["u","v","w","hvelmag","theta","mueff",
                               "theta'theta'_r","u'theta'_r","v'theta'_r","w'theta'_r",
                               "u'u'_r","v'v'_r","w'w'_r","u'v'_r","u'w'_r","v'w'_r",
                               "u'u'u'_r","v'v'v'_r","w'w'w'_r",
                               "u'theta'_sfs","v'theta'_sfs","w'theta'_sfs",
                               "u'v'_sfs","u'w'_sfs","v'w'_sfs"]

        
        # Read the NetCDF file.
        j = sortedInd[i]
        ncData = ncdf.Dataset(averagingDirectory + '/' + files[j])

        # Get time and height out of the data.
        t_ = np.array(ncData["time"][:])
        z = np.array(ncData["mean_profiles"]["h"][:])

        # Read in the bulk data variable names in the NetCDF data and store in a list.
        bulkDataVarNamesRead_ = ncData.variables.keys()
        bulkDataVarNamesRead = []
        for key in bulkDataVarNamesRead_:
            bulkDataVarNamesRead.append(key)

        # Read in the profile data variable names in the NetCDF data and store in a list.
        profileDataVarNamesRead_ = ncData["mean_profiles"].variables.keys()
        profileDataVarNamesRead = []
        for key in profileDataVarNamesRead_:
            profileDataVarNamesRead.append(key)



        # Read in all the bulk data variables in the variable list, if they exist.
        bulkData_ = np.zeros((len(t_),len(bulkDataVarNames)))

        for j in range(len(bulkDataVarNames)):
            if isStringInList(bulkDataVarNames[j],bulkDataVarNamesRead):
                bulkData_[:,j] = ncData[bulkDataVarNames[j]][:]


        # Read in all the profile data variables in the variable list.
        profileData_ = np.zeros((len(t_),len(z),len(profileDataVarNames)))
        
        for j in range(len(profileDataVarNames)):
            if isStringInList(profileDataVarNames[j],profileDataVarNamesRead):
                profileData_[:,:,j]  = ncData["mean_profiles"][profileDataVarNames[j]][:]

        
        if (i < nFiles-1):
            dtLast = t_[-1]-t_[-2]
            tNext = tStart[sortedInd[i+1]]
            index = np.argmax(np.abs(t_ - tNext) <= 1.5*dtLast)
        else:
            index = len(t_)

        if (i == 0):
            t = t_[0:index]
            bulkData = bulkData_[0:index,:]
            profileData = profileData_[0:index,:,:]
        else:
            t = np.concatenate((t,t_[0:index]),axis=0)
            bulkData = np.concatenate((bulkData,bulkData_[0:index,:]),axis=0)
            profileData = np.concatenate((profileData,profileData_[0:index,:,:]),axis=0)
    
    
    del(t_,bulkData_,profileData_)
    
    toc = time.perf_counter()
    print(f"   -Read data in {toc - tic:0.4f} s...")
    
    return z, t, bulkData, profileData
    
    
    
    
    
    
    
    
# Read in the OpenFOAM probe output data.

def probeOpenFOAM(probeDirectory,varName):
    # Find the time directories.
    [nTimes,outputTimes] = getOutputTimes(probeDirectory)
    
    
    
    # For each time directory, read the data and concatenate.
    t = []
    data = []
    
    probePosition = []
    nProbes = 0
    for i in range(nTimes):
        tInt = []
        dataInt = []
        
        fileName = probeDirectory + '/' + outputTimes[i] + '/' + varName
        
        # Handle the file header only if it is the first file in the time
        # sequence.
        if (i == 0):
            f = open(fileName,'r')
            inHeader = True
            while (inHeader):
                a = f.readline().split()
                firstChar = a[0]
                probeNum = int(a[2])
                if ((firstChar == '#') and (probeNum == nProbes)):
                    inHeader = True
                    probeX = float(a[3][1:])
                    probeY = float(a[4])
                    probeZ = float(a[5][:-1])
                    probePosition.append(np.zeros(3))
                    probePosition[nProbes][0] = probeX
                    probePosition[nProbes][1] = probeY
                    probePosition[nProbes][2] = probeZ
                    nProbes = nProbes + 1
                else:
                    inHeader = False
                    f.readline()
                    
                
        # Otherwise, skip through the header.
        else:
            f = open(fileName,'r')
            for j in range(nProbes+2):
                a = f.readline()
                
                
        # Find out the data type.
        a = f.readline().split()
        nComponents = 0
        dataType = []
        if (a[1][0] == '('):
            j = 1
            while (a[j][-1] != ')'):
               j = j + 1
            nComponents = j
            
        else:
            nComponents = 1
            
        if (nComponents == 1):
            dataType = 'scalar'
        elif (nComponents == 3):
            dataType = 'vector'
        elif (nComponents == 6):
            dataType = 'symmTensor'
        elif (nComponents == 9):
            dataType = 'tensor'
            
        
        m = 0
        while (a != []):
            tInt.append(float(a[0]))
            dataInt.append([None]*(nProbes))
            if (nComponents == 1):
                for j in range(nProbes):
                    dataInt[m][j] = float(a[j+1])
                    
            else:
                for j in range(nProbes):
                    val = np.zeros((nComponents))
                    for k in range(nComponents):
                        if (k == 0):
                            val[k] = float(a[j*nComponents + k + 1][1:])
                        elif (k==nComponents-1):
                            val[k] = float(a[j*nComponents + k + 1][:-1])
                        else:
                            val[k] = float(a[j*nComponents + k + 1])
                    
                    dataInt[m][j] = val
                           
            m = m + 1
            
            a = f.readline().split()
            
                
        # Close the file.
        f.close()
        

        tInt = np.array(tInt)
        if (i < nTimes-1):
            tNext = float(outputTimes[i+1])
            index = int((np.abs(tInt - tNext)).argmin())
            #index = np.argmax(np.abs(tInt >= tNext))
            if (index == 0):
                index = len(tInt)
        else:
            index = len(tInt)
            

        
        
        if (i == 0):
            t = tInt[0:index]
            data = dataInt[0:index]
        else:
            t = np.concatenate((t,tInt[0:index]),axis=0)
            data = np.concatenate((data,dataInt[0:index]),axis=0)
    
    
    # Convert lists to numpy arrays
    t = np.asarray(t)
    data = np.asarray(data)
    probePosition = np.asarray(probePosition)
        
    return t, data, probePosition, nComponents, nProbes












# A function to read NetCFD file information.

def ncdump(nc_fid, verb=True):
    
    from netCDF4 import Dataset
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars
