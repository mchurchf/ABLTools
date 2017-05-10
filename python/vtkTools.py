#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 09:59:42 2017

"""
import numpy as np 
import os, glob, sys
#==============================================================================
# 
#==============================================================================
def _read_vtkStructured_oneFile(vtkPath,verbose=False):
    """
    Reads in a VTK structured file.
    
    Parameters
    ----------
    vtkPath : str,
        full path to VTK file .OR. to a directory full of time sub-directories with vtk files
    verbose : bool,
        whether to print out metadata

    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw
    """    
    f = open(vtkPath)
    lines = f.readlines()
    f.close()
    
    nx, ny, nz = [ float(x) for x in lines[4].lstrip().split()[1:] ]
    xo, yo, zo = [ float(x) for x in lines[5].lstrip().split()[1:] ]
    dx, dy, dz = [ float(x) for x in lines[6].lstrip().split()[1:] ]
    npts = int(lines[7].split(" ")[1])
    
    x1d = [xo] if dx==0 else [ xo+i*dx for i in range(int(nx)) ]
    y1d = [yo] if dy==0 else [ yo+i*dy for i in range(int(ny)) ]
    z1d = np.arange(zo,zo+dz*nz,dz)           
    
    [Y,X,Z] = np.meshgrid(y1d,x1d,z1d)
    
    U = np.zeros(X.shape)
    V = np.zeros(X.shape)
    W = np.zeros(X.shape)
    
    assert(nx*ny*nz==npts)
    
    # find row index of first numeric value
    for iline,line in enumerate(lines):
        val = lines[iline].lstrip().rstrip().split()[0]    
        try:
            val = float(val)
            if isinstance(val,float):
                row = iline
                break
        except:
            1
    
    # recall that x varies first so this loop must be z->y->x            
    for iz,z in enumerate(z1d):
        for iy,y in enumerate(y1d):
            for ix,x in enumerate(x1d):
                u, v, w = [ float(x) for x in lines[row].lstrip().rstrip().split() ]
                U[ix,iy,iz] = u ; V[ix,iy,iz] = v ; W[ix,iy,iz] = w ; 
                row += 1

    data = [X,Y,Z,U,V,W]                
    meta = {}
    meta['dx'] = dx ; meta['dy'] = dy ; meta['dz'] = dz
    meta['nx'] = nx ; meta['ny'] = ny ; meta['nz'] = nz
    meta['xOrigin'] = xo ; meta['yOrigin'] = yo ; meta['zOrigin'] = zo
    meta['nPts'] = npts

    if verbose:            
        print "dx = {0}".format(dx)
        print "dy = {0}".format(dy)
        print "dz = {0}".format(dz)
        print "nx = {0}".format(nx)
        print "ny = {0}".format(ny)
        print "nz = {0}".format(nz)
        
    return data, meta
#==============================================================================
# 
#==============================================================================
def _read_vtkStructured_manyFiles(vtkPath,verbose=False):
    """
    Read in several *.vtk (VTK Structured) files.
    
    Parameters
    ----------    
    vtkPath : str,
        absolute path to directory where each time subdirectory is (and within each, a vtk structured file)
    t0 : float,
        starting time matching the name of the directory in which the *.xy file will be found
    dt : float,
        time increment in between files/directories
    nt : int,
        number of times to process
        
    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw
    """    
    pathNow = os.path.abspath(os.curdir)
    os.chdir(vtkPath)
    inp = read_inp()

    for itime,time in enumerate(inp['times']):

        timePath        = os.path.abspath(os.path.join(vtkPath,"{0:.3f}".format(time)))
        vtkFile         = glob.glob(os.path.join(timePath,'array*U*.vtk'))[0]
        data, meta      = _read_vtkStructured_oneFile(vtkFile,verbose=False)
        [X,Y,Z,U,V,W]   = data
        
        if verbose:
            print "Reading in {0}...".format(vtkFile)
    
        if itime==0:
            
            x4d = X.copy() ; y4d = Y.copy() ; z4d = Z.copy()
            u4d = U.copy() ; v4d = V.copy() ; w4d = W.copy()
    
            x4d = np.expand_dims(x4d, axis=0)
            y4d = np.expand_dims(y4d, axis=0)
            z4d = np.expand_dims(z4d, axis=0)
            u4d = np.expand_dims(u4d, axis=0)
            v4d = np.expand_dims(v4d, axis=0)
            w4d = np.expand_dims(w4d, axis=0)        
                
        else:
            
            x4d = np.append(x4d,[X],axis=0)
            y4d = np.append(y4d,[Y],axis=0)
            z4d = np.append(z4d,[Z],axis=0)
            u4d = np.append(u4d,[U],axis=0)
            v4d = np.append(v4d,[V],axis=0)
            w4d = np.append(w4d,[W],axis=0)
        
    os.chdir(pathNow)

    return [x4d,y4d,z4d,u4d,v4d,w4d], meta
#==============================================================================
# 
#==============================================================================
def read_vtkStructured(vtkPath,verbose=False):
    """
    Reads in VTK structured data (either one file or a set of files).

    Parameters
    ----------
    vtkPath : str,
        full path to VTK file .OR. to a directory full of time sub-directories with vtk files
    verbose : bool,
        whether to print out metadata

    Returns
    -------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
        
    @author: pdoubraw    
    """
    extension = os.path.splitext(vtkPath)[-1]
    if extension==".vtk":
        data, meta = _read_vtkStructured_oneFile(vtkPath=vtkPath,verbose=verbose)
    else:
        data, meta = _read_vtkStructured_manyFiles(vtkPath=vtkPath,verbose=verbose)
    return data, meta
#==============================================================================
# 
#==============================================================================
def write_vtkStructured(data,meta,fileOutPath,descStr="PLACEHOLDER",verbose=False):
    """
    Writes data in vtk structured format.
    
    Parameters
    ----------
    data : list,
        [X,Y,Z,U,V,W] each of which is a 3d array
    meta : dict,
        info about the grid
    fileOutPath : str,
        absolute path to vtk file you want to write      
    descStr : str,
        some header string describing what these data are
    """      
    f = open("top.txt", 'w')
    f.write('# vtk DataFile Version 3.0\n')  
    f.write('{0}\n'.format(descStr))  
    f.write('ASCII\n')  
    f.write('DATASET STRUCTURED_POINTS\n')  
    f.write('DIMENSIONS {0:d} {1:d} {2:d}\n'.format(int(meta['nx']),
            int(meta['ny']),int(meta['nz'])))  
    f.write('ORIGIN {0:.1f} {1:.1f} {2:.1f}\n'.format(meta['xOrigin'],meta['yOrigin'],meta['zOrigin']))  
    f.write('SPACING {0:.1f} {1:.1f} {2:.1f}\n'.format(meta['dx'],meta['dy'],meta['dz']))  
    f.write('POINT_DATA {0:d}\n'.format(meta['nPts']))  
    f.write('VECTORS vAmb float\n')  
    f.close()    
    
    [X,Y,Z,U,V,W]   = data    
    U = np.ravel(U, order='F') ; V = np.ravel(V, order='F') ; W = np.ravel(W, order='F')        
    data = np.zeros((len(U),3))
    data[:,0] = U ; data[:,1] = V ; data[:,2] = W     
    np.savetxt("bot.txt",data)
    os.system("cat top.txt bot.txt > {0}".format(fileOutPath))
    os.remove("top.txt")
    os.remove("bot.txt")    

    if verbose:
        print "Saved data to {0}".format(fileOutPath)
    return
#==============================================================================
# 
#==============================================================================
def structuredVTK(fileName):
    """
    Read in structured VTK data.  The data is returned as a list called fields that contains numpy arrays.  For example, the VTK file may contain both pressure and temperature, or velocity and vorticity, so each quantity will be contained in a separate array in the list.  Then for each array, the first index is the number of components, so for scalars, there is only one component, but for vectors there will be three.  The remaining dimensions are the number of sample points in x, y, and z.

    @author: mchurchf
    """
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
#==============================================================================
# 
#==============================================================================
def write_inp(t0,dt,nt,zref=99.0,xyname="array.1_U.xy",btsPrefix="test",theta=0.0,workingPath=None):
    """
    Write the *.inp file that is necessary to convert *.xy to *.bts using of2fast (Sang's fortran code).

    "tsConv.inp" entries
    
    0                   : theta
    107                 : numgrid_y (see TurbSim manual to understand this)
    107                 : numgrid_z (see TurbSim manual to understand this)
    1.06                : gridres_z (see TurbSim manual to understand this) 
    1.5                 : gridres_y (see TurbSim manual to understand this)
    1                   : timestep
    100                 : how many times to process
    9.0                 : reference wind speed (see TurbSim manual to understand this)
    80.0                : height of reference wind speed (see TurbSim manual to understand this)
    owez                : outfile prefix
    sampleArray1_U.xy   : name of xy files in each time directory
    15201               : first time (==name of first directory where to find *.xy file) 
    
    Parameters
    ----------
    t0 : float,
        starting time matching the name of the directory in which the *.xy file will be found
    dt : float,
        time increment in between files/directories
    nt : int,
        number of times to process
    zref : float,
        reference height that needs to be saved to *.bts file
    xyname : str,
        exact name of *.xy file which will be found within each time directory
    btsPrefix : str,
        used to name the *.bts file
    theta : float,
        angle by which to rotate the velocity components
    workingPath : str,
        absolute path to where the time directories are, which is where the *.inp file should be written to.

    Returns
    -------
    Doesn't return anything. Procedure that writes the *.inp file in the $pwd.  

    @author: pdoubraw      
    """

    workingPath     = os.path.abspath(os.path.curdir) if workingPath is None else workingPath
    sampleVTKPath   = os.path.join(workingPath,str(t0)+'*','array*U*.vtk')
    sampleVTK       = glob.glob(sampleVTKPath)
    if len(sampleVTK)==0:
        print "Could not find a sampleVTK at {0}".format(sampleVTKPath)
        sys.exit()        
    else:
        sampleVTK   = sampleVTK[0]

    data, meta          = _read_vtkStructured_oneFile(sampleVTK, verbose=False)
    [X,Y,Z,U,V,W]       = data

    iz      = np.argmin(np.abs(Z[0,0,:]-zref))
    uref    = np.mean((np.sqrt(U**2+V**2))[:,:,iz])

    inpName = os.path.join(workingPath,'tsConv.inp')
    f = open(inpName,'w+')
    f.write("{0:.1f}\n".format(theta))
    for x in [meta['ny'],meta['nz']]:
        f.write("{0:d}\n".format(int(x)))
    for x in [meta['dz'],meta['dy'],dt]:
        f.write("{0:.4f}\n".format(x))
    f.write("{0:d}\n".format(int(nt)))
    for x in [uref,zref]:
        f.write("{0:.2f}\n".format(x))
    for x in [btsPrefix,xyname,t0]:
        f.write("{0}\n".format(x))        
    f.close()
    return
#==============================================================================
# 
#==============================================================================
def read_inp(inpPath=None):
    """
    Read the *.inp file that is necessary to convert *.xy to *.bts using of2fast (Sang's fortran code).
    
    "tsConv.inp" entries
    
    0                   : theta
    107                 : numgrid_y (see TurbSim manual to understand this)
    107                 : numgrid_z (see TurbSim manual to understand this)
    1.06                : gridres_z (see TurbSim manual to understand this) 
    1.5                 : gridres_y (see TurbSim manual to understand this)
    1                   : timestep
    100                 : how many times to process
    9.0                 : reference wind speed (see TurbSim manual to understand this)
    80.0                : height of reference wind speed (see TurbSim manual to understand this)
    owez                : outfile prefix
    sampleArray1_U.xy   : name of xy files in each time directory
    15201               : first time (==name of first directory where to find *.xy file)   
    
    Parameters
    ----------
    inpPath : str,
        absolute path to the *.inp file

    Returns
    -------
    out : dict,
        keys are each line of the *.inp file
        
    @author: pdoubraw         
    """    
    inpPath = os.path.join(os.path.abspath(os.path.curdir),"tsConv.inp") if inpPath is None else inpPath   
    fPath = glob.glob(inpPath)
    if len(fPath)==0:
        print "Could not find inp file at {0}".format(inpPath)
        sys.exit()

    out = {}    
    f = open(fPath[0])
    lines = f.readlines()
    out['theta'] = float(lines[0])
    out['ny'], out['nz'] = int(lines[1]), int(lines[2])
    out['dz'], out['dy'], out['dt'] = float(lines[3]), float(lines[4]), float(lines[5])
    out['nt'] = int(lines[6])
    out['uhub'] = float(lines[7])
    out['zhub'] = float(lines[8])
    out['prefix'] = lines[9]
    out['xyname'] = lines[10]
    out['t0'] = float(lines[11])
    out['times'] = np.arange(out['t0'],out['t0']+out['nt']*out['dt'],out['dt'])
    f.close()
    return out
#==============================================================================
# 
#==============================================================================
def vtk2xy(vtkPath, verbose=True):
    """
    Based on a single *.vtk (VTK Structured) file, write out a single *.xy file in the same time directory.
    
    Parameters
    ----------
    vtkPath : str,
        full path to VTK file    
        
    @author: pdoubraw         
    """
    data, meta      = _read_vtkStructured_oneFile(vtkPath,verbose=False)
    [X,Y,Z,U,V,W]   = data

    xq          = np.mean(X[:,0,0]) 
    xi          = np.argmin(np.abs((X[:,0,0] - xq)))   
    nx, ny, nz  = X.shape
    n2d         = ny*nz

    x2d = X[xi,:,:]
    y2d = Y[xi,:,:] 
    z2d = Z[xi,:,:] 
    u2d = U[xi,:,:] 
    v2d = V[xi,:,:] 
    w2d = W[xi,:,:]     

    oPath       = os.path.split(vtkPath)[-1].split('vtk')[0]+'xy'
    oPath       = os.path.join(os.path.split(vtkPath)[0],oPath)
    dataOut     = np.zeros((n2d,6),float)

    dataOut[:,0] = x2d.T.flatten() ; dataOut[:,1] = y2d.T.flatten() ; dataOut[:,2] = z2d.T.flatten()
    dataOut[:,3] = u2d.T.flatten() ; dataOut[:,4] = v2d.T.flatten() ; dataOut[:,5] = w2d.T.flatten()
    np.savetxt(oPath,dataOut,fmt='%26.20E')
    
    if verbose:
        print "Converted {0} to {1}".format(vtkPath,oPath)
    return
#==============================================================================
#         
#==============================================================================
def vtk2bts(workingPath,t0,dt,nt,verbose=False,exePath=None,btsPrefix="prefix"):        
    """
    Given several *.vtk files (within individual time directories), generate a turbsim *.bts file.
    
    Parameters
    ----------    
    workingPath : str,
        absolute path to directory where each time subdirectory is (and within each, a vtk structured file)
    t0 : float,
        starting time matching the name of the directory in which the *.xy file will be found
    dt : float,
        time increment in between files/directories
    nt : int,
        number of times to process
    verbose : bool,
        whether to print out messages        
    exePath : str,
        absolute path to compiled executable "of2fast"
        
    @author: pdoubraw        
    """

    # a tsConv.inp file is necessary
    write_inp(t0,dt,nt,workingPath=workingPath,btsPrefix=btsPrefix)
    
    # loop over time, convert vtk=>xy at each time
    for time in np.arange(t0,t0+nt*dt,dt):

        time = int(time) if np.mod(time,int(time))==0 else time  
              
        # the fortran code was hard coded for three decimals
        timePath = os.path.abspath(os.path.join(workingPath,"{0:.3f}".format(time)))

        # if it doesn't have three decimals, then rename it (else leave it alone)
        if not(os.path.isdir(timePath)):   
            timePathNow = glob.glob(os.path.join(workingPath,"{0}".format(time)))
            if len(timePathNow)==0:
                print "Could not find this vtk file {0}".format(timePathNow)
                sys.exit()        
            else:                
                timePathNow = os.path.abspath(timePathNow[0])
            os.rename(timePathNow,timePath)
            if verbose:
                print "Renamed {0} ==> {1}".format(timePathNow,timePath)
 
        vtkPath = glob.glob(os.path.join(timePath,'array*U*.vtk'))
        
        if len(vtkPath)==0:
            print "Could not find vtk file at {0}".format(vtkPath)
            sys.exit()        
        else:
            vtkPath   = vtkPath[0]
        
        # for this time, get a *.xy file (will overwrite if already there)
        vtk2xy(vtkPath, verbose=False)
    
    # once all times are done, get ready for *.xy => *.bts
    oldPath = os.path.abspath(os.curdir)
    os.chdir(workingPath)
    os.system("ln -s {0} .".format(exePath))

    # if the executable is there, use it (it will read the tsConv.inp and do its thing)
    os.system("./of2fast")
    
    os.chdir(oldPath)
    return