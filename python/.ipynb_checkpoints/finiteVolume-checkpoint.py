import numpy as np




# Trilinearly interpolate to a point in space given a uniformly spaced, Cartesian grid of 3D data.
def interpolateTrilinearUniformCartesian(x,y,z,f,xP,yP,zP):
    
    #i = np.argmin(x <= xP) - 1
    #j = np.argmin(y <= yP) - 1
    #k = np.argmin(z <= zP) - 1
    
    i = np.argmax(x-xP >= 0)
    j = np.argmax(y-yP >= 0)
    k = np.argmax(z-zP >= 0)
    
    if (i >= len(x) - 1):
        i = len(x) - 2
    if (j >= len(y) - 1):
        j = len(y) - 2
    if (k >= len(z) - 1):
        k = len(z) - 2

    xCC = 0.5*(x[i] + x[i+1])
    yCC = 0.5*(y[j] + y[j+1])
    zCC = 0.5*(z[k] + z[k+1])

    fCC = 0.125*(f[i  ,j  ,k  ] + f[i+1,j  ,k  ] + 
                 f[i  ,j+1,k  ] + f[i+1,j+1,k  ] + 
                 f[i  ,j  ,k+1] + f[i+1,j  ,k+1] +
                 f[i  ,j+1,k+1] + f[i+1,j+1,k+1])

    fW = 0.25*(f[i  ,j  ,k  ] + f[i  ,j+1,k  ] +
               f[i  ,j  ,k+1] + f[i  ,j+1,k+1])
    fE = 0.25*(f[i+1,j  ,k  ] + f[i+1,j+1,k  ] +
               f[i+1,j  ,k+1] + f[i+1,j+1,k+1])
    fS = 0.25*(f[i  ,j  ,k  ] + f[i+1,j  ,k  ] +
               f[i  ,j  ,k+1] + f[i+1,j  ,k+1])
    fN = 0.25*(f[i  ,j+1,k  ] + f[i+1,j+1,k  ] +
               f[i  ,j+1,k+1] + f[i+1,j+1,k+1])
    fL = 0.25*(f[i  ,j  ,k  ] + f[i+1,j  ,k  ] +
               f[i  ,j+1,k  ] + f[i+1,j+1,k  ])
    fU = 0.25*(f[i  ,j  ,k+1] + f[i+1,j  ,k+1] +
               f[i  ,j+1,k+1] + f[i+1,j+1,k+1])

    spacing = (x[1]-x[0],y[1]-y[0],z[1]-z[0])
    
    dfdx = (fE - fW)/spacing[0]
    dfdy = (fN - fS)/spacing[1]
    dfdz = (fU - fL)/spacing[2]

    dx = xP - xCC
    dy = yP - yCC
    dz = zP - zCC

    fP = fCC + dfdx*dx + dfdy*dy + dfdz*dz

    #print('i,j,k = ',i,j,k)
    #print(i,x[i],x[i+1])
    #print(j,y[j],y[j+1])
    #print(k,z[k],z[k+1])
    #print(f.shape)

    #print(f[i,j,k],    f[i+1,j,k],
    #      f[i,j+1,k],  f[i+1,j+1,k],
    #      f[i,j,k+1],  f[i+1,j,k+1],
    #      f[i,j+1,k+1],f[i+1,j+1,k+1])

    #print(fW,fE,fS,fN,fL,fU)

    #print(fP)
    
    return fP
    


# Interpolate a quantity from cell centers to cell faces using 2nd-order (linear) interpolation. 
def interpolate2(zMin,z0,dz,z,f,f0,log=False):
    
    nz = len(z)
    zf = np.zeros((nz+1,))
    ff = np.zeros((nz+1,))
    
    
    # Given the cell center coordinates (z) and the grid spacing (dz) and the bottom of the mesh (zMin), 
    # we calculate the cell face coordinates (zf).  We then force the lowest cell face away from zMin to
    # the aerodynamic roughness heigh (z0) if z0 is different than zMin because we usually know the condition 
    # of f at z0.
    if (np.isscalar(dz)):
        dz = dz*np.ones((nz,))
    zf[0] = zMin
    for i in range(nz):
        zf[i+1] = zf[i] + dz[i]
    zf[0] = z0
    
    #print(zf)
    
    
    # If the interpolation is to be done with respect to log(z) instead of z.
    if (log):
        zf = np.log(zf)
        z = np.log(z)
        
    
    # Interpolate f to the cell faces using 2nd-order (linear) interpolation.  Extrapolate to the upper-most face.
    # Set the lowest face to the value at z0 (f0).
    for i in range(nz+1):
        if (i == 0):
            ff[i] = f0
        elif (i == nz):
            ff[i] = 1.5*f[-1] - 0.5*f[-2]
        else:
            dz_ip1 = z[i] - zf[i]
            dz_im1 = zf[i] - z[i-1]
            denom = dz_ip1 + dz_im1
            
            ff[i] = (dz_im1*f[i] + dz_ip1*f[i-1])/denom
            
    if (log):
        zf = np.exp(zf)
    
    
    # Return the cell faces and interpolated quantity.
    return zf,ff
    
    
    
    
    
# Calculate the 2nd-order central finite difference.  We take a finite-volume approach and interpolate to 
# cell faces, then difference between faces to get the finite-difference at cell center.
def central2(zMin,z0,dz,z,f,f0,log=False):
    
    nz = len(z)
    dfdz = np.zeros((nz,))
    
    
    # Interpolate from the cell centers to the cell faces using 2nd-order (linear) interpolation.
    zf, ff = interpolate2(zMin,z0,dz,z,f,f0,log)
    
    
    # If the finite difference is to be done with respect to log(z) instead of z.
    if (log):
        zf = np.log(zf)
        z = np.log(z)

    
    # Calculate the 2nd-order central finite difference.
    for i in range(nz):
        dz_ip1 = zf[i+1] - z[i]
        dz_im1 = z[i] - zf[i]
        df_ip1 = ff[i+1] - f[i]
        df_im1 = f[i] - ff[i]
        denom = 2.0*dz_ip1*dz_im1
        # if this is log(z) finite differencing, we are doing du/dz' * dz'/dz, where dz'/dz = d log(z)/dz = 1/z.
        # so, we must multiply the result by 1/z.  Remember that the z variable is currently log(z), hence the
        # 1/exp(z)
        if (log):
            c = (1.0/np.exp(z[i]))
        else:
            c = 1.0
        dfdz[i] = c*(dz_im1*df_ip1 + dz_ip1*df_im1)/denom
     
    
    # Return the finite-difference.
    return dfdz
        
        
        