# planarAverageTools.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 31 July 2023
#
# This is a class that contains methods to process planar-averaged
# data that comes out of atmospheric boundary layer simulations



import numpy as np
import finiteVolume as fv



class planarAverageTools:

    
    # initialize the class.
    def __init__(self,t_,zCell_,zFace_,dz_,nzCell_,nzFace_,data_,profileData_):
        self.t = t_
        
        self.zCell = zCell_
        self.zFace = zFace_
        self.dz = dz_
        self.nzCell = nzCell_
        self.nzFace = nzFace_
        
        self.data = data_
        self.profileData = profileData_

        self.nAverageWindows = 0
        self.timeWindowStart = []
        self.timeWindowEnd = []
        self.indexWindowStart = []
        self.indexWindowEnd = []
        
        self.q = []
        self.theta_surf = []
        self.uStar = []
        self.wStar = []
        self.L = []
        self.zi = []
        self.zig = []
        self.zjet = []
        self.dpdx = []
        self.dpdy = []

        self.u_cell = []
        self.v_cell = []
        self.w_cell = []
        self.hvelmag_cell =[]
        self.windSpeed_cell = []
        self.windDir_cell = []
        self.theta_cell = []
        self.mueff_cell = []
        self.thetatheta_r_cell = []
        self.utheta_r_cell = []
        self.vtheta_r_cell = []
        self.wtheta_r_cell = []
        self.uu_r_cell = []
        self.vv_r_cell = []
        self.ww_r_cell = []
        self.uv_r_cell = []
        self.uw_r_cell = []
        self.vw_r_cell = []
        self.k_r_cell = []
        self.shearStressMag_r_cell = []
        self.uuu_r_cell = []
        self.vvv_r_cell = []
        self.www_r_cell = []
        self.utheta_sfs_cell = []
        self.vtheta_sfs_cell = []
        self.wtheta_sfs_cell = []
        self.uv_sfs_cell = []
        self.uw_sfs_cell = []
        self.vw_sfs_cell = []
        self.shearStressMag_sfs_cell = []
        
        self.S_cell = []
        self.N_cell = []
        self.Ri_cell = []

        self.phi_m_cell = []
        self.phi_m_face = []
        self.phi_m_cell_log = []
        self.scriptR_cell = []
        self.scriptR_face = []
        self.nuLES = []
        self.ReLES = []
        
        self.eddyTurnOverTime = []
        
    
    
    
    # Set up the time intervals over which averaging happens.
    def setUpTimeWindowMean(self,timeAverageStart,timeAverageEnd,timeAverageWindowWidth):
        
        ### Set up the time averaging windows.
        self.nAverageWindows = 0
        self.timeWindowStart = []
        self.timeWindowEnd = []
        self.indexWindowStart = []
        self.indexWindowEnd = []

        # first, set up all the averaging windows as user specified, which may extend beyond the actual data,
        # but this gets cut down later.
        self.nAverageWindows = int(np.ceil((timeAverageEnd-timeAverageStart)/timeAverageWindowWidth))

        for i in range(self.nAverageWindows):
            self.timeWindowStart.append(timeAverageStart + i*timeAverageWindowWidth)
            self.timeWindowEnd.append(timeAverageStart + (i+1)*timeAverageWindowWidth)

        # next, cut down the list of averaging intervals to fit the actual time series we have.
        iStart = np.argmax((self.t[0] - self.timeWindowStart) >= 0)
        iEnd = np.argmax((self.t[-1] - self.timeWindowEnd) <= 0)
        #print(self.t[0],(self.t[0] - self.timeWindowStart) >= 0,np.argmax((self.t[0] - self.timeWindowStart) >= 0))
        #print(self.t[-1],(self.t[-1] - self.timeWindowEnd) <= 0,np.argmax((self.t[-1] - self.timeWindowEnd) <= 0))
        #print(iStart,iEnd)
        if (iEnd == iStart):
            iEnd = iEnd + 1
        #print(iStart,iEnd)
        self.timeWindowStart = self.timeWindowStart[iStart:iEnd]
        self.timeWindowEnd = self.timeWindowEnd[iStart:iEnd]
        self.nAverageWindows = max((iEnd - iStart),1)

        # find the time series indices for the start and end of each averaging window
        for i in range(self.nAverageWindows):
            self.indexWindowStart.append(np.argmax((self.t - self.timeWindowStart[i]) >= 0))
            self.indexWindowEnd.append(np.argmax((self.t - self.timeWindowEnd[i]) >= 0))

        print('Created', self.nAverageWindows, 'time series entries')

    
    
    
    # Assemble a series of time-averages of a given window width given the full
    # time series of planar averages.
    def timeWindowMeanAMRWind(self):

        # Reset the list of averaged variables.
        self.q = []
        self.theta_surf = []
        self.uStar = []
        self.wStar = []
        self.L = []
        self.zi = []
        self.dpdx = []
        self.dpdy = []

        self.u_cell = []
        self.v_cell = []
        self.w_cell = []
        self.hvelmag_cell =[]
        self.windSpeed_cell = []
        self.windDir_cell = []
        self.theta_cell = []
        self.mueff_cell = []
        self.thetatheta_r_cell = []
        self.utheta_r_cell = []
        self.vtheta_r_cell = []
        self.wtheta_r_cell = []
        self.uu_r_cell = []
        self.vv_r_cell = []
        self.ww_r_cell = []
        self.uv_r_cell = []
        self.uw_r_cell = []
        self.vw_r_cell = []
        self.k_r_cell = []
        self.k_sfs_cell = []
        self.shearStressMag_r_cell = []
        self.uuu_r_cell = []
        self.vvv_r_cell = []
        self.www_r_cell = []
        self.utheta_sfs_cell = []
        self.vtheta_sfs_cell = []
        self.wtheta_sfs_cell = []
        self.uu_sfs_cell = []
        self.vv_sfs_cell = []
        self.ww_sfs_cell = []
        self.uv_sfs_cell = []
        self.uw_sfs_cell = []
        self.vw_sfs_cell = []
        self.shearStressMag_sfs_cell = []


        for i in range(self.nAverageWindows):

            jS = self.indexWindowStart[i]
            jE = self.indexWindowEnd[i]

            self.q.append(np.mean(self.data[jS:jE,0],axis=0))
            self.theta_surf.append(np.mean(self.data[jS:jE,1],axis=0))
            self.uStar.append(np.mean(self.data[jS:jE,2],axis=0))
            self.wStar.append(np.mean(self.data[jS:jE,3],axis=0))
            self.L.append(np.mean(self.data[jS:jE,4],axis=0))
            self.zi.append(np.mean(self.data[jS:jE,5],axis=0))
            self.dpdx.append(np.mean(self.data[jS:jE,6],axis=0))
            self.dpdy.append(np.mean(self.data[jS:jE,7],axis=0))


            # Do the time averaging of profiles.
            self.u_cell.append(np.mean(self.profileData[jS:jE,:,0],axis=0))
            self.v_cell.append(np.mean(self.profileData[jS:jE,:,1],axis=0))
            self.w_cell.append(np.mean(self.profileData[jS:jE,:,2],axis=0))
            self.hvelmag_cell.append(np.mean(self.profileData[jS:jE,:,3],axis=0))
            self.windSpeed_cell.append(np.sqrt(np.square(self.u_cell[i]) + np.square(self.v_cell[i])))
            self.windDir_cell.append(np.arctan2(self.v_cell[i],self.u_cell[i])*(180.0/np.pi))
            self.theta_cell.append(np.mean(self.profileData[jS:jE,:,4],axis=0))
            self.mueff_cell.append(np.mean(self.profileData[jS:jE,:,5],axis=0))
            self.thetatheta_r_cell.append(np.mean(self.profileData[jS:jE:,6],axis=0))
            self.utheta_r_cell.append(np.mean(self.profileData[jS:jE,:,7],axis=0))
            self.vtheta_r_cell.append(np.mean(self.profileData[jS:jE,:,8],axis=0))
            self.wtheta_r_cell.append(np.mean(self.profileData[jS:jE,:,9],axis=0))
            self.uu_r_cell.append(np.mean(self.profileData[jS:jE,:,10],axis=0))
            self.vv_r_cell.append(np.mean(self.profileData[jS:jE,:,11],axis=0))
            self.ww_r_cell.append(np.mean(self.profileData[jS:jE,:,12],axis=0))
            self.uv_r_cell.append(np.mean(self.profileData[jS:jE,:,13],axis=0))
            self.uw_r_cell.append(np.mean(self.profileData[jS:jE,:,14],axis=0))
            self.vw_r_cell.append(np.mean(self.profileData[jS:jE,:,15],axis=0))
            self.k_r_cell.append((0.5*(self.uu_r_cell[i] + self.vv_r_cell[i] + self.ww_r_cell[i])))
            self.shearStressMag_r_cell.append(np.sqrt(np.square(self.vw_r_cell[i]) + np.square(self.uw_r_cell[i])))
            self.uuu_r_cell.append(np.mean(self.profileData[jS:jE,:,16],axis=0))
            self.vvv_r_cell.append(np.mean(self.profileData[jS:jE,:,17],axis=0))
            self.www_r_cell.append(np.mean(self.profileData[jS:jE,:,18],axis=0))
            self.utheta_sfs_cell.append(np.mean(self.profileData[jS:jE,:,19],axis=0))
            self.vtheta_sfs_cell.append(np.mean(self.profileData[jS:jE,:,20],axis=0))
            self.wtheta_sfs_cell.append(np.mean(self.profileData[jS:jE,:,21],axis=0))
            self.uv_sfs_cell.append(np.mean(self.profileData[jS:jE,:,22],axis=0))
            self.uw_sfs_cell.append(np.mean(self.profileData[jS:jE,:,23],axis=0))
            self.vw_sfs_cell.append(np.mean(self.profileData[jS:jE,:,24],axis=0))
            self.shearStressMag_sfs_cell.append(np.sqrt(np.square(self.vw_sfs_cell[i]) + np.square(self.uw_sfs_cell[i])))







    

    # Derive SGS turbulent kinetic energy.  If other than the one-equation model is used, we can
    # derive SGS k using eddy viscosity.
    def derive_k_sfs_using_mut(self,Ck,l,mu,rho):
        self.k_sfs_cell = []

        for i in range(self.nAverageWindows):
            self.k_sfs_cell.append(np.square((1.0/(Ck*l))*((self.mueff_cell[i]-mu)/rho)))
            self.uu_sfs_cell.append((2.0/3.0)*self.k_sfs_cell[i])
            self.vv_sfs_cell.append((2.0/3.0)*self.k_sfs_cell[i])
            self.ww_sfs_cell.append((2.0/3.0)*self.k_sfs_cell[i])




    


    
    
    # Compute Richardson number, shear, Brunt-Vaisalla frequency.
    def computeRi_N_S(self,z0,g,theta0,thetaSurface):
        self.S_cell = []
        self.N_cell = []
        self.Ri_cell = []

        for i in range(self.nAverageWindows):
            dTdz_ = fv.central2(0.0,z0,self.dz,self.zCell,self.theta_cell[i],thetaSurface,log=True)
            dUdz_ = fv.central2(0.0,z0,self.dz,self.zCell,self.u_cell[i],0.0,log=True)
            dVdz_ = fv.central2(0.0,z0,self.dz,self.zCell,self.v_cell[i],0.0,log=True)
            dUdzMag_ = np.sqrt(np.square(dUdz_)+np.square(dVdz_))

            self.N_cell.append(np.sqrt(((g/theta0)*dTdz_)))
            self.S_cell.append(dUdzMag_)
            self.Ri_cell.append(np.square(self.N_cell[i])/np.square(self.S_cell[i]))








    

    # Compute "high-accuracy zone" metrics.  These are metrics given in the Brasseur and Wei paper that describe the
    # the quality of an LES of a wall-bounded flow.
    def computeHAZ(self,z0,kappa):          
            
            
        self.phi_m_cell = []
        self.phi_m_face = []
        self.phi_m_cell_log = []
        self.scriptR_cell = []
        self.scriptR_face = []
        self.nuLES = []
        self.ReLES = []
        
        
        # Calculate dhvelmag/dz on cell faces and cell centers and using log-based derivatives/interpolation.
        # In the end I use the log-based one because it gives the exact answer given a log input.  I've shown
        # it to work more accurately in computing du/dz for u(z) that is log-like than standard central 
        # differencing.            
        for i in range(self.nAverageWindows):

            dhvelmagdz_face_ = np.zeros((self.nzFace,))
            dhvelmagdz_face_[1:-1] = (self.hvelmag_cell[i][1:] - self.hvelmag_cell[i][0:-1])/self.dz
            dhvelmagdz_cell_ = np.gradient(self.hvelmag_cell[i],self.dz)
            dhvelmagdz_cell_log_ = fv.central2(0.0,z0,self.dz,self.zCell,self.hvelmag_cell[i],0.0,True)
            dudz_cell_log_ = fv.central2(0.0,z0,self.dz,self.zCell,self.u_cell[i],0.0,True)
            dvdz_cell_log_ = fv.central2(0.0,z0,self.dz,self.zCell,self.v_cell[i],0.0,True)

            # Compute phi_m.
            self.phi_m_cell.append(    (kappa/self.uStar[i])*dhvelmagdz_cell_    *self.zCell)
            self.phi_m_face.append(    (kappa/self.uStar[i])*dhvelmagdz_face_    *self.zFace)
            self.phi_m_cell_log.append((kappa/self.uStar[i])*dhvelmagdz_cell_log_*self.zCell)

            # Compute script R.
            self.scriptR_cell.append(self.shearStressMag_r_cell[i][0]/self.shearStressMag_sfs_cell[i][0])

            # Compute Re_LES
            self.nuLES.append(self.shearStressMag_sfs_cell[i][0]/np.sqrt(np.square(dudz_cell_log_[0]) + np.square(dvdz_cell_log_[0])))
            self.ReLES.append(self.zi[i]*self.uStar[i]/self.nuLES[i])

   



    # Compute the eddy turn-over time.  This may be useful in determining averaging times.
    def computeEddyTurnOverTime(self):
        
        self.eddyTurnOverTime = []
        
        for i in range(self.nAverageWindows):
            self.eddyTurnOverTime.append(self.zig[i]/self.uStar[i])

   



    # Compute the inversion height using various methods as outlined in Sullivan et al. (1998).
    def compute_zi(self,zMin=0.0):
        
        self.zig = []
        
        for i in range(self.nAverageWindows):
            grad = fv.central2(0.0,0.0,self.dz,self.zCell,self.theta_cell[i],self.theta_cell[i][0],False)
            
            gradMax = -1.0E10
            indMax = 0
            for j in range(len(grad)):
                if (self.zCell[j] >= zMin):
                    if (grad[j] > gradMax):
                        gradMax = grad[j]
                        indMax = j
                    
            self.zig.append(self.zCell[indMax])

   



    # Compute the jet height using .
    def compute_zjet(self,zMin=0.0):
        
        self.zjet = []
        
        for i in range(self.nAverageWindows):
            indMax = 0
            velMax = -1.0E10
            for j in range(len(self.windSpeed_cell[i])):
                if (self.zCell[j] >= zMin):
                    if (self.windSpeed_cell[i][j] > velMax):
                        velMax = self.windSpeed_cell[i][j]
                        indMax = j
                    
            self.zjet.append(self.zCell[indMax])
                
            
      
            
            
            
            
            