# MoninObukhov.py
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 November 2022
#
# This is a class to handle Monin-Obukhov similarity theory.


import numpy as np
import scipy.special as sp




class MoninObukhov:
    
    # Initialize the class.
    def __init__(self):
        self.uStar = []
        self.TFlux = []
        self.sFlux = []
        
        self.T_z0 = []
        self.s_z0 = []
        
        self.L = []
        
        self.z0 = []
        
        self.g = []
        
        self.kappa = []
        
        self.beta_m = []
        self.beta_h = []
        self.beta_s = []
        
        self.gamma_m = []
        self.gamma_h = []
        self.gamma_s = []
        
        self.alpha_h = []
        self.alpha_s = []
        
        
    def showState(self):
        print('uStar = ',self.uStar)
        print('TFlux = ',self.TFlux)
        print('sFlux = ',self.sFlux)
        
        print('T_z0 = ',self.T_z0)
        print('s_z0 = ',self.s_z0)
        
        print('L = ',self.L)
        
        print('z0 = ',self.z0)
        
        print('g = ',self.g)
        
        print('kappa = ',self.kappa)
        print('beta_m = ',self.beta_m)
        print('beta_h = ',self.beta_h)
        print('beta_s = ',self.beta_s)
        print('gamma_m = ',self.gamma_m)
        print('gamma_h = ',self.gamma_h)
        print('gamma_s = ',self.gamma_s)
        print('alpha_h = ',self.alpha_h)
        print('alpha_s = ',self.alpha_s)
        
        
        
    # Compute the velocity profile.
    def compute_velocity_profile(self,uStar,z,L):
        u = np.zeros((len(z),))
        psi_m_z0 = self.compute_psi_m(self.z0,L)
        for i in range(len(z)):
            psi_m_z = self.compute_psi_m(z[i],L)
            u[i] = (uStar/self.kappa)*(np.log(z[i]/self.z0) - (psi_m_z - psi_m_z0))
            
        return u
        
        
        
    # Compute the temperature profile.
    def compute_temperature_profile(self,uStar,T_z0,TFlux,z,L):
        T = np.zeros((len(z),))
        psi_h_z0 = self.compute_psi_h(self.z0,L)
        for i in range(len(z)):
            psi_h_z = self.compute_psi_h(z[i],L)
            T[i] = T_z0 + ((-self.alpha_h*TFlux)/(self.kappa*uStar))*(np.log(z[i]/self.z0) - (psi_h_z - psi_h_z0))
            
        return T
        
        
        
    # Compute the scalar profile.
    def compute_scalar_profile(self,uStar,s_z0,sFlux,z,L):
        s = np.zeros((len(z),))
        psi_s_z0 = self.compute_psi_s(self.z0,L)
        for i in range(len(z)):
            psi_s_z = self.compute_psi_s(z[i],L)
            s[i] = s_z0 + ((-self.alpha_s*sFlux)/(self.kappa*uStar))*(np.log(z[i]/self.z0) - (psi_s_z - psi_s_z0))
            
        return s
        
        
        
    # Compute friction velocity and surface temperature given surface temperature flux.
    def compute_surface_conditions_given_Tflux(self,TFlux,z1,U_z1,T_z1,tol,iterMax,verbose=False):
        uStar = max(1.0E-6,(self.kappa*U_z1)/(np.log(z1/self.z0)))
        uStar1 = uStar
        uStar0 = uStar + 2.0*tol
        phi_m = 1.0
        psi_m = 0.0
        T_z0 = T_z1
        
        iter = 0
        while ((abs(uStar1 - uStar0) > tol) and (iter < iterMax)):
            uStar0 = uStar1
            L = self.compute_L(uStar,TFlux,T_z0)
            uStar = self.compute_uStar(z1,U_z1,L)
            T_z0 = self.compute_T_z0(z1,T_z1,L,TFlux,uStar)
            uStar1 = uStar
            iter = iter + 1
            
            if (verbose):
                print('-iter: ',iter)
                print(' L = ',L)
                print(' friction velocity = ',uStar)
                print(' T_surface =', T_z0)
                print(' ')
                
        return uStar,T_z0,L
        
        
        
    # Compute friction velocity and surface temperature flux given surface temperature.
    def compute_surface_conditions_given_Tsurface(self,T_z0,z1,U_z1,T_z1,tol,iterMax,verbose=False):
        uStar = max(1.0E-6,(self.kappa*U_z1)/(np.log(z1/self.z0)))
        uStar1 = uStar
        uStar0 = uStar + 2.0*tol
        phi_m = 1.0
        psi_m = 0.0
        TFlux = 0.1
        
        iter = 0
        while ((abs(uStar1 - uStar0) > tol) and (iter < iterMax)):
            uStar0 = uStar1
            L = self.compute_L(uStar,TFlux,T_z0)
            uStar = self.compute_uStar(z1,U_z1,L)
            TFlux = self.compute_TFlux(z1,T_z1,L,T_z0,uStar)
            uStar1 = uStar
            iter = iter + 1
            
            if (verbose):
                print('-iter: ',iter)
                print(' L = ',L)
                print(' friction velocity = ',uStar)
                print(' T_surface_flux =', TFlux)
                print(' ')
                
        return uStar,TFlux,L

        
        
        
    # Compute phi_m.
    def compute_phi_m(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        if ((L < 1.0E6) and (L > 0)):
            phi_m = 1.0 + self.gamma_m*(z/L) 
        elif ((L > -1.0E6) and (L <= 0)):
            phi_m = np.power((1.0 - self.beta_m*(z/L)),-0.25)
        else:
            phi_m = np.ones((nz,))
        return phi_m
    
    
    
    # Compute phi_h.
    def compute_phi_h(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        if ((L < 1.0E6) and (L > 0)):
            phi_h = self.alpha_h + self.gamma_h*(z/L) 
        elif ((L > -1.0E6) and (L <= 0)):
            phi_h = self.alpha_h*np.power((1.0 - self.beta_h*(z/L)),-0.50)
        else:
            phi_h = np.ones((nz,))
        return phi_h
    
    
    
    # Compute phi_s.
    def compute_phi_s(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        if ((L < 1.0E6) and (L > 0)):
            phi_s = self.alpha_s + self.gamma_s*(z/L) 
        elif ((L > -1.0E6) and (L <= 0)):
            phi_s = self.alpha_s*np.power((1.0 - self.beta_s*(z/L)),-0.50)
        else:
            phi_s = np.ones((nz,))
        return phi_s
    
    
    
    # Compute psi_m.
    def compute_psi_m(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        phi_m_ = self.compute_phi_m(z,L)
        if ((L < 1.0E6) and (L > 0)):
            psi_m = 1.0 -  phi_m_
        elif ((L > -1.0E6) and (L <= 0)):
            psi_m = 2.0*np.log(0.5*(1.0 + np.power(phi_m_,-1.0))) \
                  +     np.log(0.5*(1.0 + np.power(phi_m_,-2.0))) \
                  - 2.0*np.arctan(np.power(phi_m_,-1.0)) \
                  + 0.5*np.pi
        else:
            psi_m = np.zeros((nz,))
        return psi_m
    
    
    
    # Compute psi_h.
    def compute_psi_h(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        phi_h_ = self.compute_phi_h(z,L)
        if ((L < 1.0E6) and (L > 0)):
            psi_h = 1.0 - (phi_h_/self.alpha_h)
        elif ((L > -1.0E6) and (L <= 0)):
            psi_h = 2.0*np.log(0.5*(1.0 + np.power((phi_h_/self.alpha_h),-1.0)))
        else:
            psi_h = np.zeros((nz,))
        return psi_h
    
    
    
    # Compute psi_s.
    def compute_psi_s(self,z,L):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
        phi_s_ = self.compute_phi_s(z,L)
        if ((L < 1.0E6) and (L > 0)):
            psi_s = 1.0 - (phi_s_/self.alpha_s)
        elif ((L > -1.0E6) and (L <= 0)):
            psi_s = 2.0*np.log(0.5*(1.0 + np.power((phi_s_/self.alpha_s),-1.0)))
        else:
            psi_s = np.zeros((nz,))
        return psi_s
    
    
    
    # Compute Obukhov length.
    def compute_L(self,uStar,TFlux,T_z0):
        eps = 1.0E-6
        if (abs(TFlux) > eps):
            L = -np.power(uStar,3.0)/(self.kappa*self.g*TFlux/T_z0)
        else:
            L = 1.0E10  
        self.L = L
        return L
    
    
    
    # Compute friction velocity given velocity at a height z and Obukhov length.
    def compute_uStar(self,z,U,L):
        psi_m_z = self.compute_psi_m(z,L)
        psi_m_z0 = self.compute_psi_m(self.z0,L)
        uStar = (self.kappa*U) / \
                (np.log(z/self.z0) - (psi_m_z - psi_m_z0))
        uStar = max(uStar,1.0E-6)
        self.uStar = uStar
        return uStar
    
    
    
    # Compute theta at z0 given theta at a height z and Obukhov length.
    def compute_T_z0(self,z,T,L,TFlux,uStar):
        psi_h_z = self.compute_psi_h(z,L)
        psi_h_z0 = self.compute_psi_h(self.z0,L)
        T_z0 = T + ((self.alpha_h*TFlux)/(self.kappa*uStar))*(np.log(z/self.z0) - (psi_h_z - psi_h_z0))
        return T_z0
    
    
    
    # Compute temperature flux given theta at a height z and Obukhov length.
    def compute_TFlux(self,z,T,L,T_z0,uStar):
        psi_h_z = self.compute_psi_h(z,L)
        psi_h_z0 = self.compute_psi_h(self.z0,L)
        TFlux = (-(T - T_z0)*self.kappa*uStar) / \
                (self.alpha_h*(np.log(z/self.z0) - (psi_h_z - psi_h_z0)))
        return TFlux
    
    
    
    # Compute scalar at z0 given scalar at a height z and Obukhov length.
    def compute_s_z0(self,z,s,L,sFlux,uStar):
        psi_s_z = self.compute_psi_s(z,L)
        psi_s_z0 = self.compute_psi_s(self.z0,L)
        s_z0 = s + ((self.alpha_s*sFlux)/(self.kappa*uStar))*(np.log(z/self.z0) - (psi_s_z - psi_s_z0))
        return s_z0
    
    
    
    # Compute temperature flux given theta at a height z and Obukhov length.
    def compute_sFlux(self,z,s,L,s_z0,uStar):
        psi_s_z = self.compute_psi_s(z,L)
        psi_s_z0 = self.compute_psi_s(self.z0,L)
        sFlux = (-(s - s_z0)*self.kappa*uStar) / \
                (self.alpha_s*(np.log(z/self.z0) - (psi_s_z - psi_s_z0)))
        return sFlux
    
    
    # Add wiggles to the non-dimensional velocity gradient profile.
    def add_wiggles(self,z,wigglePeak,zA):
        nz = len(z)
        
        #wigglePeak*z[0]^b/z^b = 0.01*wigglePeak
        
        #z^b = z[0]^b/0.01
        
        #z^b / z[0]^b = 1/0.01
        
        #(z/z[0])^b = 1/0.01
        #b*ln(z/z[0]) = ln(1/0.01)
        #b = ln(1/0.01)/ln(z/z[0])
        
        # background attenuation function
        b = 1.5
        b = np.log(1.0/0.01)/np.log(zA/z[0])
        a = wigglePeak*np.power(z[0],b)
        f = np.zeros((nz,))
        g = np.zeros((nz,))

        for i in range(nz):
            f[i] = a/np.power(z[i],b)
            if (i % 2):
                g[i] = -1.0
            else:
                g[i] = 1.0
            
        fg = f*g
        
        return fg

        
    
    # An analytical velocity gradient "overshoot" function.  It is the product of a Gaussian and a linear function.
    def compute_phi_m_overshoot(self,z,zPeak,overshootStrength):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
            
        # function 1, the Gaussian
        a = np.sqrt(2.0*np.square(zPeak))
        b = (overshootStrength-1.0)*(1.0 / (zPeak * np.exp(-0.5)))
        f1 = b * np.exp(-np.square(z/a))
        f1P = -f1*2.0*(z/(np.square(a)))

        # function 2, linear
        f2 = z
        f2P = np.ones((nz,))

        # the complete function, f = f1*f2
        f = 1.0 + f1*f2
        fP = f1P*f2 + f1*f2P
        
        return f
    
    
    
    # An analytical velocity profile with the analytical velocity gradient "overshoot" function built in.
    # Note, I have only found an analytical function for the neutral situation.
    def compute_velocity_profile_overshoot(self,z,uStar,zPeak,overshootStrength):
        nz = 1
        if not (np.isscalar(z)):
            nz = len(z)
           
        a = np.sqrt(2.0*np.square(zPeak))
        b = (overshootStrength-1.0)*(1.0 / (zPeak * np.exp(-0.5)))
        
        u = (uStar/self.kappa)*(0.5*np.sqrt(np.pi)*a*b*(sp.erf(z/a) - sp.erf(self.z0/a)) + np.log(z/self.z0))
        
        return u
    
    
    
    # Integrate a profile given phi_m.
    def integrate_velocity_profile(self,L,uStar,zPeak,overshootStrength,wigglePeak,zA,z,zf,verbose=False):
        tol = 1.0E-3
        diff = 1.0E5
        
        z0 = zf[0]
        zMax = zf[-1]
        dz = zf[-1] - zf[-2]
        nz = len(z)
        
        zz = z
        w = self.add_wiggles(zz,wigglePeak,zA)
        
        # Integrate the velocity profile using the given resolution
        os = self.compute_phi_m_overshoot(z,zPeak,overshootStrength)
        phi_m = w + os*self.compute_phi_m(z,L)
        dudz = (uStar/(z*self.kappa))*phi_m
        dudzP = (uStar/self.kappa)*phi_m
        
        # The actual integration is done with respect to z' = log(z)
        uf = np.zeros(len(zf))
        zfP = np.log(zf)
        for i in range(1,len(zf)):
            dzP = zfP[i] - zfP[i-1]
            uf[i] = uf[i-1] + dzP*dudzP[i-1]
            
        ufTopOld = uf[-1]
        
        # The integration is repeated over and over using consecutively doubled
        # grid resolution until the value at the top of the profile converges
        # to within the tolerance
        dz_ = dz
        nz_ = nz
        m = 1
        while (diff > tol):
            dz_ = dz_/2.0
            nz_ = 2*nz_
            zf_ = np.linspace(0.0,zMax,nz_+1)
            zf_[0] = z0
            z_ = np.linspace(0.5*dz_,zMax-0.5*dz_,nz_)
            ww = np.interp(z_,zz,w)
            
            os = self.compute_phi_m_overshoot(z_,zPeak,overshootStrength)
            phi_m = ww + os*self.compute_phi_m(z_,L)
            dudz = (uStar/(z_*self.kappa))*phi_m
            dudzP = (uStar/self.kappa)*phi_m
     
            del(uf)
            uf = np.zeros(len(zf_))
            zfP = np.log(zf_)
            for i in range(1,len(zf_)):
                dzP = zfP[i] - zfP[i-1]
                uf[i] = uf[i-1] + dzP*dudzP[i-1]
               
            diff = np.abs(uf[-1] - ufTopOld)
            ufTopOld = uf[-1]
            
            if (verbose):
                print('iter: ',m,'diff: ',diff)
                
            m = m+1
            
        
        # Interpolate from the refined profile back to the given resolution.
        u = np.interp(z,zf_,uf)
        return u
        