# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 08:54:48 2017

@author: pdoubraw
"""
import os, glob, struct
from math import gamma as gammaFunc
import numpy as np
from subprocess import call
from shutil import rmtree
from scipy import fftpack
import matplotlib.pyplot as plt
from pandas import read_csv

def parse(string):
    """
    Convert strings to float.
    """
    try:
        return float(string)
    except:
        return np.nan

power_law = lambda z, zref, uref, alpha: uref*(z/zref)**alpha

class stochasticTurbulence:

    def __init__(self,prefix='prefix',D=154.0):
        """
        Instantiate the object.

        Parameters
        ----------
        prefix : string,
            can be used to read existing files or generate new ones.
        """
        self.prefix = prefix
        self.D = D
        self.R = D/2.0

    def initialize_mannBox(self, seed=1209, flagHF=True, gamma=0,
                        nX=64, nY=32, nZ=32, dX=10.0, dY=5.0, dZ=5.0,
                    target_TI=0.1, target_uHub=8.0, zHub=90, T=3900.):
        """
        In case you don't have a Mann box yet and want to generate one, first
        you must initialize it with the desired parameters.

        Parameters
        ----------
        target_TI : float,
            hub height turbulence intensity [-]
        U : float,
            hub height mean wind speed [m/s]
        z : float,
            hub height [m]
        gamma : float,
            non-dimensional shear distortion parameter, set to zero for isotropic turbulence [-]
        nX : int,
            number of points along streamwise direction [-]
        nY : int,
            number of points along cross-stream direction [-]
        nZ : int,
            number of points along vertical direction [-]
        dX : float,
            resolution along streamwise direction [m]
        dY : float,
            resolution along cross-stream direction [m]
        dZ : float,
            resolution along vertical direction [m]
        """
        check = np.array([ (np.log2(x)%1) for x in [nX,nY,nZ] ])
        if any(check!=0):
            raise Exception("nX, nY, and nZ must be a power of 2")

        self.seed   = seed
        self.flagHF = flagHF
        self.gamma  = gamma
        self.nX = int(nX) ; self.nY = int(nY) ; self.nZ = int(nZ)
        self.dX = dX ; self.dY = dY ; self.dZ = dZ

        if target_TI>1:
            target_TI = target_TI/100.0

        self.target_TI = target_TI ; self.target_uHub = target_uHub
        self.sigma_1_NTM = self.target_TI * (0.75 * self.target_uHub + 5.6)
        self.sigma_1 = self.target_uHub * self.target_TI
        self.sigma_iso = 0.55 * self.sigma_1
        self.lambda_1 = 42.0 if zHub>60 else 0.7*zHub
        self.L = 0.8 * self.lambda_1
        constants = (55./9.) * (gammaFunc(5./6.)) / ((np.sqrt(np.pi))*gammaFunc(1./3.))
        self.alphaEpsilon = self.sigma_iso**2 * constants * self.L**(-2./3.)
        self.length = self.nX * self.dX
        self.width = self.nY * self.dY
        self.height = self.nZ * self.dZ
        self.zBot = 0.0
        self.T = T
        self.dT = self.dX / self.target_uHub
        self.zHub = zHub

    def generate_mannBox(self,pathToExecutable,outPath=None, \
                             justCreateBat=True,overwriteDir=False):
        """
        Generate a Mann Box.

        Parameters
        ----------
        pathToExecutable : string,
            where to find executable of Mann command line tool, which can be downloaded from http://www.hawc2.dk/download/pre-processing-tools
        """
        if outPath is None:
            outPath = pathToExecutable

        if not(hasattr(self,'alphaEpsilon')):
            raise Exception("First you must initialize the box parameters.")

        if overwriteDir:
            if os.path.isdir(outPath):
                print("Deleting and re-creating: {0}".format(outPath))
                rmtree(outPath)
            os.mkdir(outPath)

        callstr = "{0} {1} {2:.5f} {3:.1f} {4:.1f} {5:d} {6:d} {7:d} {8:d} {9:.2f} {10:.2f} {11:.2f} {12}".format(
                "mann_turb_x64.exe", self.prefix, self.alphaEpsilon, self.L,
                       self.gamma, self.seed, self.nX, self.nY, self.nZ,
                       self.dX, self.dY, self.dZ, self.flagHF)

        if justCreateBat:
            outBat = os.path.join(outPath,"run_{0}.bat".format(self.prefix))
            with open(outBat, "w") as text_file:
                text_file.write(callstr)
            print ("Successfully generated {0}".format(outBat))
        else:
            os.chdir(pathToExecutable)
            exePath = glob.glob(os.path.join(pathToExecutable,'mann*.exe'))
            if len(exePath)==0:
                raise Exception("Could not find executable.")
            call(callstr)
            outFilePaths = glob.glob(os.path.join(os.getcwd(),"{0}*.bin".format(self.prefix)))
            outFilePaths.append(glob.glob(os.path.join(os.getcwd(),"{0}*.txt".format(self.prefix)))[0])
            for outFilePath in outFilePaths:
                os.rename(outFilePath,os.path.join(outPath,os.path.split(outFilePath)[-1]))

            print ("Successfully generated turbulence." )
            self.readMannBox(outPath)

    def readMannTxt(self,pathToBinaries):
        """
        Reads in the text file that is generated as an output of the CL tool.

        Parameters
        ----------
        pathToBinaries : string,
            where to find the .bin / .txt files created by the Mann CL Tool
        """
        txtPath = glob.glob(os.path.join(pathToBinaries,
             '{0}_u.txt'.format(self.prefix)))
        if len(txtPath)==0:
            raise Exception("Could not find text file")
        with open(txtPath[0]) as f:
            content = f.readlines()
        content = np.array([x.strip() for x in content])

        self.alphaEpsilon = float(content[np.argmax([ x.startswith('Alfa') for x in content ])].split(' ')[-1])
        self.L = float(content[np.argmax([ x.startswith('L_mann') for x in content ])].split(' ')[-1])
        self.gamma = float(content[np.argmax([ x.startswith('Gamma') for x in content ])].split(' ')[-1])

        uLine = content[np.argmax([ x.startswith('Length') for x in content ])].replace(':',' ').split(' ')
        vLine = content[np.argmax([ x.startswith('Width') for x in content ])].replace(':',' ').split(' ')
        wLine = content[np.argmax([ x.startswith('Height') for x in content ])].replace(':',' ').split(' ')

        self.length, self.dX, self.nX = [ parse(x) for x in uLine if ~np.isnan(parse(x)) ]
        self.width, self.dY, self.nY = [ parse(x) for x in vLine if ~np.isnan(parse(x)) ]
        self.height, self.dZ, self.nZ = [ parse(x) for x in wLine if ~np.isnan(parse(x)) ]

        self.nX = int(self.nX) ; self.nY = int(self.nY) ; self.nZ = int(self.nZ)

    def readMannBox(self,pathToBinaries,zHub=99.3,target_uHub=8.0,target_TI=0.103,T=3900.0):
        """
        Read a previously generated Mann box.

        Parameters
        ----------
        pathToBinaries : string,
            where to find binary outputs
        zHub : float,
            optional argument used to identify which point represents hub location, needed once a grid is generated
        """
        self.readMannTxt(pathToBinaries)
        self.fileDir = pathToBinaries
        components  = ['u','v','w']
        data = {}
        nNumbers = self.nX*self.nY*self.nZ
        for component in components:
            componentPath = glob.glob(os.path.join(pathToBinaries,
             '{0}_{1}.bin'.format(self.prefix,component)))
            if len(componentPath)==0:
                raise Exception("Could not find file for component {0}".format(component))
            print ("Opening file {0}...".format(componentPath[0]))
            with open(componentPath[0], mode='rb') as file:
                nBytes = os.path.getsize(componentPath[0])
                if (nBytes/4.0 != nNumbers):
                    raise Exception("File size does not match expectations")
                fileContent = file.read()

                data[component] = np.array(struct.unpack('f'*int(nNumbers),fileContent)).reshape(
                                    (self.nX,self.nY,self.nZ)).astype('float')

        self.u = data['u'] ; self.v = data['v'] ; self.w = data['w']
        self.zBot = 0.0
        self.target_uHub = target_uHub
        self.target_TI = target_TI
        self.zHub = zHub
        self.update()
        self.generateGrid()
        self.T = T
        self.dT = self.T / self.nX

    def readMannBoxSameAsAbove(self,pathToBinaries,zHub=99.3):
        """
        I was cross-checking whether a different method of reading the binary file would yield the same data. It does.

        Parameters
        ----------
        pathToBinaries : string,
            where to find binary outputs
        zHub : float,
            optional argument used to identify which point represents hub location, needed once a grid is generated
        """
        self.readMannTxt(pathToBinaries)
        components  = ['u','v','w']
        data = {}
        nNumbers = self.nX*self.nY*self.nZ
        for component in components:
            componentPath = glob.glob(os.path.join(pathToBinaries,
             '{0}_{1}.bin'.format(self.prefix,component)))
            if len(componentPath)==0:
                raise Exception("Could not find file for component {0}".format(component))
            with open(componentPath[0], mode='rb') as file:
                nBytes = os.path.getsize(componentPath[0])
                if (nBytes/4.0 != nNumbers):
                    raise Exception("File size does not match expectations")
                fileContent = file.read()

                data[component] = np.zeros((self.nZ,self.nY,self.nX),float)
                byteStart = 0
                for ix in range(self.nX):
                    for iy in range(self.nY):
                        for iz in range(self.nZ):
                            data[component][iz,iy,ix] = struct.unpack('f',fileContent[byteStart:byteStart+4])[0]
                            byteStart += 4

        for component in components:
            data[component] = np.einsum('kji->ijk',data[component])
        self.u = data['u'] ; self.v = data['v'] ; self.w = data['w']
        self.U = np.sqrt(self.u**2+self.v**2)
        self.zBot = 0.0
        self.zHub = zHub
        self.generateGrid()

    def readBTS(self,pathToBinary):
        """
        Read TurbSim Binary FF.

        Parameters
        ----------
        pathToBinaries : string,
            where to find binary outputs
        """
        fPath = os.path.join(pathToBinary,"{0}.bts".format(self.prefix))

        filePath = glob.glob(fPath)

        if len(filePath)==0:
            raise Exception("Could not find file at {0}.".format(fPath))

        print ("Opening file {0}...".format(filePath[0]))
        self.filePath = filePath
        self.fileDir  = pathToBinary
        components = ['u','v','w']

        with open(filePath[0], mode='rb') as file:
            fileContent = file.read()

            self.nZ, self.nY, self.nTower, self.nTimeSteps = \
                                    struct.unpack('i'*4,fileContent[2:18])
            self.dZ, self.dY, self.dT, self.uHub, self.zHub, self.zBot = \
                                    struct.unpack('f'*6,fileContent[18:42])
            self.nSeconds = self.nTimeSteps * self.dT

            vSlope = {} ; vIntercept = {} ; byteStart = 42
            for component in components:
                vSlope[component], vIntercept[component] = struct.unpack('f'*2,
                          fileContent[byteStart:byteStart+8])
                byteStart += 8

            nChar = struct.unpack('i',fileContent[byteStart:byteStart+4])
            byteStart += 4

            vNumber = ""
            for i in range(int(nChar[0])):
                vNumber += str(chr(struct.unpack('B',fileContent[byteStart:byteStart+1])[0]))
                byteStart += 1
            self.info = vNumber

            data = np.zeros((3,self.nY,self.nZ,self.nTimeSteps),float)
            for it in range(self.nTimeSteps):
                for iz in range(self.nZ):
                    for iy in range(self.nY):
                        for iComponent,component in enumerate(components):
                            data[iComponent,iy,iz,it] = struct.unpack('h',fileContent[byteStart:byteStart+2])[0]
                            byteStart += 2

            for iComponent,component in enumerate(components):
                data[iComponent,:,:,:] = (data[iComponent,:,:,:] - vIntercept[component])/vSlope[component]

            data = np.einsum('ljki->lijk',data)
            self.u = data[0,:,:,:]
            self.v = data[1,:,:,:]
            self.w = data[2,:,:,:]
            self.U = np.sqrt(self.u**2+self.v**2)

            if self.nTower>0:
                dataTower = np.zeros((self.nTower,self.nZ,self.nTimeSteps),float)
                for iz in range(self.nTower):
                    for iComponent,component in enumerate(components):
                        dataTower[iComponent,iz,it] = struct.unpack('h',fileContent[byteStart:byteStart+2])[0]
                        byteStart += 2
                for iComponent,component in enumerate(components):
                    dataTower[iComponent,:,:] = (dataTower[iComponent,:,:] - vIntercept[component])/vSlope[component]

        self.generateGrid()

    def stdDict(self):
        """
        Calculate standard deviation (over entire time series) of each component.

        Returns
        -------
        stdDictOut : dict,
           keys are 'u','v','w' and give standard deviation at approximate hub location
        """
        stdDictOut = {}
        for component in ['u','v','w','U']:
            stdDictOut[component] = np.std(getattr(self,component)[:,self.jHub,self.kHub])
        return stdDictOut

    def scaleStd(self,stdDictIn,verbose=False):
        """
        In case the Mann box does not have the desired std. dev. (to eventually have the correct TI), then this function can be used to scale it based on desired values.

        Parameters
        ----------
        stdDictIn : dict,
            keys are 'u','v','w' and give target standard deviation for each component
        """
        for component in stdDictIn.keys():
            current = getattr(self,component)
            oldStd = np.std(current[:,self.jHub,self.kHub])
            if verbose:
                print ("Scaling Std. Dev. of {0}: {1:.2f} to {2:.2f}".format(component,oldStd,stdDictIn[component]))
            setattr(self,component, current * stdDictIn[component] / oldStd)
        self.update()

    def update(self):
        self.U = np.sqrt(self.u**2+self.v**2)

    def scaleU(self,targetU,alpha=None,zRef=99.3):
        """
        Parameters
        ----------

        targetU : float,
            desired mean of 'u' component wind at approximate hub location

        alpha : float,
            power law exponent

        zRef : float,
            height of targetU
        """
        if alpha is not None:
            u = power_law(self.z,zRef,targetU,alpha)
        else:
            u = targetU
        offset = u - np.mean(self.u[:,self.jHub,self.kHub])
        self.u = self.u + offset
        self.update()

    def generateGrid(self):
        """
        Generates mesh attributes.
        """
        self.y = np.array([ 0 + i*self.dY for i in range(self.nY) ])
        self.z = np.array([ self.zBot + i*self.dZ for i in range(self.nZ) ])
        [self.Y,self.Z] = np.meshgrid(self.y,self.z)
        self.kHub = int(self.z2k(self.zHub))
        self.yHub = int(np.mean(self.y))
        self.jHub = int(self.nY/2)
        self.tiHub = self.TI(j=self.jHub,k=self.kHub)
        self.getRotorPoints()

    def y2j(self,y):
        """
        Computes j grid index for a given y.
        """
        return np.argmin(np.abs(self.y-y))

    def z2k(self,z):
        """
        Computes k grid index for a given z.
        """
        return np.argmin(np.abs(self.z-z))

    def TI(self,y=None,z=None,j=None,k=None):
        """
        If no argument is given, compute TI over entire grid and return array of size (nY,nZ). Else, compute TI at the specified point.

        Parameters
        ----------
        y : float,
            cross-stream position [m]
        z : float,
            vertical position AGL [m]
        j : int,
            grid index along cross-stream
        k : int,
            grid index along vertical
        """
        if ((y==None) & (j==None)):
            return np.std(self.U,axis=0) / np.mean(self.U,axis=0)
        if ((y==None) & (j!=None)):
            return (np.std(self.U[:,j,k])/np.mean(self.U[:,j,k]))
        if ((y!=None) & (j==None)):
            uSeries = self.U[:,self.y2j(y),self.z2k(z)]
            return np.std(uSeries)/np.mean(uSeries)

    def visualize(self,component='U',time=0):
        """
        Quick peak at the data for a given component, at a specific time.
        """
        data    = getattr(self,component)[time,:,:]
        plt.figure() ;
        plt.imshow(data) ;
        plt.colorbar()
        plt.show()

    def spectrum(self,component='u',y=None,z=None):
        """
        Calculate spectrum of a specific component, given time series at ~ hub.

        Parameters
        ----------
        component : string,
            which component to use
        y : float,
            y coordinate [m] of specific location
        z : float,
            z coordinate [m] of specific location

        """
        if y==None:
            k = self.kHub
            j = self.jHub
        data    = getattr(self,component)
        data    = data[:,j,k]
        N       = data.size
        freqs   = fftpack.fftfreq(N,self.dT)[1:N/2]
        psd     = (np.abs(fftpack.fft(data,N)[1:N/2]))**2
        return [freqs, psd]

    def readFastOut(self,pathToFastOut=None):
        """
        Parameters
        ----------
        pathToFastOut : str,
            complete absolute path to the *.out file produced by FAST v8
        """
        if pathToFastOut is None:
            listFiles = glob.glob(os.path.join(self.fileDir,'*.out'))
            if len(listFiles)!=1:
                'Check what is in {0} and try again'.format(self.fileDir)
            else:
                pathToFastOut = listFiles[0]
        self.fastOut = read_csv(pathToFastOut, delim_whitespace=True,
                                skiprows=[0,1,2,3,4,5,7], index_col=[0])
        self.fastOut["U"] = np.sqrt(self.fastOut.Wind1VelX**2+self.fastOut.Wind1VelY**2)

    def getRotorPoints(self):
        """
        In the square y-z slice, return which points are at the edge of the rotor in the horizontal and vertical directions.

        Returns
        -------
        jLeft : int,
            index for grid point that matches the left side of the rotor (when looking towards upstream)
        jRight : int,
            index for grid point that matches the right side of the rotor (when looking towards upstream)
        kBot : int,
            index for grid point that matches the bottom of the rotor
        kTop : int,
            index for grid point that matches the top of the rotor
        """
        self.zBotRotor      = self.zHub - self.R
        self.zTopRotor      = self.zHub + self.R
        self.yLeftRotor     = self.yHub - self.R
        self.yRightRotor    = self.yHub + self.R
        self.jLeftRotor  = self.y2j(self.yLeftRotor)
        self.jRightRotor = self.y2j(self.yRightRotor)
        self.kBotRotor   = self.z2k(self.zBotRotor)
        self.kTopRotor   = self.z2k(self.zTopRotor)
