# spectra.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 26 Aug 2024
#
# This is a module that contains methods to compute spectra.

# tukeyWindow
# hannWindow
# spectrum1D

import numpy as np



# Tukey window.
def tukeyWindow(t,alpha=0.5):
    import numpy as np

    nt = len(t)
    T = t[-1] - t[0]
    w = np.zeros((nt,))

    for i in range(nt):
        if ((t[i]-t[0]) < 0.5*alpha*T):
            w[i] = 0.5*(1.0-np.cos(2.0*np.pi*(t[i]-t[0])/(alpha*T)))
        elif (((t[i]-t[0]) >= 0.5*alpha*T) and ((t[i]-t[0]) <= T/2.0)):
            w[i] = 1.0
        else:
            w[i] = w[-(i+1)]

    return w





# Hann window.  It is just the Tukey window with an alpha parameter of 1.0
def hannWindow(t):
    import numpy as np
    w = tukeyWindow(t,1.0)

    return w
    



# This takes the 1D spectrum of a single series of data.  It can apply 
# a Tukey or Hann window to the data if desired.
def spectrum1D(t,f,windowType='tukey',tukeyAlpha=0.5):
    import numpy as np

    # Get the number of samples in the series
    nt = len(t)

    # Get the length of the series.
    T = t[-1] - t[0]

    # Compute the window function.
    if ((windowType == 'tukey') or (windowType == 'Tukey')):
        w = tukeyWindow(t,tukeyAlpha)
    elif ((windowType == 'hann') or (windowType == 'Hann')):
        w = hannWindow(t)
    else:
        w = np.ones((nt,))

    # Take the FFT of the windowed data and normalize by number of samples
    Sf = (np.sqrt(2.0)/nt) * np.fft.fft(w*f - np.mean(w*f))

    # Get the complex conjugate of the spectrum with itself, which should just
    # be the magnitude as a real number, but because of round-off, we are left
    # with a tiny imaginary part, so remove that with the 'real' function.
    Sf = np.real(Sf * np.conjugate(Sf))

    # Compute the frequency array.  It is (1/sample length) times wave number
    freq = (1.0/T) * np.linspace(0.0,float(nt),nt+1)

    # Only keep the half range of the data.
    Sf = Sf[0:int(nt/2)]
    freq = freq[0:int(nt/2)]

    # Compute the windowed function and return it to have access to the actual
    # time series input to FFT.
    g = w*f

    return Sf,freq,g
                        
    