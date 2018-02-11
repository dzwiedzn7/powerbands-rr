import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import boxcar,blackman,hamming,hann,hanning,resample,periodogram


def anotations(signal,flags):

    """
    replace values from original signal with flag anotation,with nan values
    """
    anotation = np.where(flags==2)
    signal[anotation] = np.nan
    return signal



def interpolation(signal):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(signal.shape[0])
    good = np.where(np.isfinite(signal))
    f = interp1d(inds[good], signal[good],kind='linear',bounds_error=False)
    out = np.where(np.isfinite(signal),signal,f(inds))
    return out


def resampling(signal,freq,kind='cubic',resolution = 2048):
    """
    resampling of linear interpolated singal in order to get signal of desired frequency
    """
    x = np.arange(len(signal))
    new_length = 60 * freq
    new_x = np.linspace(x.min(), x.max(), new_length) #start,stop,resolution=378
    out = interp1d(x, signal,kind)(new_x)

    #nextPower = np.ceil(np.log2(new_length))
    #deficit = int(np.power(2, nextPower) - len(new_y))
    #out = np.zeros(deficit+len(new_y),dtype=new_y.dtype)
    #out[:len(new_y)] = new_y
    return out

def zero_padding(signal,resolution=2048):
    """
    zero padd right end of signal to obtaing signal of length 2**n to perform fft
    """
    out = np.zeros(resolution,dtype=signal.dtype)
    out[:len(signal)] = signal
    return out




def signal_fft_plot(signal,freq,window_func):
    """
    plot zero centered fast fourier transform of signal
    """
    N = signal.shape[0]
    window = window_func(len(resampled))
    frq = np.fft.fftfreq(len(signal),1./freq)
    frq = frq[range(N/2)]
    Y = np.fft.fft(signal*window)/N
    Y = Y[range(N/2)]
    plt.title("Frequency Spectrum of Heart Rate Variability")
    plt.plot(frq, abs(Y)) #Plot it
    plt.xlabel("Frequencies in Hz")
    plt.show()
    #plt.savefig("FFT_Signal_hz" + str(freq) + )









if __name__ == "__main__":
    windows =[boxcar,blackman,hamming,hann,hanning] #chosed window function,we can add more
    frequencies = [2,3,4,5,6,7] #chosed frequencies,we can add more

    data_files = glob.glob('*.rea')
    result = np.zeros([len(data_files),len(frequencies),len(windows),3])

    """
    here are three for loops,one for every window function,one for every frequencies and one foer every data file
    """
    for name in data_files:
        signal = np.loadtxt(name,skiprows=1, unpack=True)[1] * 0.001 #convert ms to seconds
        flags = np.loadtxt(name,skiprows=1, unpack=True)[2]
        nan = anotations(signal,flags)
        inter = interpolation(nan)

        for  idx,freq in enumerate(frequencies):
            resampled = resampling(inter,freq)
            #signal_fft_plot(resampled,freq)  #just to show example of how to use it
            total_variance = np.std(resampled)
            resampled = zero_padding(resampled)

            for jdx,window_func in enumerate(windows):
                N = resampled.shape[0]
                window = window_func(len(resampled))
                frq = np.fft.fftfreq(len(resampled),1./freq) #becouse retarded integer division
                frq = frq[range(N/2)]
                fourier = np.fft.fft(resampled*window)/N
                fourier = fourier[range(N/2)]
                PSD = sum(abs(fourier)**2)/len(resampled)
                print "PSD: ",PSD,"VAR: ",total_variance #the total power of the periodogram should be equal the variance of the original rr sequence
                """
                below it will calculate powerbands for very low,low and hight frequencies.Wasn't sure if I should wrap it in new function/method
                """
                vlf = np.trapz(abs(fourier[(frq<=0.04)]))
                print "VLF:", vlf
                lf = np.trapz(abs(fourier[(frq>=0.04) & (frq<=0.15)]))
                print "LF:", lf
                hf = np.trapz(abs(fourier[(frq>=0.16) & (frq<=0.5)]))
                print "HF:", hf
                result[idx,jdx] = vlf,lf,hf
    result = np.sum(result,axis=0)/len(data_files)
    print result,result.shape
    np.savetxt("partial2.csv", result, delimiter=",",fmt='%s')
    """
    conver to latex table with tably.py,command$ python3 tably.py [filename.csv] > [filename.tex]
    """
