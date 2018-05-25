import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import boxcar,blackman,hamming,hann,hanning,resample,periodogram,kaiser,bartlett
from scipy.stats import kruskal
import gc

def table_prod_4freqrange(data,windows=1,frequencies=1):
    """
    this function will produce numpy array with psd of diffrent frequency bands
    of rr intervals signal with respect to windows function and frequency
    """
    import warnings
    warnings.filterwarnings("error")
    if windows == 1:
        windows =[boxcar,blackman,hamming,hann,bartlett] #chosed window function,we can add more
    if frequencies ==1:
        frequencies = [2.0,3.0,3.5,4.0,4.5,5.0,6.0,10.0] #chosed frequencies,we can add more
    result = np.zeros([len(data),len(frequencies),len(windows),3])
    couter = 0
    for zdx,name in enumerate(data):
        try:
            powers = 1
            data1 = 1
            del powers,data1
            data1 = Frequency_domain(name)
            for jdx,freq in enumerate(frequencies):
                for idx,win in enumerate(windows):
                    powers = data1.power_for_band(freq,win(data1.N))
                    if np.isnan(powers).any():
                        continue
                    else:
                        lol = [int(i) for i in powers]
                        result[zdx,jdx,idx] = lol
        except RuntimeWarning:
            gc.collect()
            couter += 1
            continue
    #print(couter)
    return result


def table_prod_4total(data,windows=1,frequencies=1):
    """
    this function will produce numpy array with total PSD
    of rr intervals signal with respect to windows function and frequency
    """
    import warnings
    warnings.filterwarnings("error")
    if windows == 1:
        windows =[boxcar,blackman,hamming,hann,bartlett] #chosed window function,we can add more
    if frequencies ==1:
        frequencies = [2.0,3.0,3.5,4.0,4.5,5.0,6.0,10.0] #chosed frequencies,we can add more
    result = np.zeros([len(data),len(frequencies),len(windows)])
    couter = 0
    for zdx,name in enumerate(data):
        try:
            powers = 1
            data1 = 1
            del powers,data1
            data1 = Frequency_domain(name)
            for jdx,freq in enumerate(frequencies):
                for idx,win in enumerate(windows):
                    powers = data1.fourier_windowing_energy_norm(win(data1.N))
                    powers = sum(abs(powers)**2)
                    if np.isnan(powers).any():
                        continue
                    else:
                        result[zdx,jdx,idx] = powers
        except RuntimeWarning:
            gc.collect()
            couter += 1
            continue
    #print(couter)
    return result




class Load():
    """
    here our signal is born

    read_signal read column of rr intervals [in ms] -our data
    read_flags read coulumn of annotations.Annotations inform us which rr intervals
    are 'unnormal' and need to be interpolated
    """
    def __init__(self,name):
        #self.data_files = glob.glob('*.rea')
        #self.one = self.data_files[0]
        self.one = name
        self.signal = self.read_signal(self.one)
        self.flags = self.read_flags(self.one)


    def read_signal(self,name):
        signal = np.loadtxt(name,skiprows=1, unpack=True)[1]# * 0.001 #convert ms to seconds
        return signal

    def read_flags(self,name):
        flags = np.loadtxt(name,skiprows=1, unpack=True)[2]
        return flags


class Processing(Load):
    """
    Here our signal is processed as long as it will be ready for FFT
    which mean, it will be resampled with constant frequency and zero padded
    becouse FFT requires input signal to be 2**N lenght
    """
    def __init__(self,name,freq):
        super(Processing, self).__init__(name)
        self.anon = Processing.anotations(self.signal,self.flags)
        self.inter = Processing.interpolation(self.anon)
        self.resampled = Processing.resampling(self.inter,freq)
        self.resampled_padded = Processing.zero_padding(self.resampled)

    @staticmethod
    def anotations(signal,flags):
        """
        replace values from original signal with flag anotation,with nan values
        """
        anotation = np.where(flags==2)
        signal[anotation] = np.nan
        return signal

    @staticmethod
    def interpolation(signal):
        '''
        interpolate to fill nan values
        '''
        inds = np.arange(signal.shape[0])
        good = np.where(np.isfinite(signal))
        f = interp1d(inds[good], signal[good],kind='linear',bounds_error=False)
        out = np.where(np.isfinite(signal),signal,f(inds))
        return out

    def add_one_by_one_gen(l):
        cumsum = 0
        for elt in l:
            cumsum += elt
            yield cumsum

    @staticmethod
    def resampling(signal,freq,kind='cubic',resolution = 2048):
        """
        resampling of linear interpolated singal in order to get signal of desired frequency
        """
        x = np.array([i for i in Processing.add_one_by_one_gen(signal)])
        new_length = 60 * freq
        new_x = np.linspace(x.min(), x.max(), new_length) #start,stop,resolution=378
        out = interp1d(x, signal,kind)(new_x)
        return out

    @staticmethod
    def zero_padding(signal,resolution=2048):
        """
        zero padd right end of signal to obtaing signal of length 2**n to perform fft
        """
        out = np.zeros(resolution,dtype=signal.dtype)
        out[:len(signal)] = signal
        return out



class Frequency_domain(Processing):
    """
    core class
    here our can be finnaly transormed to Frequency_domain w/o window function

    main point of this whole project was calculate PSD of given frequency band
    of rr intervals signal which is accomplished bellow
    """
    def __init__(self,name):
        super(Frequency_domain,self).__init__(name,freq=4)
        self.fsignal = self.resampled_padded
        self.N = self.fsignal.shape[0]

    def fourier_freq(self,freq):
        """
        additional method,easy way to obtain ordered x-axis values for FFT
        """
        self.freq = freq
        frq = np.fft.fftfreq(self.N, 1/freq)
        frq = frq[list(range(int(self.N / 2)))]
        return  frq

    def fourier_no_windowing(self):
        """
        simple FFT
        """
        bc = np.fft.fft(self.fsignal) / self.N
        bc = bc[list(range(int(self.N / 2)))]
        return  bc

    def fourier_windowing(self,window):
        """
        FFT with window
        """
        Y = np.fft.fft(self.fsignal*window) / self.N
        Y = Y[list(range(int(self.N / 2)))]
        return Y

    def fourier_windowing_energy_norm(self,window):
        """
        FFT windowed with compensation of energy loss due to windowing
        """
        bc = self.fourier_no_windowing()**2
        Y = self.fourier_windowing(window)**2
        VALID_FINAL = abs(Y) * (sum(abs(bc)) / sum(abs(Y)))
        return VALID_FINAL

    def power_for_band(self,freq,window):
        """
        calculate PSD for very low,low and high frequency band
        """
        frq = self.fourier_freq(freq)
        fourier = self.fourier_windowing_energy_norm(window)
        vlf = sum(fourier[(frq <= 0.04)] ** 2)
        lf = sum(fourier[(frq >= 0.04) & (frq <= 0.15)] ** 2)
        hf = sum(fourier[(frq >= 0.16) & (frq <= 0.5)] ** 2)
        return [vlf,lf,hf]

    def Parsheval_theorem(self,window):
        """
        important test!

        energy of signal before FFT must be equal energy of signal after FFT
        """
        total_variance = np.sum((window*self.fsignal)**2)
        fourier = np.fft.fft(self.fsignal*window)/np.sqrt(self.N)
        PSD = np.sum(np.abs(fourier)**2)
        print("Total variance is equal {} and PSD of signal is equal {}".format(total_variance,PSD))

    def Energy_loss_check(self,window):
        """
        2 important test!

        check if energy of FFT signal afrer windowing is equal energy of signal after
        """
        left = np.sum(abs((self.fourier_no_windowing()**2)))
        fourier = self.fourier_windowing_energy_norm(window)
        right = np.sum(abs(fourier))
        print("Energy of Boxcar window = {} and Energy of diffrent window = {}".format(left,right))
