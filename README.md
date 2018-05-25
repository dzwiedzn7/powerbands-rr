# powerbands-rr 
## python project for rr-intervals signal analysis in frequency domain

![alt text](https://github.com/dzwiedzn7/powerbands-rr/blob/master/Frequency%20of%206.0Hz.png)

core of this projects is **freqbands** package, containing methods for calculating power spectral density **PSD** for given frequency band:
1. very low frequency frq <= 0.04
2. low frequency frq >= 0.04 and frq <= 0.15
3. high frequency frq >= 0.16 and frq <= 0.5)
4. total PSD of full FFT signal

It contains two pre-made run files,running them will give you full analysis with hypothesis kruskal test and boxplots.

In maybe not so distant future this project will be extended to time-domain also.
