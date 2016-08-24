import numpy as np
from numpy import fft
import pandas as pd
# from scipy import fftpack
from scipy import interpolate
from scipy import integrate
from scipy import signal
import matplotlib.pyplot as plt
# import os

def dataload_lanl_csv(osc_file):
    """ Load csv file for demodulation.
    Returns dictionary containing timebase info, oscillations, and field.
    """
    df_tdo = pd.read_csv(osc_file, '\t', usecols=[5, 9], skiprows=4)
    df_time = pd.read_csv(osc_file, '\t', usecols=[5, 9], skiprows=2, nrows=1, header=None)
    osc_dict = {'file': osc_file, 'dt_fast': df_time[5][0], 'tdo': df_tdo['PDO']}
    osc_dict.update({'dt_slow': df_time[9][0], 'Field': df_tdo['Field_fixed']})
    return osc_dict

def dataload_spf_bin(osc_file):
    """ Load byte data taken w/ Labview on small pulsed field.
    Assumes text file w/ scaling factors named osc_file+'Scale.txt' in same directory as osc_file.
    Returns dictionary containing timebase info, oscillations, and field.
    """
    scale_file = osc_file+'Scale.txt'
    raw_data = np.fromfile(osc_file, dtype=np.uint8) # load byte data as uint8
    osc_data = raw_data[:len(raw_data)//2] # first half of file is tdo signal
    pu_data = raw_data[len(raw_data)//2:] # second half of file is pickup voltage
    scale_data = pd.read_csv(scale_file, '\t') # load scaling data
    osc = scale_data.yincrement[0]*(osc_data-128.)+scale_data.yorigin[0] # convert to Voltage
    pu = scale_data.yincrement[1]*(pu_data-128.)+scale_data.yorigin[1]
    df_osc = pd.DataFrame({'tdo': osc, 'pu': pu})
    osc_dict = {'file': osc_file, 'dt': scale_data.dt[0], 'df_osc': df_osc}
    return osc_dict

def pu_to_field(osc_dict, Nturns, diameter, plot=False):
    """ Integrates pickup signal to obtain field. Takes output of dataload_spf_bin as input.
    Nturns: number of turns in pickup coil.
    diameter: of pickup coil in meters.
    """
    coilArea = Nturns*np.pi*(diameter/2)**2
    osc_dict['df_osc'].pu -= np.mean(osc_dict['df_osc'].pu[:100000]) # subtract offset
    field = (1/coilArea)*integrate.cumtrapz(osc_dict['df_osc'].pu, x=None, dx=osc_dict['dt'])
    # data is oversampled, so decimate it:
    field = signal.decimate(field, 40)
    tdo = signal.decimate(osc_dict['df_osc'].tdo, 10)
    # get new time steps after decimation:
    dt_fast = osc_dict['dt']*osc_dict['df_osc'].tdo.size/len(tdo)
    dt_slow = osc_dict['dt']*osc_dict['df_osc'].pu.size/len(field)
    osc_dict = {'file': osc_dict['file'], 'dt_fast': dt_fast}
    osc_dict.update({'tdo': tdo, 'dt_slow': dt_slow, 'Field': field})
    if plot:
        plt.plot(field)
        plt.xlabel('Index')
        plt.ylabel('Field (T)')
    return osc_dict


def demod(osc_dict, L=None, nstep=None, plot=True):
    """ Performs demodulation of data loaded via dataload_lanl_csv or dataload_spf_bin+pu_to_field.
    """
    dt_fast = osc_dict['dt_fast']
    dt_slow = osc_dict['dt_slow']
    tdo = np.array(osc_dict['tdo'])
    field = np.array(osc_dict['Field'])
    file = osc_dict['file']
    npnts = tdo.size # number of points in tdo signal
    npnts_slow = int(npnts*dt_fast/dt_slow) # number of points in everything else
    time_fast = np.arange(0, dt_fast*npnts, dt_fast) # timebase for tdo
    time_slow = np.arange(0, dt_slow*npnts_slow, dt_slow) # timebase for everything else
    Fs = 1/dt_fast # sampling frequency
    if L is None:
        init_fft = fft.rfft(tdo[:1024])[:-1]
        f = Fs*np.arange(0, 512)/(512)
        fmax = f[np.argmax(abs(init_fft)[1:])]
        pts_per_T = int(2/(fmax*dt_fast))
        L = 30*pts_per_T
        if nstep is None:
            nstep = 20*pts_per_T
    npad = np.power(2, 16)-L # zero-padding (for speed should be 2^N-L, but probably doesn't matter)
    pad = np.zeros(npad)
    N = int(npnts/nstep)
    freq = np.empty(N)
    amp = np.empty(N)
    time = np.empty(N)
    # Perform sliding FFT
    for i in range(N-nstep):
        fftdata = np.concatenate([tdo[i*nstep:i*nstep+L], pad])
        #tdo_fft = fft.rfft(df_tdo.tdo[i*nstep:i*nstep+L]) # if no zero-padding
        #f = Fs*np.arange(0,L)/L; # if no zero-padding
        f = Fs*np.arange(1000, L+npad)/(L+npad)
        freq[i] = f[np.argmax(abs(fft.rfft(fftdata)[1000:]/L))]
        amp[i] = np.sqrt(2*np.mean(np.square(fftdata[:L])))
        time[i] = dt_fast*nstep*i
    time = time[:-nstep]
    freq = freq[:-nstep]
    amp = amp[:-nstep]
    demod_dict = {'file': file, 'time_tdo': time_fast, 'tdo': tdo}
    demod_dict.update({'time_Field': time_slow, 'Field': field[:len(time_slow)]})
    demod_dict.update({'time_freq': time, 'freq': freq, 'amp': amp})
    if plot:
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(time, freq*1e-6, 'b', label='Frequency')
        ax2.plot(time, amp*1e3, 'r', label='Amplitude')
        ax2.set_ylabel('Amplitude (mV)')
        ax1.set_ylabel('Frequency (MHz)')
        ax1.set_xlabel('Time (s)')
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc=0)
        ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
        ax1.patch.set_visible(False) # hide the 'canvas'
        plt.show()
    return demod_dict

def spline_plot(demod_dict, plot=True, updown=True, delta=True, amp=True):
    """ Performs cubic spline and plots freq vs. field
    Takes output of demod as input.
    """
    ncs = 32768//2
    t_max = np.min([demod_dict['time_freq'][-1], demod_dict['time_Field'][-1]])
    time_cs = np.linspace(0, t_max, ncs)
    # find spline representation
    tck = interpolate.splrep(demod_dict['time_Field'], demod_dict['Field'])
    # evaluate spline representation
    field_cs = interpolate.splev(time_cs, tck, der=0)
    if amp:
        tck = interpolate.splrep(demod_dict['time_freq'], demod_dict['amp'])
        amp_cs = interpolate.splev(time_cs, tck, der=0)
    freq = demod_dict['freq']
    if delta:
        freq -= freq[0]
    tck = interpolate.splrep(demod_dict['time_freq'], freq) # find spline representation
    freq_cs = interpolate.splev(time_cs, tck, der=0) # evaluate spline representation
    if amp:
        df_spline = pd.DataFrame({'time': time_cs, 'Field': field_cs,
                                  'freq': freq_cs, 'amp': amp_cs})
        df_spline = df_spline[['time', 'Field', 'freq', 'amp']]
    else:
        df_spline = pd.DataFrame({'time': time_cs, 'Field': field_cs, 'freq': freq_cs})
        df_spline = df_spline[['time', 'Field', 'freq']]
    if plot and updown:
        maxidx = df_spline.Field.argmax()
        if not delta:
            plt.plot(df_spline.Field.iloc[:maxidx], df_spline.freq.iloc[:maxidx]*1e-6,
                     label='Up sweep')
            plt.plot(df_spline.Field.iloc[maxidx:], df_spline.freq.iloc[maxidx:]*1e-6,
                     label='Down sweep')
            plt.ylabel('Frequency (MHz)')
        else:
            plt.plot(df_spline.Field.iloc[:maxidx], df_spline.freq.iloc[:maxidx]*1e-3,
                     label='Up sweep')
            plt.plot(df_spline.Field.iloc[maxidx:], df_spline.freq.iloc[maxidx:]*1e-3,
                     label='Down sweep')
            plt.ylabel('$\Delta$Freq. (kHz)')
        plt.xlabel('Field (T)')
        ax = plt.gca()
        ax.set_xlim(-0.1)
        plt.legend(loc=0)
        plt.show()
    if plot and not updown:
        if not delta:
            plt.plot(df_spline.Field, df_spline.freq*1e-6, label='Full sweep')
            plt.ylabel('Frequency (MHz)')
        else:
            plt.plot(df_spline.Field, df_spline.freq*1e-3, label='Full sweep')
            plt.ylabel('$\Delta$Freq. (kHz)')
        plt.xlabel('Field (T)')
        ax = plt.gca()
        ax.set_xlim(-0.1)
        plt.legend(loc=0)
        plt.show()
    if plot and amp:
        if updown:
            maxidx = df_spline.Field.argmax()
            plt.plot(df_spline.Field.iloc[:maxidx], df_spline.amp.iloc[:maxidx]*1e3,
                     label='Up sweep')
            plt.plot(df_spline.Field.iloc[maxidx:], df_spline.amp.iloc[maxidx:]*1e3,
                     label='Down sweep')
            plt.ylabel('Amplitude (mV)')
        plt.xlabel('Field (T)')
        ax = plt.gca()
        ax.set_xlim(-0.1)
        plt.legend(loc=0)
        plt.show()
    demod_dict.update({'df_spline': df_spline})
    return demod_dict

def data_save(demod_dict, file_root=None, spline=True, tdo_pu=False):
    """ Saves .txt files. Options:
    -spline (default): Saves file_root_demod.txt containing time_cs, field_cs, freq_cs
    -tdo_pu: Saves file_root_raw.txt containing tdo signal (first row after header is time step dt)
    """
    if file_root is None:
        file_root = demod_dict['file']
        if file_root.endswith('.txt'):
            file_root = file_root[:-4]
    demod_file = file_root+'_demod.txt'
    raw_file = file_root+'_raw.txt'
    if spline:
        demod_dict['df_spline'].columns = ['time_s', 'Field_T', 'freq_Hz', 'amp_V']
        demod_dict['df_spline'].to_csv(demod_file, sep='\t', header=True, index=False)
    if tdo_pu:
        data_dict = dataload_spf_bin(demod_dict['file'])
        df_data = pd.DataFrame({'tdoV': data_dict['df_osc'].tdo, 'pickupV': data_dict['df_osc'].pu})
        df_data.loc[-1] = [data_dict['dt'], data_dict['dt']]
        df_data.index += 1
        df_data = df_data.sort_index()
        df_data.to_csv(raw_file, sep='\t', header=True, index=False)