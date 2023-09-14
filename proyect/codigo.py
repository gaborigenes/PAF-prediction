import wfdb
from pyrqa.time_series import TimeSeries
from pyrqa.settings import Settings
from pyrqa.analysis_type import Classic
from pyrqa.neighbourhood import FixedRadius
from pyrqa.metric import EuclideanMetric
from pyrqa.computation import RQAComputation
from pyrqa.computation import RPComputation
from pyrqa.image_generator import ImageGenerator
import sys
import numpy as np
import matplotlib.pyplot as plt
from wfdb.processing import gqrs_detect
import pyhrv.time_domain as td
import pyhrv.nonlinear as nl
import pyhrv.frequency_domain as fd


def sample(ecg, duracion, inicio, fin, n):
    #print(inicio, fin)
    registro = wfdb.rdrecord('Tesis/proyect/archivos-CSV/mit-bih-paf-database/p0'+str(ecg),  sampfrom = inicio, sampto= fin, channels = [0])
   
    #extraer muestra de 60 segundos
    longitud = duracion*registro.fs
    sample = registro.p_signal[:longitud]
    muestra = sample.flatten()
    fin = (n+1)*duracion*registro.fs

    '''
    #graficar
    fig , ax = plt.subplots()
    x = np.linspace(0, len(muestra), num=10)
    x_val = [int(x[i]/registro.fs) for i in range(len(x))]
    plt.xticks(x, x_val)
    plt.plot(muestra)
    '''
    #print('nuevo fin', fin)
   
    return muestra , fin, registro.fs
def graf(ecg , peaks, n_de_muestras):
    fig , ax = plt.subplots()
    fig.set_dpi(100)
   
    if n_de_muestras < 6:
        n_de_muestras = n_de_muestras
    else:
        n_de_muestras = 6
       
    for i in range(n_de_muestras):
           
        ax1=fig.add_subplot(6,1,i+1)
        ax1.plot(ecg[i],  label='ECG Signal')
        ax1.set_xticks(np.arange(0, len(ecg), 10))
        plt.scatter(peaks[i], ecg[i][peaks[i]]+1.5, color='red', label='RR Peaks')
        #plt.vlines(peaks[i], ymin=min(ecg[i])+0.2, ymax= max(ecg[i])-0.2, color = 'r', label = 'axvline - full height')
        ax1.set_title('Original ECG')
        ax1.set_xlabel('Time (ms)')
        ax1.set_ylabel('mV')

for i in range(1,2):
    signal = i
    print('#: ',signal)
    duracion = 900
    np.set_printoptions(suppress=True, threshold = sys.maxsize, linewidth = 500)
    aux = wfdb.rdrecord('Tesis/proyect/archivos-CSV/mit-bih-paf-database/p0'+str(signal))
    tiempo_total = len(aux.p_signal)/aux.fs
    n_de_muestras = 1 #int(tiempo_total/duracion)
    print('cantidad de muestras: ', n_de_muestras)
       
       
   
    ecg = np.zeros((n_de_muestras, duracion*aux.fs))
    peaks = []
    RR = []
    
    inicio = 0
    fin = duracion*aux.fs
    for n in range(n_de_muestras):
        ecg[n, :] , fin , fs = sample(signal, duracion, inicio , fin, n+1)
        inicio = int(fin - duracion*aux.fs)
        qrs_inds = gqrs_detect(ecg[n,:], fs=fs)
        peaks.append(qrs_inds)
           
    #print(ecg)  
    for i in range(n_de_muestras):
        RR.append( np.diff(peaks[i]) / fs  )
   
    #print(RR)
    graf(ecg, peaks, n_de_muestras)
    

    
  

    desviacion_absoluta_media = np.zeros(2)
    
    DMA = []
    for i in range(n_de_muestras):
        RR_mean = np.round(np.mean(RR[i]), decimals=3)
       
        for j in range(len(RR[i])):
            resta = RR[i][j]-RR_mean
            DMA.append((np.abs(resta)))
        desviacion_absoluta_media[i] = np.round(np.sum(DMA)/len(RR[i]), decimals=3)
        
    
        


    vector = np.zeros((n_de_muestras,21))
       
    for i in range(n_de_muestras):
        #Variables temporales
        RR_mean = np.round(np.mean(RR[i]), decimals=3)
        rmssd = np.round(td.rmssd(RR[i]), decimals = 3)
        SDSD = np.round(td.sdsd(RR[i])  , decimals = 3)
        SDRR = np.round(td.sdnn(RR[i])  , decimals = 3)
        nn_50 = np.round(td.nn50(RR[i]) , decimals = 3)
        #desviacion_absoluta_media = np.round(np.sum(np.abs(RR[i] - RR_mean))/len(RR[0]) , decimals = 3)
       
       
        #variables no lineales
        sd = nl.poincare(RR[i] , show=False)
        sd1 = np.round(sd[1], decimals = 3)
        sd2 = np.round(sd[2], decimals = 3)
       
        #variables de frecuencia
        welch = fd.welch_psd(RR[i], show=True)
        AR = fd.ar_psd(RR[i], show = False)
        abs_psd = np.round(welch['fft_abs'], decimals = 3)
        LF_HF = np.round(welch['fft_ratio'], decimals = 3)
        ar_psd = np.round(AR['ar_abs'], decimals = 3)
        LF_HF_ar = np.round(AR['ar_ratio'], decimals = 3)
        
        
        #recurrence analisis
        time_series = TimeSeries(RR[i],
                             embedding_dimension= 7,
                             time_delay = 2)
   
        settings = Settings(time_series,
                        analysis_type=Classic,
                        neighbourhood=FixedRadius(0.65),
                        similarity_measure=EuclideanMetric,
                        theiler_corrector=1)
        computation = RQAComputation.create(settings,verbose =False)
   
        result = computation.run()
        result.min_diagonal_line_length = 2
        result.min_vertical_line_length = 2
        result.min_white_vertical_line_length = 2
        
        rec_rate = np.round(result.recurrence_rate, decimals = 3)
        l_min = np.round(result.average_diagonal_line, decimals = 3)
        l_max = np.round(result.longest_diagonal_line, decimals = 3)
        entropia = np.round(result.entropy_diagonal_lines, decimals = 3)
        trapping_time = np.round(result.trapping_time, decimals = 3)
       
       
           
        vector[i] = [RR_mean,
                     SDRR[0],
                     SDSD[0],
                     rmssd[0],
                     nn_50[1],
                     sd1,
                     sd2,
                     desviacion_absoluta_media[i],
                     abs_psd[0],
                     abs_psd[1],
                     abs_psd[2],
                     LF_HF,
                     ar_psd[0],
                     ar_psd[1],
                     ar_psd[2],
                     LF_HF_ar,
                     rec_rate, 
                     l_min,
                     l_max,
                     entropia,
                     trapping_time]
        
        
                     
                     
       

    matriz = str(vector).replace(' [', '').replace('[', '').replace(']', '')
    print(matriz)
  
    with open('Tesis/proyect/base_predictor/lejanos/p'+str(signal)+'.txt', 'w') as file :
        file.write(str(matriz) + '\n')
        print('n'+str(signal)+'.txt creado')
        file.close()
