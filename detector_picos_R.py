# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 20:19:00 2023

@author: Ricardo j Quintero
"""



import wfdb
#WFDB es un paquete de python. Sirve para leer, escribir y procesar bases de datos en forma de onda (waveform database)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
#from scipy import signal as sg
import pan_tompkins_QRS 
from ritmo_cardiaco import*


archivo = 'paf-prediction-challenge-database-1.0.0/p04'
registro = wfdb.rdrecord(archivo, sampfrom=180, sampto=4000, )    
anotacion = wfdb.rdann(archivo, 'dat', sampfrom=180, sampto=4000,shift_samps=True)
#wfdb.plot_wfdb(record=registro, annotation=anotacion,figsize=(20,10))

#============================================================================================

fs = anotacion.fs     
ecg = pd.DataFrame(np.array([list(range(len(registro.adc()))),registro.adc()[:,0]]).T,columns=['TimeStamp','ecg'])
ecg_integrada = pan_tompkins_QRS.pan_tompkins_qrs(ecg,fs)


#===========================================================================================

#convierto ECG a un arreglo numpy
signal = ecg.iloc[:,1].to_numpy()

rc = ritmo_cardiaco(signal, fs)

result = rc.hallar_picos_r()
result = np.array(result)

# Clip the x locations less than 0 (Learning Phase)
result = result[result > 0]

# Calculate the heart rate
heartRate = (60*fs)/np.average(np.diff(result[1:]))
print("Heart Rate",heartRate, "BPM")

# Plotting the R peak locations in ECG signal
plt.figure(figsize = (20,4), dpi = 100)
plt.xticks(np.arange(0, len(signal)+1, 250))
plt.plot(signal, color = 'blue')        
plt.scatter(result, signal[result], color = 'red', s = 50, marker= '*')
plt.xlabel('Samples')
plt.ylabel('MLIImV')
plt.title("R Peak Locations")

wfdb.plot_wfdb(record=registro, annotation=anotacion,figsize=(20,10))