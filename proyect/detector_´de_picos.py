# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 16:42:38 2023

@author: Ricardo j Quintero
"""
import wfdb
import numpy as np
import pandas as pd
import pt_QRS2_numpy
import time

start = time.time()
numero_signal = 16272
archivo =  f'archivos-CSV/files/MIT-BIH-NSR_760_MUESTRAS/16265_muestra_40'  #f'mit-bih-af-database/0{str(numero_signal)}'
registro = wfdb.rdrecord(archivo, sampfrom=0 )
#anotacion = wfdb.rdann(archivo, 'dat', sampfrom=0, sampto=45*250, shift_samps=True)
wfdb.plot_wfdb(record=registro, time_units='seconds', figsize=(15,8) )

fs = 0.35*128     #anoto el valor de la frecuencia
#convierto señal ecg en un dataframe de pandas, siendo esta el registro transpuesto
#con columnas timestamp y el valor del ecg
#dataframe(ndarray, columns)

ecg = pd.DataFrame(np.array([list(range(len(registro.adc())))  ,  registro.adc()[:,0]] ).T,columns=['TimeStamp','ecg'])
print((ecg))

ecg_integrada = pt_QRS2_numpy.solve(ecg,fs,numero_signal)
#print(ecg_integrada)

end = time.time()
print(end-start,' s')