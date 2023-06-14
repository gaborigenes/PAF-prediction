# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 16:42:38 2023

@author: Ricardo j Quintero
"""
import wfdb
import numpy as np
import pandas as pd
import pt_QRS


archivo = f'mit-bih-arrhythmia-database-1.0.0/115'
registro = wfdb.rdrecord(archivo, sampfrom=0, sampto = 3000)
anotacion = wfdb.rdann(archivo, 'atr', sampfrom=0, sampto=3000, shift_samps=True)
wfdb.plot_wfdb(record=registro, time_units='seconds', figsize=(15,8) )

fs = anotacion.fs     #anoto el valor de la frecuencia
#convierto señal ecg en un dataframe de pandas, siendo esta el registro transpuesto
#con columnas timestamp y el valor del ecg
#dataframe(ndarray, columns)

ecg = pd.DataFrame(np.array([list(range(len(registro.adc())))  ,  registro.adc()[:,0]] ).T,columns=['TimeStamp','ecg'])
print((ecg))

ecg_integrada = pt_QRS.solve(ecg,fs)