# -*- coding: utf-8 -*-
"""
Created on Wed May 17 01:30:02 2023

@author: lapto thinkpad
"""

import numpy as np
import csv
import wfdb
import os
import sys
numero_ecg = 19830
ecg = str(numero_ecg)
with open('files/csv/NSR/'+ecg+'.csv', 'r') as f:
    data = list(csv.reader(f, delimiter=";"))

data = np.array(data[1:])
A = np.zeros(data.shape, dtype=int)
#print(data[:,1])

A[:,0] = [int(i) for i in data[:,0]]
A[:,1] = [int(i) for i in data[:,1]]


archivo = 'mit-bih-nsr-database/'+ecg
for i in range(len(A)):
    inicio = A[i,0]
    fin = A[i,1]+5506
    cuenta = str(i+1) 
    print(inicio,fin)
    #print('from: ',sampfrom , 'to: ', sampto)    
    
    registro, fields = wfdb.rdsamp(archivo, sampfrom=inicio, sampto=fin, )  

    file = 'files/MIT-BIH-NSR_760_MUESTRAS/'+ecg+'_muestra_'+cuenta
    
    existe = os.path.isfile(file + str('.dat'))
    #print(existe)
    if(existe == False):
        #print('no existe')
        wfdb.wrsamp(file, 250, units= ['mV','mV'], sig_name=['I', 'II'], p_signal=registro, fmt =['16','16'], adc_gain=[200,200], baseline=[0,0])
        print(file + str('.dat')+' creado')
    else:
        print('el archivo existe\nexit')
        sys.exit() 

    #ecg = pd.DataFrame(np.array([list(range(len(registro.adc()))),registro.adc()[:,0]]).T,columns=['TimeStamp','ecg'])
    
