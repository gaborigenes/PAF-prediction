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
import pyhrv.time_domain as td
import pyhrv.nonlinear as nl
#import math
import os

start = time.time()
numero_signal = 16265
archivo =  f'archivos-CSV/files/MIT-BIH-NSR/'+str(numero_signal)+'/'+str(numero_signal)+'_muestra_40'
registro = wfdb.rdrecord(archivo, sampfrom=0)
#anotacion = wfdb.rdann(archivo, 'dat', sampfrom=0, sampto=45*250, shift_samps=True)
#wfdb.plot_wfdb(record=registro, time_units='seconds', figsize=(15,8) )

fs = 0.35*128   #anoto el valor de la frecuencia
#convierto señal ecg en un dataframe de pandas, siendo esta el registro transpuesto
#con columnas timestamp y el valor del ecg
#dataframe(ndarray, columns)

ecg = pd.DataFrame(np.array([list(range(len(registro.adc())))  ,  registro.adc()[:,0]] ).T,columns=['TimeStamp','ecg'])
print((ecg))

op_x, op_y = pt_QRS2_numpy.solve(ecg,fs,numero_signal)

print('\n=================================================================\n')
  
#print(op_x)
RR = np.round(1000*np.diff(op_x)/128, decimals=2) #hago la diferencia entre las posiciones de los picos R
#print('RR(ms): ', RR)  #RR


'''
rr_n = RR[1:-1]
rr_n1 = RR[2:]
tupla = [(rr_n[i] , rr_n1[i]) for i in range(len(rr_n))]
li = [math.dist(tupla[i],tupla[i+1]) for i in range(len(tupla)-1)]
#print(li)
li_mean = (np.mean(li))
'''

RR_mean = np.round(np.mean(RR), decimals=3)
rmssd = np.round(td.rmssd(RR), decimals = 3)
SDSD = np.round(td.sdsd(RR)  , decimals = 3)
SDRR = np.round(td.sdnn(RR)  , decimals = 3)
nn_50 = np.round(td.nn50(RR) , decimals = 3)
entropy = np.round(nl.sample_entropy(RR), decimals=3)
desviacion_absoluta_media = np.round(np.sum(np.abs(RR - RR_mean))/len(RR) , decimals = 3)
#VLI = np.round(np.square(np.sum(np.square(li - li_mean)))/len(li) , decimals = 3)

#grafica de poincare
sd = nl.poincare(RR) 
sd1 = np.round(sd[1], decimals = 3)
sd2 = np.round(sd[2], decimals = 3)

vector_entrada = [RR_mean,SDRR[0],SDSD[0],rmssd[0],nn_50[1],sd1,sd2,entropy[0], desviacion_absoluta_media]


print('RR_mean: ', RR_mean)
print('sdnn: ', SDRR[0])
print('sdsd: ', SDSD[0])
print('rmssd: ', rmssd[0])
print('pnn50: ', nn_50[1])
print('sd1: ', sd1)
print('sd2: ',sd2)
print('entropia: ', entropy[0])
print('desviacion absoluta media: ', desviacion_absoluta_media)
#print('VLI: ', VLI)

file = open('vector_entrada.txt', 'w')
file.write(str(vector_entrada))
file.close()

end = time.time()
#print(end-start,' s')


























'''
for i in range(len(op_x)-8):
    
    RR_average = np.mean(RR1[i:i+8]) #saco la media de los ultimos 8 picos
    print('RR1_average',i,': ', RR_average)
    
    RR_average2 = RR_average
    RR_low_limit = 0.92*RR_average2
    RR_high_limit = 1.16*RR_average2
    #print(RR_low_limit, RR_high_limit)
     
    for i in range(1, len(RR1)):
        if RR_low_limit < RR1[i]<RR_high_limit:
            RR2.append(RR1[i])
            #print('RR2: ', RR2)   
            if(len(RR2) > 8):
                RR2.remove(RR2[0])
                RR_average2 = np.mean(RR2)
        
    print('RR2: ', RR2)
       
    RR_avg_op = 60*1000/(2*RR_average) #bpm
    if RR_average<=RR_high_limit and RR_average >= RR_low_limit: #If each of the eight most-recent sequential RR intervals that are calculated from RR AVERAGE1 is between the RR LOW LIMIT and the RR HIGH LIMIT,
        print("Normal sinus with average heart beat (bpm)", RR_avg_op)
    else:
        print("Not a normal sinus with average heart beat (bpm)", RR_avg_op)
        
'''  
'''
FC= 60000/RR

RR_mean = np.mean(RR)
SDNN = np.std(RR)
RMSSD = np.sqrt(np.mean(np.square(np.diff(RR)))/(len(RR)-1))
FC_mean =  np.mean(FC)



STD_FC = np.std(FC)

NN50 = np.sum(np.abs(np.diff(RR)) > 50)*1
pNN50 = ((NN50)/len(RR))*100
print('RR_mean:', RR_mean)

print('SDNN:', SDNN)

print('RMSSD:', RMSSD)
print('NN50: ',NN50)
print('pNN50: ',pNN50)
'''

