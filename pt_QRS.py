# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 17:08:08 2023

@author: Ricardo j Quintero
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def pasa_bajos(signal):
    
    y = signal.copy()  #copio la señal de entrada para trabajarla
    for n in signal.index: #itero por los indices del dataframe
        if(n<12):  #criterio de longitud minima de procesamiento
            continue
        #ecuación de diferencias de la función de transferencia
        ' '
        y.iloc[n,1] = 2*y.iloc[n-1,1] - y.iloc[n-2,1] + signal.iloc[n,1] - 2*signal.iloc[n-6,1] + signal.iloc[n-12,1]
    
    return y

def pasa_altos(signal):
    
    y = signal.copy()
    for n in signal.index:
        if(n<32):
            continue
        y.iloc[n,1] = y.iloc[n-1,1] - signal.iloc[n,1]/32 + signal.iloc[n-16,1] - signal.iloc[n-17,1] + signal.iloc[n-32,1]/32 #----> difference equation found in Biomedical Signal Analysis, By Rangaraj M. Rangayyan
        
    return y
    


def filtro_pasa_bandas(signal):
    
    signal_pasabajo = pasa_bajos(signal)
    signal_pasabandas = pasa_altos(signal_pasabajo)
    
    return signal_pasabandas

#================================================================================================================================================================================
    
def derivativa(signal):
    
    y = signal.copy()
    for n in signal.index:
        if (n<4):
            continue
        y.iloc[n, 1] = (- signal.iloc[n-2-2,1] - 2*signal.iloc[n-1-2,1] + 2*signal.iloc[n+1-2,1]  + signal.iloc[n+2-2,1])/8
    return y
    
#================================================================================================================================================================================
def cuadratica(signal):
    y = signal.copy() # creating a replica, equivalent to creating an array of same dimensions
    for n in signal.index:
      y.iloc[n,1] = y.iloc[n,1]**2
    return y

    

#================================================================================================================================================================================
def ventana_movil_integradora(signal):

    ventana = 30 # tamaño de la ventana
    media_ventana = int(ventana/2) # media ventana
    y = signal.copy() # copiando la señal
    
    for n in signal.index:

      suma = 0 # inicializo la suma, siendo 0 para cada n
      if(n > len(signal)-ventana): # no se puede procesar luego de este indice 
        break
      if(n < media_ventana): # indice menos que la mitad de la ventana
        continue
      for i in range(n-media_ventana,n+media_ventana+1): # i in range of one sample for each n
        suma = suma + signal.iloc[i,1] # equation as described in paper
      y.iloc[n,1] = suma/(media_ventana+1)	# dividing by N = number of samples in integration window

    return y




def solve(signal,fs):
    
    #recibo la señal tipo pandas dataframe
    signal_pasa_bandas = filtro_pasa_bandas(signal)
    signal_derivativa = derivativa(signal_pasa_bandas)
    signal_cuadratica = cuadratica(signal_derivativa)
    signal_integrada = ventana_movil_integradora(signal_cuadratica)
    
    '''
    metodología
    
    Se utiliza dos conjuntos de umbrales para detectar los complejos QRS. Un umbral para la señal ECG filtrada y el otro 
    para la señal producida por la ventana movil integradora
    Usando umbrales en dos señales, se mejora la confianza de la deteccion compara con usar una sola forma de onda, hay dos niveles
    de umbrales separados en cada uno de los dos conjuntos de umbrales, uno es la mitad de otro, los umbrales se adaptan de manera
    continua a las caracteristicas de la señal, ya que estan basados en la señal mas reciente y los picos de ruido que son detectados
    en las señales procesadas en curso.
    
    Si el programa no halla un complejo QRS en el intervalo de tiempo correspondiente al 166% del promedio del intevarlo RR actual
    el pico máximo detectado en ese intervalo de tiempo que permanece entre estos dos umbrales es considerado un posible complejo QRS
    y el mas bajo de los dos umbrales es aplicado. Para el caso de ritmos cardiacos irregulares, los dos umbrales son reducidos a la 
    mitad con el objetivo de incrementar la sensitividad de deteccion y evadir la perdida de latidos válidos.
    
    Una vez un complejo QRS valido es reconocido, hay un periodo refractario de 200ms antes de que el proximo pueda ser detectado
    ya que los complejos QRS no pueden ocurrir antes de esto desde el punto de vista fisiologico
    
    Este periodo refractario elimina la posibilidad de una faltsa deteccion como los es el multiple detonante en el mismo QRS durante
    este intervalo de tiempo
    
    El umbral bajo es usado si no se detecta QRS en un cierto intervalo de tiempo de tal manera que se usa una tecnica de busqueda en 
    retroceso para revisar algun complejo QRS
    
    para ser un pico, este debe exceder el valor UMBRAL I1 cuando la señal es inicialmente analizada o UMBRAL I2 si la busqueda en
    retroceso es requerida para hallar el QRS
    '''
    
    #paso 1: hallar la señal de ruido al substraer op del pasabandas del ECG original ( se aplicaron filtros pasa altos y pasa bajos de ganancia de 32 y 36)
    substraer_ruido = signal - (signal_pasa_bandas/(36*32))
    
    
#=======================================================================================================
    
    #paso 2: UMBRALES Y DETECION EN FORMA DE ONDA INTEGRADA: el conjunto de umbrales aplicados a la forma de onda de integracion es calculado
        #2.1: hallar los picos de la señal  filtrada y actualizando los umbrales
        
    signal_pico_i = 0 #inicializo pico de señal total inicial, REVISAR SI LO PUEDO QUITAR
    pico_ruido_i = 0 #inicializo pico de ruido total inicial REVISAR SI LO PUEDO QUITAR
        
    spk_i = signal_pico_i #Corriendo la señal pico estimada
    npk_i = pico_ruido_i #corriendo la señal ruido estimada
        
    
    index_pico_i = [0] #un arreglo para almacenar el indice de los picos
    umbral1_i = spk_i #valor inicial de umbral1 = 0
    detalles_pico_i = [] #arreglo para almacenar los detalles de voltaje y posicion de los picos
    
    for i in range(10, len(signal_integrada.iloc[:,1])):
        if signal_integrada.iloc[i,1]>signal_pico_i:
            signal_pico_i = signal_integrada.iloc[i,1] # actualizando maximo total de la señal
       
        if substraer_ruido.iloc[i,1] > pico_ruido_i:
            pico_ruido_i = substraer_ruido.iloc[i,1] #actuyalizo el maximo del ruido
        
        #actualizo umbrales
        npk_i = 0.125*pico_ruido_i + 0.875*npk_i
        spk_i = 0.125*signal_pico_i + 0.875*spk_i
        
        
        umbral1_i = npk_i + 0.25*(spk_i-npk_i)
        umbral2_i = 0.5*umbral1_i
        
        #ida
        if (signal_integrada.iloc[i,1] >= umbral1_i): #si la señal es mayor o igual al umbral
            if(index_pico_i[-1]+200 <i): # periodo refractario de 200ms antes del proximo
                index_pico_i.append(i) #almaceno indice del pico
                detalles_pico_i.append([signal_integrada.iloc[i,0], signal_integrada.iloc[i,1]]) #voltaje y posicion del pico
                
        #vuelta
        #2.2: chequeando la clasificacion como una señal pico usando umbral1 y luego umbral2
        else:
            if(signal_integrada.iloc[i,1]>=umbral2_i):
                if(index_pico_i[-1]+200 <i): # periodo refractario de 200ms antes del proximo
                    index_pico_i.append(i) #almaceno indice del pico
                    detalles_pico_i.append([signal_integrada.iloc[i,0], signal_integrada.iloc[i,1]]) #voltaje y posicion del pico
                    spk_i = 0.25*signal_pico_i + 0.75*spk_i #actualizo spk_i con el searchback
    
    print("STEP 2: thresholds and detection on integration waveform (s4)")
    print("Step 2: npk_i: ",npk_i)
    print("Step 2: spk_i: ",spk_i)
    print("Step 2: Threshold 1_i: ",umbral1_i)
    print("Step 2: Threshold 2_i: ",umbral2_i)
    print("Step 2: Number of peaks is: ",len(index_pico_i))
                
            
        
#=======================================================================================================
    #paso 3 : conjunto de umbrales aplicados al ECG filtrado determinado usando RR_miss limit
        #3.1: hallar los picos filtrados y actualizando los umbrales
        #3.2: chequear la clasificacion como un pico usando umbral1 y umbral2
    '''
    signal_pasa_bandas.iloc[:,1] = signal_pasa_bandas.iloc[:,1] + abs(min(signal_pasa_bandas.iloc[:,1])) # BPF op is in negative volts due to phase shift, makig it positive in order to find the peaks
    signal_pico_f = 0
    pico_ruido_f= 0
    spk_f = signal_pico_f
    npk_f = pico_ruido_f
    index_pico_f = [0]
    umbral1_f = spk_f
    detalles_pico_f = []
    
    for i in range(1, len(signal_integrada.iloc[:,1])):
        if signal_integrada.iloc[i,1] > signal_pico_f:
            signal_pico_f = signal_integrada.iloc[i,1] # actualizando maximo total de la señal
       
        if substraer_ruido.iloc[i,1] > pico_ruido_f:
            pico_ruido_f = substraer_ruido.iloc[i,1] #actuyalizo el maximo del ruido
        
        #actualizo umbrales
        spk_f = 0.125*signal_pico_f + 0.875*spk_f
        npk_f = 0.125*pico_ruido_f + 0.875*npk_f
        
        umbral1_f = npk_f + 0.25*(spk_f-npk_f)
        umbral2_f = 0.5*umbral1_f
        
        #ida
        if (signal_integrada.iloc[i,1] >= umbral1_f): #si la señal es mayor o igual al umbral
            if(index_pico_f[-1]+200 <i): # periodo refractario de 200ms antes del proximo
                index_pico_f.append(i) #almaceno indice del pico
                detalles_pico_f.append([signal_integrada.iloc[i,0], signal_integrada.iloc[i,1]]) #voltaje y posicion del pico
                
        #vuelta
        #2.2: chequeando la clasificacion como una señal pico usando umbral1 y luego umbral2
        else:
            if(signal_integrada.iloc[i,1]>=umbral2_f):
                if(index_pico_f[-1]+200 <i): # periodo refractario de 200ms antes del proximo
                    index_pico_f.append(i) #almaceno indice del pico
                    detalles_pico_f.append([signal_integrada.iloc[i,0], signal_integrada.iloc[i,1]]) #voltaje y posicion del pico
                    spk_f = 0.25*signal_pico_f + 0.75*spk_f #actualizo spk_i con el searchback
                   
    print("STEP 3: thresholds and detection on integration waveform (s4)")
    print("Step 3: npk_i: ",npk_f)
    print("Step 3: spk_i: ",spk_f)
    print("Step 3: Threshold 1_i: ",umbral1_f)
    print("Step 3: Threshold 2_i: ",umbral2_f)
    print("Step 3: Number of peaks is: ",len(index_pico_f))
    '''
    
    
    #paso 4 chequear los picos faltantes con rr limit
    
    umbral1_f = 0
    umbral2_f = 0
    signal_pico_f = signal_integrada.iloc[0,1]
    spk_f = 0
    npk_f = 0
    index_pico_f = index_pico_i
    avg = np.average(np.diff(index_pico_i)) #RR average
    detalles_pico_f = detalles_pico_i
    index_picos_fallidos = []
    
    
    for i in range(1, len(index_pico_f)):
        pico_f = signal_integrada.iloc[i-1, 1] #el maximo más reciente
        if(index_pico_f[i]-index_pico_f[i-1]>=1.66*avg): #buscando picos perdidos en la señal integrada
            for j in range(int(index_pico_f[i-1]+1, int(index_pico_f[i]))):
                if signal_pasa_bandas.iloc[j,1]>signal_pico_f:
                    signal_pico_f = signal_pasa_bandas.iloc[j,1] #actualizo el maximo
                
                if substraer_ruido.iloc[j,1] > pico_ruido_f:
                    pico_ruido_f = substraer_ruido.iloc[j,1] #actualizando el maximo del ruido
                
                #actualizar umbrales
                spk_f = 0.125*signal_pico_f + 0.875*spk_f
                npk_f = 0.125*pico_ruido_f + 0.875*npk_f
        
                umbral1_f = npk_f + 0.25*(spk_f-npk_f)
                umbral2_f = 0.5*umbral1_f
                
                #ida
                if (signal_integrada.iloc[j,1] >= umbral1_f): #si la señal es mayor o igual al umbral
                    if(index_pico_f[-1]+200 <j): # periodo refractario de 200ms antes del proximo
                        index_pico_f.append(j) #almaceno indice del pico
                        detalles_pico_f.append([signal_integrada.iloc[j,0], signal_integrada.iloc[j,1]]) #voltaje y posicion del pico
                
                #vuelta
                #2.2: chequeando la clasificacion como una señal pico usando umbral1 y luego umbral2
                else:
                    if(signal_integrada.iloc[j,1]>=umbral2_f):
                        if(index_pico_f[-1]+200 <j): # periodo refractario de 200ms antes del proximo
                            index_pico_f.append(j) #almaceno indice del pico
                            detalles_pico_f.append([signal_integrada.iloc[j,0], signal_integrada.iloc[j,1]]) #voltaje y posicion del pico
       
                
    print("STEP 5: Checking missed peaks using RR_MISSED LIMIT")             
    print("Step 5: npk_f: ",npk_f)
    print("Step 5: spk_f: ",spk_f)
    print("Step 5: Threshold 1_f: ",umbral1_f)
    print("Step 5: Threshold 2_f: ",umbral2_f)
    print("Step 5: Number of missed peaks is: ",len(index_picos_fallidos))
                
        

  
    
    
    
    
    
    
    
    
    #GRAFICAS
    fig = plt.figure(figsize=(20,20))
    fig.set_dpi(100)
    
    
    ax1=fig.add_subplot(6,1,1)
    ax1.plot(signal.iloc[:,0], signal.iloc[:,1])
    ax1.set_title('Original ECG')
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('mV')
    
    '''
    ax2=fig.add_subplot(6,1,2)
    ax2.plot(signal_pasa_bandas.iloc[:,0], signal_pasa_bandas.iloc[:,1])
    ax2.set_title('filtro pasabandas')
    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('mV')
    
    ax3=fig.add_subplot(6,1,3)
    ax3.plot(signal_derivativa.iloc[:,0],signal_derivativa.iloc[:,1])
    ax3.set_title('filtro derivativo')
    ax3.set_xlabel('Time (ms)')
    ax3.set_ylabel('mV')
    
    ax4=fig.add_subplot(6,1,4)
    ax4.plot(signal_cuadratica.iloc[:,0],signal_cuadratica.iloc[:,1])
    ax4.set_title('cuadratico')
    ax4.set_xlabel('Time (ms)')
    ax4.set_ylabel('mV')
    
    ax5=fig.add_subplot(6,1,5)
    ax5.plot(signal_integrada.iloc[:,0],signal_integrada.iloc[:,1])
    ax5.set_title('señal filtrada')
    ax5.set_xlabel('Time (ms)')
    ax5.set_ylabel('mV')
    
    plt.show()
    '''
    
    
    
    
    
    '''
    inicial_x = []
    inicial_y = []
    
    for i in range(1, len(detalles_pico_i)):
        inicial_x.append(detalles_pico_[i][0]
        inicial_y.append(detalles_pico_[i][1]
    final_x = []
    final_y = []
    
    for in range (1, len(detalles_pico_f)):
        final_x.append(detalles_pico_[i][0]
        final_y.append(detalles_pico_[i][1]
    #comparo 
    if len(inicial_x)<len(final_x):
        op1 = inicial_x
        op2 = final_x
        op3 = inicial_y
        op4 = final_y
    else:
        op1 = final_x
        op2 = inicial_x
        op3 = final_y
        op4 = inicial_y
    
    resultante_x = []
    resultante_y = []
    for i in range (1,len(op1)):
    '''
    