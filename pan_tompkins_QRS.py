# -*- coding: utf-8 -*-
"""
Created on Sat May  6 13:33:36 2023

@author: Ricardo j Quintero
"""

import numpy as np
import matplotlib.pyplot as plt

def filtro_pasa_bandas(entrada):

          salida = None
          signal = entrada.copy()
        
          #FILTRO PASA BAJOS
          #y(nT) = 2y(nT-T) - y(nT - 2T) + x(nT) - 2x(nT-6T) + x(nT - 12T)
          '''
          ¿por que los if?
          '''
          for i in range(len(entrada)):
              signal[i] = entrada[i]   
        
              if (i >=1): 
                signal[i] += 2*signal[i - 1]  
              if(i>=2):
                signal[i] -= signal[i - 2]   
              if(i >=6):
                signal[i] -= 2*entrada[i-6]  
              if(i>=12):
                signal[i] +=entrada[i-12] 
        
            
            
          #tomo el resultado y le copio la señal resultante
          salida = signal.copy()
        
          #FILTRO PASA ALTOS
          #y = 32x(nT-16T) - y(nT-T) - x(nT) + x(nT - 32T)
          for i in range(len(entrada)):
                 
                salida[i] = -1*signal[i]
                if(i>=1):
                 salida[i] -= salida[i-1]
                if(i>=16):
                  salida[i] += 32*signal[i-16]
                if(i >= 32):
                  salida[i] += signal[i-32]
          '''
          #se normaliza el valor REVISARRRR
          '''
          max_val = max(max(salida), -min(salida))
          salida = salida/max_val
        
          return salida

#=============================================================================
  
def derivada(entrada,fs):

        #inicializo el resultado, la entrada es la salida del anterior
        salida = entrada.copy()
        T = 1/fs
    
        #aplico el filtro de manera recursiva
        for i in range(len(entrada)):
    
          salida[i] = 0
    
          if(i>=1):
            salida[i] -= 2*entrada[i - 1]
    
          if (i>=2):
            salida[i] -= entrada[i-2]
          
          
          #---->averigua esta linea que viene <<<<------
          if(i>= 2 and i <= len(entrada)-2):
            salida[i] += 2*entrada[i+1]
    
          if(i>= 2 and i <= len(entrada)-3):
            salida[i] += entrada[i+2]
    
          salida[i] = (salida[i])/(8*T)
    
        return salida
  #====================================================  
  
def sqr(entrada):

        return np.square(entrada)
#====================================================================
def ventana_integradora(entrada,fs):

        salida = entrada.copy()
        win_size = round(0.15*fs)
        sum = 0
    
        for j in range(win_size):
          sum += entrada[j]/win_size
          salida[j] = sum
    
        for i in range(win_size, len(salida)):
          sum += entrada[i]/win_size
          sum -= entrada[i-win_size]/win_size
          salida[i] = sum
   
        return salida
    
#============================================================================================================
        
def pan_tompkins_qrs(signal,fs):
      entrada = signal.iloc[:,1].to_numpy()

      #filtro pasa_bandas
      global bpass
      bpass = filtro_pasa_bandas(entrada.copy())
      """
      #imprimo la gráfica del filtro pasabandas
      plt.figure(figsize = (20,4), dpi = 100)
      plt.xticks(np.arange(0, len(bpass)+1, 100))
      plt.plot(bpass[32:len(bpass)-2])
      plt.xlabel('Samples')
      plt.ylabel('MLIImV')
      plt.title("Bandpassed Signal")
      """
      
      
      
      #derivativo
      global derivativo
      derivativo = derivada(bpass.copy(),fs)
      """
      plt.figure(figsize = (20,4), dpi = 100)
      plt.xticks(np.arange(0, len(derivativo)+1, 150))
      plt.plot(derivativo[32:len(derivativo)-2])
      plt.xlabel('Samples')
      plt.ylabel('MLIImV')      
      plt.title("Derivative Signal")
      """
      
    
      global cuadratica
      cuadratica = sqr(derivativo.copy())
      """
      plt.figure(figsize = (20,4), dpi = 100)
      plt.xticks(np.arange(0, len(cuadratica)+1, 150))
      plt.plot(cuadratica[32:len(cuadratica)-2])
      plt.xlabel('Samples')
      plt.ylabel('MLIImV')
      plt.title("Squared Signal")
      """
      
      global ventana_movil
      ventana_movil = ventana_integradora(cuadratica.copy(),fs)
      """
      plt.figure(figsize = (20,4), dpi = 100)
      plt.xticks(np.arange(0, len(ventana_movil)+1, 150))
      plt.plot(ventana_movil[100:len(ventana_movil)-2])
      plt.xlabel('Samples')
      plt.ylabel('MLIImV')
      plt.title("Moving Window Integrated Signal")
      """
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    