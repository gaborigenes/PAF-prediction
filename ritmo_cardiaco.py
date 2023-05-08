# -*- coding: utf-8 -*-
"""
Created on Sat May  6 15:51:40 2023

@author: Ricardo j Quintero
"""

from scipy import signal as sg
import numpy as np
import pan_tompkins_QRS

class ritmo_cardiaco():

  def __init__(self,signal,freq):

    # inicializan las variables
    #SPKI estimación de la señal pico corriendo (es un pico que el algoritmo ha estblecido anteriormente como QRS)
    #NPKI #estimación del pico de ruido (con QRS establecido)
    #SPKF #estimación de la señal pico corriendo (cuando el complejo se consigue usando el segundo umbral)
    #NPKF estimación del pico de ruido (usando el segundo umbral)
    #RR1 primer promedio RR
    #THRESHOLD_I1 primer umbral inicialmente aplicado a la forma de ondas integrada
    #THRESHOLD_F1 primer umbral aplicado a la forma de ondas filtrada
    #RR2 segundo promedio RR
    #THRESHOLD_I2 segundo umbral inicialmente aplicado a la forma de ondas integrada
    #THRESHOLD_F2 segundo umbral aplicado a la forma de ondas filtrada

    self.RR1, self.RR2, self.probable_peaks, self.r_locs, self.peaks, self.result = ([] for i in range(6))
    self.SPKI, self.NPKI, self.THRESHOLD_I1, self.THRESHOLD_I2, self.SPKF, self.NPKF, self.THRESHOLD_F1, self.THRESHOLD_F2 = (0 for i in range(8))
   
    #bandera de onda T
    self.T_wave = False
    self.m_win = pan_tompkins_QRS.ventana_movil 
    self.b_pass = pan_tompkins_QRS.bpass
    self.samp_freq = freq
    self.signal = signal
    self.win_150ms = round(0.15*self.samp_freq)

    self.RR_low_limit = 0
    self.RR_high_limit = 0
    self.RR_missed_limit = 0
    self.RR_average1= 0


#===============================================================================
  def pico_aprox(self):

    #FFT convolucion para filtrar la señal y suavizarla (entrada, ventana, tamño de salida)
    slopes = sg.fftconvolve(self.m_win, np.full((25,), 1) / 25, mode='same')

    # Finding approximate peak locations
    #no entiendo bien esta linea, despues de hallar el punto de origen y fin de la onda, se calcula el ancho y si es menor o igual a 64 ms (la mitad de numero de muestras)
    for i in range(round(0.5*self.samp_freq) + 1,len(slopes)-1):      
        #cambio de signo
        if (slopes[i] > slopes[i-1]) and (slopes[i+1] < slopes[i]):
            self.peaks.append(i) 
  
#===============================================================================

  def ajustar_intervalo_rr(self,ind):

    #entrada: indice actual en el arreglo de picos

    #hallando los 8 RR más recientes desde el indice actual hasta el 7 anteriores
    #¿por que se divide entre la frecuencia?
    self.RR1 = np.diff(self.peaks[max(0,ind-8): ind+1])/self.samp_freq
    

    #calculo de promedio de RR
    self.RR_average1 = np.mean(self.RR1)
    RR_average2 = self.RR_average1


    #Hallando los 8 RR mas recientes entre RR low limit y RR high limit
    if(ind >= 8):
      for i in range(0,8):
        if(self.RR_low_limit < self.RR1[i] < self.RR_high_limit):
           self.RR2.append(self.RR1[i])

           #movimiento de la ventana
           if(len(self.RR2) > 8):
            self.RR2.remove(self.RR2[0])
            RR_average2 = np.mean(self.RR2)

    '''
    ajustando RR low limit y RR high limit
    los limites RR son
    RR LOW LIMIT = 92% RR AVERAGE2
    RR HIGH LIMIT = 116% RR AVERAGE 2
    RR MISSED LIMIT = 166% RR AVERAGE 2
    '''  
    if(len(self.RR2)>7 or ind < 8):
      self.RR_low_limit = 0.92*RR_average2
      self.RR_high_limit = 1.16*RR_average2
      self.RR_missed_limit = 1.66*RR_average2


#============================================================================================================  

  def searchback(self, peak_val, RRn, sb_win):
    '''
    Searchback
    parametros:
    peak_Val  localización del pico a considerar
    RRn: intervalo más reciente
    sb_win: searchback window = RRn*frecuencia PREGUNTARRRR

    '''

    #chequear si el intervalo RR más reciente es mayor que RR_missed_limit
    if(RRn > self.RR_missed_limit):
      #aqui se inicia la ventana de retroceso
      win_rr = self.m_win[peak_val - sb_win +1 : peak_val + 1]

      #hallar las localizaciones de x dentro de la ventana teniendo valores y mayores al umbral I1
      coord = np.asarray(win_rr > self.THRESHOLD_I1).nonzero()[0]


      #hallar la localización x del maximo pico en la ventana de busqueda
      if(len(coord) > 0):
        for posicion in coord:
          if(win_rr[posicion] == max(win_rr[coord])):
            x_max = posicion
            break

      else:
        x_max = None

      #Si se ha conseguido un valor pico
      if (x_max is not None):   
            # Update the thresholds corresponding to moving window integration
            self.SPKI = 0.25 * self.m_win[x_max] + 0.75 * self.SPKI                         
            self.Threshold_I1 = self.NPKI + 0.25 * (self.SPKI - self.NPKI)
            self.Threshold_I2 = 0.5 * self.Threshold_I1         

            # Initialize a window to searchback 
            win_rr = self.b_pass[x_max - self.win_150ms: min(len(self.b_pass) -1, x_max)]  

            # Find the x locations inside the window having y values greater than Threshold F1                   
            coord = np.asarray(win_rr > self.Threshold_F1).nonzero()[0]

            # Find the x location of the max peak value in the search window
            if (len(coord) > 0):
              for posicion in coord:
                  if (win_rr[posicion] == max(win_rr[coord])):
                      r_max = posicion
                      break
            else:
              r_max = None

            # If the max peak value is found
            if (r_max is not None):
              # Update the thresholds corresponding to bandpass filter
              if self.b_pass[r_max] > self.Threshold_F2:                                                        
                  self.SPKF = 0.25 * self.b_pass[r_max] + 0.75 * self.SPKF                            
                  self.Threshold_F1 = self.NPKF + 0.25 * (self.SPKF - self.NPKF)
                  self.Threshold_F2 = 0.5 * self.Threshold_F1      

                  # Append the probable R peak location                      
                  self.r_locs.append(r_max)     


#============================================================================================================  

  def hallar_onda_t(self, peak_val, RRn, ind, prev_ind):
    '''
    parametros:
    peak_val: el valor maximo
    RRn: el intervalo RR más reciente
    ind: indice actual
    prev_ind: indice anterior 

    '''
    #cuando un RR está entre 200ms<x<360ms se realiza un juego donde se ha hallado correctamnte el QRS es una onda T


    #si la máxima pendiente ocurre durante esa onda es menos a la mitad del QRS anterior, es una onda T, de lo contrario es un QRS
  
    if(self.m_win[peak_val] >= self.THRESHOLD_I1):
      if(ind >0 and 0.20 <RRn< 0.36):
        #se halla las mendiente de la onda actual y la anterior
        curr_slope = max(np.diff(self.m_win[peak_val-round(self.win_150ms/2): peak_val + 1]))
        last_slope = max(np.diff(self.m_win[self.peaks[prev_ind] - round(self.win_150ms/2) : self.peaks[prev_ind] + 1]))


        #si la pendiente de la onda actual es menor que la mitad que la pendiente de la onda anterior
        if(curr_slope< 0.5*last_slope):
          #la onda te es hallada y se actualiza el umbral de ruido  
          self.T_wave = True
          self.NPKI = 0.125*self.m_win[peak_val] + 0.875*self.NPKI

          #si no se halla una una onda T actualiza los umbrales
      if(not self.T_wave):
          if (self.probable_peaks[ind] > self.THRESHOLD_F1):   
                self.SPKI = 0.125 * self.m_win[peak_val]  + 0.875 * self.SPKI                                         
                self.SPKF = 0.125 * self.b_pass[ind] + 0.875 * self.SPKF 
              
                #agrega la localización del pico r probable
                self.r_locs.append(self.probable_peaks[ind])
          
          else:

                self.SPKI = 0.125 * self.m_win[peak_val]  + 0.875 * self.SPKI
                self.NPKF = 0.125 * self.b_pass[ind] + 0.875 * self.NPKF 
    
    #actualiza los umbrales de ruido
    elif(self.m_win[peak_val] < self.THRESHOLD_I1) or (self.THRESHOLD_I1 < self.m_win[peak_val] < self.THRESHOLD_I2):
      self.NPKI = 0.125 * self.m_win[peak_val]  + 0.875 * self.NPKI  
      self.NPKF = 0.125 * self.b_pass[ind] + 0.875 * self.NPKF



#============================================================================================================  

  def ajustar_umbrales(self, peak_val, ind):
    '''
    ajustar el ruido y la señal umbral durante la fase de aprendizaje
    '''
    if (self.m_win[peak_val] >= self.THRESHOLD_I1): 
        # Update signal threshold
        self.SPKI = 0.125 * self.m_win[peak_val]  + 0.875 * self.SPKI

        if (self.probable_peaks[ind] > self.THRESHOLD_F1):                                            
            self.SPKF = 0.125 * self.b_pass[ind] + 0.875 * self.SPKF 

            # Append the probable R peak location
            self.r_locs.append(self.probable_peaks[ind])  

        else:
            # Update noise threshold
            self.NPKF = 0.125 * self.b_pass[ind] + 0.875 * self.NPKF                                    
        
    # Update noise thresholds    
    elif (self.m_win[peak_val] < self.THRESHOLD_I2) or (self.THRESHOLD_I2 < self.m_win[peak_val] < self.THRESHOLD_I1):
        self.NPKI = 0.125 * self.m_win[peak_val]  + 0.875 * self.NPKI  
        self.NPKF = 0.125 * self.b_pass[ind] + 0.875 * self.NPKF

#============================================================================================================  

  def actualizar_umbrales(self):
    '''
    actualizar el umbral de ruido y señal para la proxima iteración
    
    '''
    self.THRESHOLD_I1 = self.NPKI + 0.25 * (self.SPKI - self.NPKI)
    self.THRESHOLD_F1 = self.NPKF + 0.25 * (self.SPKF - self.NPKF)
    self.THRESHOLD_I2 = 0.5 * self.THRESHOLD_I1 
    self.THRESHOLD_F2 = 0.5 * self.THRESHOLD_F1
    self.T_wave = False 



#============================================================================================================  

  def ecg_searchback(self):
    '''
    retroceso en el ECG para incrementar la eficiencia

    '''

    #filtra las localizaciones de los picos R
    self.r_locs = np.unique(np.array(self.r_locs).astype(int))

    #inicializa una ventana para searchback

    win_200ms = round(0.2*self.samp_freq)

    for r_val in self.r_locs:
      coord = np.arange(r_val - win_200ms, min(len(self.signal), r_val + win_200ms + 1), 1)

      # Find the x location of the max peak value
      if (len(coord) > 0):
        for pos in coord:
            if (self.signal[pos] == max(self.signal[coord])):
                x_max = pos
                break
      else:
          x_max = None

      # Append the peak location
      if (x_max is not None):   
         self.result.append(x_max)


#============================================================================================================  

  def hallar_picos_r(self):
    '''
    R Peak Detection
    '''

    # Find approximate peak locations
    self.pico_aprox()

    # Iterate over possible peak locations
    for ind in range(len(self.peaks)):

        # Initialize the search window for peak detection
        peak_val = self.peaks[ind]
        win_300ms = np.arange(max(0, self.peaks[ind] - self.win_150ms), min(self.peaks[ind] + self.win_150ms, len(self.b_pass)-1), 1)
        max_val = max(self.b_pass[win_300ms], default = 0)

        # Find the x location of the max peak value
        if (max_val != 0):        
          x_coord = np.asarray(self.b_pass == max_val).nonzero()
          self.probable_peaks.append(x_coord[0][0])
        
        if (ind < len(self.probable_peaks) and ind != 0):
            # Adjust RR interval and limits
            self.ajustar_intervalo_rr(ind)
            
            # Adjust thresholds in case of irregular beats
            if (self.RR_average1 < self.RR_low_limit or self.RR_average1 > self.RR_missed_limit): 
                self.THRESHOLD_I1 /= 2
                self.THRESHOLD_F1 /= 2

            RRn = self.RR1[-1]

            # Searchback
            self.searchback(peak_val,RRn,round(RRn*self.samp_freq))

            # T Wave Identification
            self.hallar_onda_t(peak_val,RRn,ind,ind-1)

        else:
          # Adjust threholds
          self.ajustar_umbrales(peak_val,ind)

        # Update threholds for next iteration
        self.actualizar_umbrales()

    # Searchback in ECG signal 
    self.ecg_searchback()

    return self.result