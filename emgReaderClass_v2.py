# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 17:14:36 2018

@author: Paulo Augusto
"""
import numpy as np
from numpy import fft
import matplotlib.pyplot as plt
import csv
import scipy.signal as sig
import os

# Class for reading EMG signals 
class emgReader:
    
    lastValues=400
    topThs=20
    botThs=10
    fs=2000
    
    # Gets the 1,2,3...n'th more predominant frequencies as an array
    def nBest(self,vector,n):

        i=[]
        maxi=0.0
        temp=[x for x in vector]
        
        for j in range(0,n):
            maxi=np.nanargmax(temp)
            temp[maxi]=np.nan
            i.append(maxi*self.fs/len(vector))
    
        return i
        
    # Get the n more predominant frequencies for all archives
    def getFreqs(self,fftv,n):
    #    for fftsig in emgFft:
    #        for i in range(0,len(fftsig)/2):
    #            if fftsig[i]==max(fftsig[0:len(fftsig)/2]):
    #                freqs.append(float(i)*1000/len(fftsig))
        freqs=[]
        i=0
        for fftsignal in fftv:
            freqs.append([])
            freqs[i]=self.nBest(fftsignal[0:len(fftsignal)/2],n)
            i+=1
        return freqs
    
    # This function separates different EMG signals in the same archive. The 
    # variables lastValues,botThs and topThs can (and must) be changed in order
    # to detect different signals with more precision
    # To detect different EMG signals in one archive, the code calculates the 
    # mean value of the last lastValues elements. If that goes above topThs
    # threshold, the beggining of a new signal is detected, and it is ended
    # when that mean goes below botThs.
    def separateEmg(self,vector,origin):  
        
        realEmg=[]
        numberOfEmg=0
        temp=0
        flag=0
        i=0
        for i in range(0,len(vector)):
            temp+= vector[i]
            if i>=self.lastValues:
                temp-=vector[i-self.lastValues]
            media=temp/self.lastValues
#            if i%100==0:
#                print media
            if media>self.topThs:        
                if flag==0:
#                    print 'Signal: ',i
                    realEmg.append([])
                    for j in range(i-self.lastValues,i):
                        realEmg[numberOfEmg].append(origin[j])  
                flag=1
        #        realEmg[numberOfEmg].append(emgAbs[i])
                  
                realEmg[numberOfEmg].append(origin[i])
            if flag==1 and media<self.botThs:
                numberOfEmg+=1
                flag=0
            
        return realEmg
    
    # Apply a 4 order butterworth filter to the signal
    def filterData(self,vector,fs):
        high=240.0
        low=3.0
        b,a=sig.butter(4,[low/(fs/2) , high/(fs/2)],btype='band')
        zi = sig.lfilter_zi(b, a)
        z, _ = sig.lfilter(b, a, vector, zi=zi*vector[0])
        vector=z
    
    # Unused 
    def mean(self,vector,first,last):
        temp=0
        for i in range(0,last):
            temp+=int(vector[first-i])
        return temp/last
              
    # Read the EMG signal Files. This signals are acquired using the bioPLUX 
    # software
    
    def getCsvData(self,muscle,folder):
        global fs
        def getArqs():
            arqVec=[]
            for arq in os.listdir('./'+folder+'/'):
                if os.path.splitext(arq)[1]=='.csv':
                    arqVec.append(arq)
            return arqVec
 
        real=[]
        fftvec=[]
        arqVec=getArqs()


        for arq in arqVec:
            if arq[0:3]==muscle:
                while len(real)<int(arq[4:6]):
                    fftvec.append([])
                    real.append([])
                with open('./'+folder+'/'+arq) as csvfile:
                    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                    data = [float(c) for c in spamreader.next()[0].split(',')]
                    self.filterData(data,2000)
                    real[int(arq[4:6])-1].append(data)
                    fftvec[int(arq[4:6])-1].append([abs(c) for c in fft.fft(data)])

        
        return real,fftvec

    
    def getData(self,arq):
        data=[]
        


        with open(arq) as txtfile:
            line=[['1'],]
            while True:
                line=txtfile.readline();
                if line=='':
                    break
                if line[0][0]=='#':
                    continue
                line=line.split('\t')
                line[3]=line[3][0:4]
                data.append(line)

        return data
                  
    # This is the core function of this class. The DC signal is removed, then 
    # then the signal is rectified, and the EMG signal are separated. A fft 
    # is performed.
    def analyzeEmg(self,arq,fs):
        


        data = self.getData(arq)
            
        emgValues=[float(line[3]) for line in data];     
        
        self.filterData(emgValues,fs)
        
        emgValuesMean= np.mean(emgValues)
        emgAbs = [abs(x-emgValuesMean) for x in emgValues]
        
        
        realEmg= self.separateEmg(emgAbs,emgValues)    
        
        realEmgAbs=[signal-np.mean(signal) for signal in realEmg]
        emgFft=[abs(fft.fft(signal)) for signal in realEmgAbs]
        
        freqs= self.getFreqs(emgFft,5)
        return emgValues,realEmg,emgFft
        
    def plotAllEmg():
        for signal in realEmg:
            plt.plot(signal)
    
    def plotAllFft():
        for signal in emgFft:
            plot.plot(signal)

#z,p,k=sig.butter(4,[3.0/fs, 200.0/fs],btype='bandpass',analog=False)
#x=[1,2,3,0,-9,2,3,1,3,19,1,12,3,-1,12,4,12,34]
#b,a=sig.butter(4,[3.0/(fs/2) , 200.0/(fs/2)],btype='band')
#zi = sig.lfilter_zi(b, a)
#z, _ = sig.lfilter(b, a, x, zi=zi*x[0])
#emgValues=z


#t = np.linspace(-1, 1, 201)
#x = (np.sin(2*np.pi*0.75*t*(1-t) + 2.1) + 0.1*np.sin(2*np.pi*1.25*t + 1) + 0.18*np.cos(2*np.pi*3.85*t))
#xn = x + np.random.randn(len(t)) * 0.08
#b, a = sig.butter(3, 0.05)
#zi = sig.lfilter_zi(b, a)
#z, _ = sig.lfilter(b, a, xn, zi=zi*xn[0])
#plt.figure
#plt.plot(t, xn, 'b', alpha=0.75)
#plt.plot(t, z, 'r--')
#plt.grid(True)
#plt.show()