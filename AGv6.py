# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 09:40:28 2018

@author: Paulo Augusto
"""

import numpy as np
#from numpy import fft
import matplotlib.pyplot as plt
#import scipy.signal as sig
import os
import random 
import emgReaderClass as erc
import threading
import multiprocessing
import dataPlotter

bias=0
maxGen=2000
startOver=True
tamPop=100

maxFreq=180     #240
freqStep=3      #3

taxaMut=0.01
taxaMutMin=0.01
taxaMutMax=10.0
chanceMut=4

bestTypes=[]

continuous=False
binaryFit=False
multFac=1.5
binaryCrossChance=0.5
vectorialMutationChance=0.5

taxaMutMult=4.0

##############################################################################
guid=0
real=[]
origin=[]
fv=[]
frv=[]
nArq=0


parameters={'bicepsinteiro.txt': [400,20,10],\
            'bicepsmetade.txt': [400,20,10],\
            'emgwk.txt': [400,20,10],\
            'emgmed.txt':[400,20,10],\
#            'xoxoxo.txt':[300,40,30],\
            'emgabrindo.txt':[500,20,20],\
            'emgapertando.txt':[400,20,20]}


def countGuesses(indiv):
    
    arqVec=getArqs()
    nArq=len(arqVec)
    score=0
    
    for arq in range(0,nArq):
      
        for i in range(0,len(real[arq])):
            tam=len(real[arq][i])
            
            x= getFreqVector(fv[arq][i])
            x=np.array(x)
            
            pont=x*indiv.cromo.freqFactor
#            test.append(pont)
            
            if np.argmax(pont[0]) == arq:

                score+=1
    
    return score

def sayWho(indiv,real,fv):
    tam=len(fv)
    x= getFreqVector(fv)
    x=np.array(x)
    pont=x*indiv.cromo.freqFactor
#    pont+=tam*indiv.cromo.tamFactor
#    pont+=np.mean(real)*indiv.cromo.meanFactor
    return pont
    
def getArqs():
    arqVec=[]
    for arq in os.listdir('.'):
        if os.path.splitext(arq)[1]=='.txt':
            arqVec.append(arq)
    arqVec.reverse()
    return arqVec


class cromossome:
   
    def getRandomVec(self,n):
        v=[]
        for i in range(0,n):
           v.append(random.random()*2-1) 
        return v
    
    def __init__(self):
        self.freqFactor=[]
        n=len(getArqs())
        for i in range(0,maxFreq/freqStep+bias):    
            self.freqFactor.append(self.getRandomVec(n))
        self.freqFactor=np.matrix(self.freqFactor)

        
class ind:
    def __init__(self):
        global guid
        self.uid=guid
        guid+=1
        self.cromo = cromossome()
        self.fit=0
        self.marker='none'

def getFreqVector(fv):
    x=[]
    tam=len(fv)
    for j in range(0,tam/2-5):
            k=int(round(float(j)*1000/float(tam)))
            
            if(k % 3 == 0):
                if len(x)==maxFreq/freqStep:
### REMOVE essa linha
                    if bias==1:
                        x.append(-1)
###
                    break
                x.append(sum(fv[k:k+freqStep])*2/tam)
    return x


def readArqs():
    arqVec=getArqs()
    nArq=len(arqVec)
    reader=erc.emgReader()

    for arq in range(0,nArq):      
        origin.append([])
        real.append([])
        fv.append([])
        frv.append([])
        
        reader.lastValues=parameters[arqVec[arq]][0]
        reader.topThs=parameters[arqVec[arq]][1]
        reader.botThs=parameters[arqVec[arq]][2]
        
        origin[arq],real[arq],fv[arq] = reader.analyzeEmg(arqVec[arq],1000)
    for arq in range(0,nArq):      
        for i in range(0,len(fv[arq])):
            frv[arq].append(getFreqVector(fv[arq][i]))

def fitness(indiv):
    
    global nArq
    score=0
    for arq in range(0,nArq):
      
        for i in range(0,len(fv[arq])):
            tam=len(real[arq][i])
            
            pont=np.array(frv[arq][i])*indiv.cromo.freqFactor
#            print pont
            test=pont
            
            if np.argmax(pont) == arq:
                if not binaryFit:
    ###########################################################################
                    if continuous:
                        score+=(np.max(pont[0])-np.min(pont[0]))/np.mean(pont[0]-np.min(pont[0]))
    ###########################################################################
                    else:
                        if np.max(np.array(pont)) >=multFac*np.mean(np.array(pont)):
                            score+=1.5
                        else:
                            score+=1
    ###########################################################################
                else:
                    score+=1
    return score
    
class population:
    def __init__(self):
        self.population=[]

    def initPop(self,tamPop):
        for i in range(0,tamPop):
            self.population.append(ind())
    def evaluateAll(self):
        for ind in self.population:
            ind.fit=fitness(ind)        
    def getBest(self):
        return self.population[np.argmax(self.population)]
            

def mutate(indiv):

    global taxaMut,chanceMut

    if random.random()<vectorialMutationChance:
        vec=ind().cromo.freqFactor
        amp=np.sqrt(np.sum(pow(i,2) for i in vec.A1))
        vec/=amp
        vec*=taxaMut*random.random()
        indiv.cromo.freqFactor+=vec
        indiv.marker='vectorial'
#    for line in indiv.cromo.freqFactor:
#        for i in range(0,len(np.array(line)[0])):
#            if random.random()*1000<chanceMut:
#                line[0,i]+=mut*random.random()
    else:
        for line in indiv.cromo.freqFactor:
            for i in range(0,len(np.array(line)[0])):
                if random.random()*1000<chanceMut:
                    if random.random()<0.5:
                        mut =  taxaMut
                    else:
                        mut = -taxaMut
                    line[0,i]+=mut*random.random()
    
        indiv.marker='absolute'
def meanCrossover(pais):
    filho= ind()
    
    somaFreqs = sum([pai.cromo.freqFactor for pai in pais])
    
    tam= len(pais)
    filho.cromo.freqFactor=somaFreqs/tam
    mutate(filho)
    filho.marker+=' meaned '
    
    return filho

def binaryCrossover(pais):
    filho=ind()
    for i in range(0,len(filho.cromo.freqFactor)):
        for j in range(0,len(filho.cromo.freqFactor[0].A1)):
            if random.random()<0.5:
                filho.cromo.freqFactor[i,j]=pais[0].cromo.freqFactor[i,j]
            else:
                filho.cromo.freqFactor[i,j]=pais[1].cromo.freqFactor[i,j] 
    mutate(filho)
    filho.marker+=' binerized '                  
    return filho          

def weightedCrossover(pais):
    if random.random()<binaryCrossChance:
        return binaryCrossover(pais)
    else:
        return meanCrossover(pais)
def torneio(pop):
    
    bestIndiv=pop.population[0]
    for indiv in pop.population:
        if indiv.fit>=bestIndiv.fit:
            bestIndiv=indiv
    return bestIndiv
        
def genNewPop(best,pop):   
    newpop=population()
    for indiv in pop.population:
        if indiv == best:
            newpop.population.append(indiv)
            continue
        else:
            temp=weightedCrossover([best,indiv])
            newpop.population.append(temp)
    return newpop

def removeSuckers(pop,n):
    
    def getFit(indiv):
        return indiv.fit
    pop.population.sort(reverse=False,key=getFit)
    for i in range(0,n):
        pop.population[i]=ind()
        
        
def getPopMean(pop):
    temp=0.0
    tam=len(pop.population)
    for indiv in pop.population:
        temp+=indiv.fit
    return temp/tam

def normalizePop(pop):
    for indiv in pop.population:    
        maxF=0    
        for line in indiv.cromo.freqFactor:
            for i in range(0,len(np.array(line)[0])):        
                if abs(line[0,i]) > maxF:
                    maxF=abs(line[0,i])
                
        for line in indiv.cromo.freqFactor:
            for i in range(0,len(np.array(line)[0])):        
                line[0,i]/=maxF
    
def plotGens(best,mean):
    plt.plot(best,'go')
    plt.plot(mean,'b-')

class populationControl():
    global  tamPop,\
            taxaMut,\
            chanceMut,\
            bestAll,\
            bias,\
            maxGen,\
            tamPop,\
            taxaMut,\
            taxaMutMax,\
            chanceMut,\
            continuous,\
            binaryFit,\
            multFac,\
            binaryCrossChance,\
            taxaMutMult,\
            taxaMutMin

    def __init__(self):
        self._tamPop=tamPop
        self._taxaMut=taxaMut
        self._chanceMut=chanceMut
        self._bias=bias
        self._maxGen=maxGen
        self._tamPop=tamPop            
        self._taxaMutMin=taxaMutMin
        self._taxaMutMax=taxaMutMax
        self._chanceMut=chanceMut
        self._continuous=continuous
        self._binaryFit=binaryFit
        self._multFac=multFac
        self._binaryCrossChance=binaryCrossChance
        self._taxaMutMult=taxaMutMult
        self._counter=0
        self._expansion=False
    def control(self,gen,counter,best,last):
        global taxaMut
#        taxaMut=self._taxaMutMax
        ascendingCounter=0
        if gen>25:
            if best.fit<=last.fit*1.001:
                self._counter+=1
            else:
    #            taxaMut=self._taxaMut
                chanceMut=self._chanceMut
                self._expansion=False
                self._counter=0
                ascendingCounter=0
                
            
            if self._counter==8:
                if self._expansion:
                    if taxaMut<self._taxaMutMax:
                        taxaMut*=self._taxaMutMult
                    else:
    #                    if ascendingCounter>3:
                        self._expansion=False
    #                        ascendingCounter+=1
    
                else:
                    if taxaMut>self._taxaMutMin:
                        taxaMut/=self._taxaMutMult
                    else:
                        self._expansion=True    
                
                self._counter=0  
                

def main():
    
    global  maxFreq,\
            freqStep,\
            tamPop,\
            taxaMut,\
            chanceMut,\
            nArq,\
            bestAll,\
            startOver,\
            bestTypes
            
    nArq=len(getArqs())
    gen=0
    counter=0
    last=ind()
    bestVec=[]
    meanVec=[]
    taxaVec=[]
    taxaMut=taxaMutMax
    if startOver:
        pop = population()
        pop.initPop(tamPop)
    else:
        pop=bestAll


#    plotter=dataPlotter.dataPlotter('Geracao','Melhor de Todos',bestVec)
#    threading.Thread(target=plotter.start).start()
    controller=populationControl()
    readArqs()
    while gen<maxGen:
        gen+=1
        pop.evaluateAll()
        best=torneio(pop)
        
        if not last.uid==best.uid:
            bestTypes.append(best.marker)
            
        print gen,best.fit,':',best.marker,tamPop,taxaMut,chanceMut,maxGen#,':', [p.fit for p in population]
        pop=genNewPop(best,pop)
    ###########################################################################
        controller.control(gen,counter,best,last)
        last=best
        
        taxaVec.append(20*np.log(taxaMut))
        bestVec.append(last.fit)
        meanVec.append(getPopMean(pop))
    ###########################################################################
#        createSuckers(pop.tamPop/3)
        removeSuckers(pop,tamPop/5)
#        normalizePop(pop)
    
    plotGens(bestVec,meanVec)
    plotGens(bestVec,taxaVec)
    pop.evaluateAll()
    print [p.fit for p in pop.population]
    
    return pop
    
    
bestAll=main()