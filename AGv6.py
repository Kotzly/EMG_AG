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

bias=0              # If bias = 1, every cromossome will have a non frequency dependant DNA
maxGen=2000         # The max number of generations
startOver=True      # If True, the code will not consider the last simulation
tamPop=100          # Population number

maxFreq=180         # This is the max Frequency to consider #240
freqStep=3          # For freqStep=3 -> The code will consider [1,2,3],[3,4,5], etc# 3

taxaMut=0.01        # The mutation rate
taxaMutMin=0.01     # Minimum mutation rate
taxaMutMax=10.0     # Maximum mutation rate
chanceMut=4         # The chance of mutation (only for the "absolute" mutation)

bestTypes=[]        # Logging variable

continuous=False    # If True, the code will use a continuous fitness function (not recommended)
binaryFit=False     # If True, the fitness of each individual will be 1 for each right guess
                    # If False, it will be continuous if "continuous" is True, or 1 point if
                    # it guesses correctly, and 1.5 if it guesses with an confidence above 
                    # a "multFac" threshold
multFac=1.5         # 
binaryCrossChance=0.5           # The chance of ocurring a binary cross. 1 minus this 
                                # is the chance of ans mean crossing
vectorialMutationChance=0.5     # The chance of vectorial mutation. 1 minus this is 
                                # chance of an absolute mutation

taxaMutMult=4.0                 # The factor by which taxaMut will be multiplied

##############################################################################
guid=0                          # Individual ID (logging variable)

real=[]                         # DATA
origin=[]                       # DATA
fv=[]                           # DATA
frv=[]                          # DATA
nArq=0                          # DATA


# lastValues, botThs and topThs to be used in each archive
parameters={'bicepsinteiro.txt': [400,20,10],\
            'bicepsmetade.txt': [400,20,10],\
            'emgwk.txt': [400,20,10],\
            'emgmed.txt':[400,20,10],\
#            'xoxoxo.txt':[300,40,30],\
            'emgabrindo.txt':[500,20,20],\
            'emgapertando.txt':[400,20,20]}


# Method that return the number of right guesses of and individual
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

# This function just multiplies the chromossome of an individual by the frequency
# vector of an signal, return the result. The position that gets the higher
# number represent from which archive it thinks this signal belongs
def sayWho(indiv,real,fv):
    tam=len(fv)
    x= getFreqVector(fv)
    x=np.array(x)
    pont=x*indiv.cromo.freqFactor
    return pont
    
# Gets the *.txt files
def getArqs():
    arqVec=[]
    for arq in os.listdir('.'):
        if os.path.splitext(arq)[1]=='.txt':
            arqVec.append(arq)
    arqVec.reverse()
    return arqVec

# Chromossome class each chromossome mainly consists of an nArqs x (maxFreq/freqStep)
# matrix. Each column represent an archive, and each line represent a set of 
# freqStep frequencies
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


# Individual class
class ind:
    def __init__(self):
        global guid
        self.uid=guid
        guid+=1
        self.cromo = cromossome()
        self.fit=0
        self.marker='none'

# This function takes the fft data od an signal, and returns a similar vector,
# but instead of getting one element per frequency it take a number of freqStep
# frequencies, sum it and divide by freqStep
def getFreqVector(fv):
    x=[]
    tam=len(fv)
    for j in range(0,tam/2-5):
            k=int(round(float(j)*1000/float(tam)))
            
            if(k % 3 == 0):
                if len(x)==maxFreq/freqStep:
                    ##### BIAS ######
                    if bias==1:
                        x.append(-1)
                    #################
                    break
                x.append(sum(fv[k:k+freqStep])*2/tam)
    return x

# Read the data archives. The original signal is stored in origin. Each signal
# Is stored in real. real[arq][5] will contain the 5th signal of the arq'th file
# (as read by getArqs). The fft data will be stored at "fv" (indexes works the
# the same as for "real"). The frequency vector as got by getFrequencyVector
# is stored at frv
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

# Fitness method. Each signal frequency vector is multiplied by indiv
# chromossome. The numbers got are reconized as the score of each archive.
# Let's say that the 0th element gets the largest number. That mean this 
# individual "thinks" that that signal belongs to archive 4 (getArqs()[0])
# The fitnnes is then calculated by the number of right guesses of each
# individual
            
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
###############################################################################
                    if continuous:
                        score+=(np.max(pont[0])-np.min(pont[0]))/np.mean(pont[0]-np.min(pont[0]))
###############################################################################
                    else:
                        if np.max(np.array(pont)) >=multFac*np.mean(np.array(pont)):
                            score+=1.5
                        else:
                            score+=1
    ###########################################################################
                else:
                    score+=1
    return score
    
# Population class
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
            
# Mutation method. The mutation can be vetorial or absolute.
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
        
# Crossover by adding different chromossomes and dividing by the number of
# fathers
def meanCrossover(pais):
    filho= ind()
    
    somaFreqs = sum([pai.cromo.freqFactor for pai in pais])
    
    tam= len(pais)
    filho.cromo.freqFactor=somaFreqs/tam
    mutate(filho)
    filho.marker+=' meaned '
    
    return filho

# Crossover by replacing the sons genes by his mother's or his father's, with 
# 50% chance
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

# Mixed crossover
def weightedCrossover(pais):
    if random.random()<binaryCrossChance:
        return binaryCrossover(pais)
    else:
        return meanCrossover(pais)

# Tournament. Returns the best fitted individual
def torneio(pop):
    
    bestIndiv=pop.population[0]
    for indiv in pop.population:
        if indiv.fit>=bestIndiv.fit:
            bestIndiv=indiv
    return bestIndiv
        
# Generate a new population by performing crossovers with best and the reminder
# population
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

# Remove the n less fitted individuals, replacing them by new ones
def removeSuckers(pop,n):
    
    def getFit(indiv):
        return indiv.fit
    pop.population.sort(reverse=False,key=getFit)
    for i in range(0,n):
        pop.population[i]=ind()
        
# Returns the mean fitness of poppulation in pop
def getPopMean(pop):
    temp=0.0
    tam=len(pop.population)
    for indiv in pop.population:
        temp+=indiv.fit
    return temp/tam

# Not used. Divide all chromossomes of a population by the highest number 
# amongst them
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
    
# Plot a graph
def plotGens(best,mean):
    plt.plot(best,'go')
    plt.plot(mean,'b-')

# Class for controlling the GA variables
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
            if best.fit<=last.fit*1.001: #If the fitness doesnt grow by 0.1%
                self._counter+=1
            else:
    #            taxaMut=self._taxaMut
                chanceMut=self._chanceMut
                self._expansion=False
                self._counter=0
                ascendingCounter=0
                
            
            if self._counter==8:    # If the fitness doesnt grow in n generations
                if self._expansion: # If it the taxaMut is increasing 
                    if taxaMut<self._taxaMutMax:    # If taxaMut is less than the maximum
                        taxaMut*=self._taxaMutMult
                    else:           # If taxaMut bigger than the maximum
                        self._expansion=False
    
                else:               # If taxaMut is decreasing
                    if taxaMut>self._taxaMutMin:    # If it is bigger than the minimum
                        taxaMut/=self._taxaMutMult
                    else:                           # If it is less than the minimum
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
            
        print(gen,best.fit,':',best.marker,tamPop,taxaMut,chanceMut,maxGen)#,':', [p.fit for p in population]
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
    print([p.fit for p in pop.population])
    
    return pop
    
    
bestAll=main()