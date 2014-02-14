from __future__ import division

## Self-Organizing Maps by Paras Chopra
## www.paraschopra.com
## paras1987@gmail.com
##
## Please give credit if you use my work.

## Enhancements done by Kyle Dickerson - Jan 2008
## kyle.dickerson@gmail.com
##
## Added division import from future, else all division is integer division
## Made batch update optional
## Cleaned up unused variables
## Fixed bug in assignment of self.radius in SOM class
## Various alterations for speed optimization
## Speed increase of ~ 10%

from random import *
from math import *
import sys

class Node:

    def __init__(self, FV_size=10, PV_size=10, Y=0, X=0):
        self.FV = [random() for i in range(FV_size)] # Feature Vector
        self.PV = [random() for i in range(PV_size)] # Prediction Vector
        self.X=X # X location
        self.Y=Y # Y location

class SOM:

    #Let radius=False if you want to autocalculate the radius
    def __init__(self, height=10, width=10, FV_size=10, PV_size=10, radius=False, learning_rate=0.005):
        self.height=height
        self.width=width
        self.radius = radius if radius else (height+width)/2
        self.total=height*width
        self.learning_rate=learning_rate
        self.nodes=[0]*(self.total)
        self.FV_size=FV_size
        self.PV_size=PV_size
        for i in range(self.height):
            for j in range(self.width):
                self.nodes[(i)*(self.width)+j]=Node(FV_size, PV_size,i,j)

    # Train_vector format: [ [FV[0], PV[0]],
    #                        [FV[1], PV[1]], so on..
    
    def train(self, iterations=1000, train_vector=[[[0.0],[0.0]]], batch_update=True):
        time_constant=iterations/log(self.radius)
        radius_decaying=0.0
        learning_rate_decaying=0.0
        influence=0.0
        for i in range(1,iterations+1):
            if batch_update:
                FV_delta = [[0.0]*self.FV_size]*self.total
                PV_delta = [[0.0]*self.PV_size]*self.total
            else:
                shuffle(train_vector) # Randomize the training data, so it gets presented in a different order each time
                
            radius_decaying=self.radius*exp(-1.0*i/time_constant)
            learning_rate_decaying=self.learning_rate*exp(-1.0*i/time_constant)            
            sys.stdout.write("\rTraining Iteration: " + str(i) + "/" + str(iterations))

            for  j in range(len(train_vector)):
                if not batch_update:
                    FV_delta = [[0.0]*self.FV_size]*self.total
                    PV_delta = [[0.0]*self.PV_size]*self.total
                
                input_FV=train_vector[j][0]
                input_PV=train_vector[j][1]
                best=self.best_match(input_FV)
                
                for k in range(self.total):
                    dist=self.distance(self.nodes[best],self.nodes[k])
                    if dist < radius_decaying:
                        influence=exp((-1.0*(dist**2))/(2*radius_decaying*i))
                    
                        inf_lrd = influence*learning_rate_decaying
                        FV_delta[k] = [FV_delta[k][l] + inf_lrd*(input_FV[l]-self.nodes[k].FV[l]) for l in range(self.FV_size)]
                        PV_delta[k] = [PV_delta[k][l] + inf_lrd*(input_PV[l]-self.nodes[k].PV[l]) for l in range(self.PV_size)]
                    
                if not batch_update:
                    for k in range(self.total):
                        self.nodes[k].FV = [self.nodes[k].FV[l] + FV_delta[k][l] for l in range(self.FV_size)]
                        self.nodes[k].PV = [self.nodes[k].PV[l] + PV_delta[k][l] for l in range(self.PV_size)]
            if batch_update:
                for k in range(self.total):
                    self.nodes[k].FV = [self.nodes[k].FV[l] + FV_delta[k][l] for l in range(self.FV_size)]
                    self.nodes[k].PV = [self.nodes[k].PV[l] + PV_delta[k][l] for l in range(self.PV_size)]
            
        sys.stdout.write("\n")


    #Returns prediction vector
    def predict(self, FV=[0.0]):
        best=self.best_match(FV)
        return self.nodes[best].PV
        
    #Returns best matching unit's index
    def best_match(self, target_FV=[0.0]):
        minimum=sqrt(self.FV_size) #Minimum distance
        minimum_index=1 #Minimum distance unit
        for i in range(self.total):
            temp=self.FV_distance(self.nodes[i].FV,target_FV)
            if temp<minimum:
                minimum=temp
                minimum_index=i
        return minimum_index

    def FV_distance(self, FV_1=[0.0], FV_2=[0.0]):
        temp=0.0
        for j in range(self.FV_size):
            temp+=(FV_1[j]-FV_2[j])**2
        temp=sqrt(temp)
        return temp

    def distance(self, node1, node2):
        return sqrt((node1.X-node2.X)**2+(node1.Y-node2.Y)**2)

if __name__ == "__main__":
    print "Initialization..."
    a=SOM(5,5,2,1,False,0.05)

    print "Training for the XOR function -- using Batch Update..."
    a.train(100,[[[1,0],[1]],[[1,1],[0]],[[0,1],[1]],[[0,0],[0]]])
    print "Predictions for the XOR function..."
    print "Prediction 0 0,", round(a.predict([0,0])[0])
    print "Prediction 1 0,", round(a.predict([1,0])[0])
    print "Prediction 0 1,", round(a.predict([0,1])[0])
    print "Prediction 1 1,", round(a.predict([1,1])[0])
    
    print "Training for the XOR function -- using Single-Step Update..."
    a.train(100,[[[1,0],[1]],[[1,1],[0]],[[0,1],[1]],[[0,0],[0]]], False)
    print "Predictions for the XOR function..."
    print "Prediction 0 0,", round(a.predict([0,0])[0])
    print "Prediction 1 0,", round(a.predict([1,0])[0])
    print "Prediction 0 1,", round(a.predict([0,1])[0])
    print "Prediction 1 1,", round(a.predict([1,1])[0])

