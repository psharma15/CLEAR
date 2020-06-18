# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:46:52 2020

@author: Pragya Sharma

Code for AI Imaging - v1, using pre-calculated 8 feature vectors
"""

import os
import numpy as np
import matplotlib
import torch
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import csv
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, random_split
from torchvision import transforms, utils
import time
import pandas as pd
import scipy.io
import sklearn.metrics
import seaborn as sns
import random

random.seed(1)
torch.manual_seed(1)
torch.cuda.manual_seed(1)
np.random.seed(1)

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# Writing custom dataset
class AI_Image_dataset(Dataset):
    def __init__(self, file_name):
        data = scipy.io.loadmat(file_name)
        self.X = data['featVec']
        # Only use with object features
        self.X = self.X[:,[0,4],:]
        self.Y = data['labelVec']

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return (self.Y[idx],self.X[idx])

#Network definition
class Net(nn.Module):
    def __init__(self):
        super(Net,self).__init__()
        n_CW = 120 # This directly affects accuracy if value is low.
        self.conv1 = nn.Conv1d(2,n_CW,kernel_size=100, stride=20, padding=49)
        self.conv2 = nn.Conv1d(n_CW,int(n_CW/2) , kernel_size=50, stride=10, padding=31)
        self.conv3 = nn.Conv1d(int(n_CW/2),6,kernel_size=16, stride=10, padding=0)
        self.maxPool1 = nn.MaxPool1d(kernel_size=5,stride=1,padding=2)
        self.avgPool1 = nn.AvgPool1d(kernel_size=5)       
        self.drop1 = nn.Dropout(p=0.2)
        
    def forward(self, x):
        x = (F.relu(self.conv1(x)))
        print(x.shape)
        #x = self.maxPool1(F.relu(self.conv2(x)))
        x = (F.relu(self.conv2(x)))
        print(x.shape)
        x = self.avgPool1(F.relu(self.conv3(x)))
        print(x.shape)
        x = self.drop1(x)
        print(x.shape)
        #x = x.view(x.shape[0],-1)              
        return x

#Network definition

if __name__ == '__main__':
    
    file_name=r"E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\3-D imaging results Guoyi\TrueModel2019\ResultsFeaturesDict\feat_label_setup1_v2.mat"

    dataset = AI_Image_dataset(file_name)
    print(len(dataset))
    #print(dataset[100])
    #print(dataset[122:124])

    totaldata_loader = DataLoader(dataset, batch_size=len(dataset), shuffle=True)
    #print(next(iter(dataloader)))
    batchsize_train = 50
    fractionTrain = 0.75
    trainset, testset = random_split(dataset, [int(fractionTrain*len(dataset)), len(dataset)-int(fractionTrain*len(dataset))])
    train_loader = DataLoader(trainset, batch_size=batchsize_train,shuffle=True)
    batchsize_test = len(testset)
    test_loader = DataLoader(testset, batch_size=batchsize_test, shuffle=True)
            
    
    #hyperparameter definition    
    learning_rate = 0.01
    momentum = 0.1
    random_seed=1
    torch.backends.cudnn.enabled = False
    torch.manual_seed(random_seed)
    
    train_loss_epoch = []
    test_acc_epoch =[]
    total_acc_epoch = []
    #epoch_range = np.arange(25,225,25)
    epoch_range = np.array([100])
    for EPOCHS in epoch_range:
    
        network = Net()
        
        #training
        optimizer = optim.SGD(network.parameters(), lr=learning_rate, momentum=momentum)
        
        
        plt.figure()
        
        #network.train()
        Training_Loss = []
        start_time = time.time()
        for epoch in range(EPOCHS):
            train_loss = 0
            for Y,X in train_loader:
                 #X = torch.rand((3, 2,1024))
                X = X.view(-1,X.shape[1],X.shape[2])
               # X = X.view(-1,1,X.shape[1])
                Y = Y.view(-1,Y.shape[1])
                X = X.float()
                Y = Y.long()            
                current_batchsize = X.shape[0]
                optimizer.zero_grad()
                output = network(X)
                loss = F.cross_entropy(output,Y)
                train_loss = train_loss + (loss*current_batchsize/len(train_loader.dataset))
                loss.backward()                     #calculate the gradient decent
                optimizer.step()                    #update the weight
            Training_Loss.append(train_loss)
            if epoch%5==0:
                plt.clf()
                symbol = output.detach()
                symbol =  symbol.max(dim=1).indices
                plt.plot(Y[0], 'bo')
                plt.plot(symbol[0], 'rx')
                plt.title('epoch %i: training loss is ' %epoch + '%f'%loss)
                plt.show()
                plt.pause(0.5)
        
        train_loss_epoch.append(Training_Loss[-1].data.numpy())
        training_time=time.time()-start_time
        print('total training time is',training_time)
        
        plt.figure()
        plt.plot(Training_Loss)
        plt.title("training loss")
        plt.xlabel("epoch")
        plt.show()
        
        #evaluation
        network.eval()
        correct = []
        test_loss = 0
        
        for Y, X in test_loader:
           X = X.view(-1,X.shape[1],X.shape[2])
           # X = X.view(-1,1,X.shape[1])
           Y = Y.view(-1,Y.shape[1])
           X = X.float()
           Y = Y.long()      
           current_batchsize = X.shape[0]
           output = network(X)
           test_loss = test_loss + (F.cross_entropy(output,Y) * current_batchsize/len(test_loader.dataset))
           symbol = output.detach()
           symbol = symbol.max(dim=1).indices
           correct.append(Y.eq(symbol).numpy())
        correct = [ item for sublist in correct for item in sublist]   
        correct = np.asarray(correct)  
        correct = correct.reshape(-1)
        accuracy = correct.sum()/len(test_loader.dataset)
            
        print('test_lost is',(test_loss/(len(test_loader.dataset)/batchsize_test)))
        print('accuracy is: ',accuracy)
        test_acc_epoch.append(accuracy)

    plt.figure()
    plt.plot(epoch_range,train_loss_epoch,'b')
    plt.plot(epoch_range,test_acc_epoch,'r')
    plt.show()
    
    y = Y.data.numpy()
    y_pred = symbol.data.numpy()
    a = np.concatenate((y,y_pred),axis=1)
    
    # Only looking at non zero points
    idxYnoZero = y!=0
    a1 = np.concatenate((np.reshape(y[idxYnoZero],(-1,1)),np.reshape(y_pred[idxYnoZero],(-1,1))),axis=1)
    
    # Plot confusion matrix
    confusion_mat = sklearn.metrics.confusion_matrix(y,y_pred)
    #df_cm = pd.DataFrame(confusion_mat, range(6),range(6))
    plt.figure()
    fig, ax = plt.subplots()
    sns.set(font_scale=1.2) # For label size
    sns.heatmap(confusion_mat, annot=True, cmap="YlGnBu") #Font size
    ax.set_ylim([6,0])
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()
    
    # -----------------------------------------------------------------------
    # Checking results for entire dataset
    correct2 = []
    totaldata_loss = 0
    
    for Y, X in totaldata_loader:
       X = X.view(-1,X.shape[1],X.shape[2])
       # X = X.view(-1,1,X.shape[1])
       Y = Y.view(-1,Y.shape[1])
       X = X.float()
       Y = Y.long()      
       current_batchsize = X.shape[0]
       output = network(X)
       totaldata_loss = totaldata_loss + (F.cross_entropy(output,Y) * current_batchsize/len(totaldata_loader.dataset))
       symbol = output.detach()
       symbol = symbol.max(dim=1).indices
       correct2.append(Y.eq(symbol).numpy())
    correct2 = [ item for sublist in correct2 for item in sublist]   
    correct2 = np.asarray(correct2)  
    correct2 = correct2.reshape(-1)
    accuracy2 = correct2.sum()/len(totaldata_loader.dataset)
        
    print('total data loss is',(totaldata_loss/(len(totaldata_loader.dataset)/batchsize_test)))
    print('total data accuracy is: ',accuracy2)

    y_totaldata = Y.data.numpy()
    y_pred_totaldata = symbol.data.numpy()
    a2 = np.concatenate((y_totaldata,y_pred_totaldata),axis=1)
    
    # Only looking at non zero points
    idxYnoZero = y_totaldata!=0
    a3 = np.concatenate((np.reshape(y_totaldata[idxYnoZero],(-1,1)),np.reshape(y_pred_totaldata[idxYnoZero],(-1,1))),axis=1)
    
    # Plot confusion matrix
    confusion_mat2 = sklearn.metrics.confusion_matrix(y_totaldata,y_pred_totaldata)
    #df_cm = pd.DataFrame(confusion_mat, range(6),range(6))
    plt.figure()
    fig, ax = plt.subplots()
    sns.set(font_scale=1.2) # For label size
    sns.heatmap(confusion_mat2, annot=True, cmap="YlGnBu") #Font size
    ax.set_ylim([6,0])
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()