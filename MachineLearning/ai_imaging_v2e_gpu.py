
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:39:17 2020
@author: Pragya Sharma
Update of ai_imaging_v1: Wed Mar  4 13:46:52 2020
------------------------------------------------------------------------------
Code for AI Imaging - v3, using pre-calculated 8 feature vectors - 
feat_label_setup1_v3
------------------------------------------------------------------------------
This code is dividing data using location tags such as there is no overlap in 
locations of the training and testing data.
13 May: Adding data normalization through custom dataset

15 May 2020: This is the only correctly working code compatible with estimating 
test accuracies at intermediate epochs.

17 May 2020: Adding option to work with intermediate tags.

12 June 2020: Adding GPU capability
"""

# %% All Imports
import sys
import math
import heapq
import numpy as np
import matplotlib
from matplotlib import style
import torch
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, random_split
from sklearn import preprocessing 
import time
import scipy.io
import sklearn.metrics
import seaborn as sns
import random

np.set_printoptions(precision=2)
random.seed(1)
torch.manual_seed(1)
torch.cuda.manual_seed(1)
np.random.seed(1)

print('__Python VERSION:', sys.version)
print('__pyTorch VERSION:', torch.__version__)
print('__CUDNN VERSION:', torch.backends.cudnn.version())
print('__Number CUDA Devices:', torch.cuda.device_count())
print('__Devices')
print('Active CUDA Device: GPU', torch.cuda.current_device())
print ('Available devices ', torch.cuda.device_count())
print ('Current cuda device ', torch.cuda.current_device())


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('\n*****************************************************')
print('Using device: ',(device))
print('*****************************************************')

#plt.style.use("ggplot")
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

# %% All classes and functions
# ------------------------------------------------------------------------------------------------------
# Writing custom dataset
class RF_dataset(Dataset):
    def __init__(self, dataX,dataY):
        self.X = dataX
        self.Y = dataY
    def __len__(self):
        return len(self.Y)

    def __getitem__(self, idx):
        #print(self.Y[idx])
        return (self.Y[idx],self.X[idx])

# ------------------------------------------------------------------------------------------------------
#Network definition
class occupancyCounting(nn.Module):
    def __init__(self):
        super(occupancyCounting,self).__init__()
        n_CW = 120 # This directly affects accuracy if value is low.

        nFeat = 4
        
        
        # tag skip = 1
        self.conv1 = nn.Conv1d(nFeat,n_CW,kernel_size=5, stride=4, padding=1) #featvec
        self.conv2 = nn.Conv1d(n_CW,int(n_CW/2) , kernel_size=64, stride=35, padding=1) #featvec
        self.conv3 = nn.Conv1d(int(n_CW/2),6,kernel_size=24, stride=10, padding=0) #featvec
        self.maxPool1 = nn.MaxPool1d(kernel_size=8,stride=4,padding=0)
        self.avgPool1 = nn.AvgPool1d(kernel_size=5)       
        self.drop1 = nn.Dropout(p=0.2)
        """
        # Tag skip = 2
        self.conv1 = nn.Conv1d(nFeat,n_CW,kernel_size=4, stride=2, padding=1) #featvec
        self.conv2 = nn.Conv1d(n_CW,int(n_CW/2) , kernel_size=8, stride=4, padding=1) #featvec
        self.conv3 = nn.Conv1d(int(n_CW/2),6,kernel_size=24, stride=16, padding=0) #featvec
        self.maxPool1 = nn.MaxPool1d(kernel_size=16,stride=8,padding=0)
        self.avgPool1 = nn.AvgPool1d(kernel_size=5)       
        self.drop1 = nn.Dropout(p=0.2)
        
        
        # tag skip = 4
        self.conv1 = nn.Conv1d(nFeat,n_CW,kernel_size=2, stride=1, padding=1) #featvec
        self.conv2 = nn.Conv1d(n_CW,int(n_CW/2) , kernel_size=8, stride=4, padding=1) #featvec
        self.conv3 = nn.Conv1d(int(n_CW/2),6,kernel_size=24, stride=16, padding=0) #featvec
        self.maxPool1 = nn.MaxPool1d(kernel_size=16,stride=8,padding=0)
        self.avgPool1 = nn.AvgPool1d(kernel_size=5)       
        self.drop1 = nn.Dropout(p=0.2)
        """
        
    def forward(self, x):
        
        # tag skip 1
        x = (F.relu(self.conv1(x)))
        #print(x.shape)
        x =  (F.relu(self.conv2(x))) 
        #print(x.shape)
        x = self.avgPool1(F.relu(self.conv3(x))) # ReLU and pooling operations kind of commute, so order change has very less effect.
        #print(x.shape)
        x = self.drop1(x)
        """
        # Tag skip 2
        x = (F.relu(self.conv1(x)))
        #print(x.shape)
        x =  self.maxPool1(F.relu(self.conv2(x))) 
        #print(x.shape)
        x = self.avgPool1(F.relu(self.conv3(x)))
        #print(x.shape)
        x = self.drop1(x)
        
        # Tag skip 4
        x = (F.relu(self.conv1(x)))
        #print(x.shape)
        x =  self.maxPool1(F.relu(self.conv2(x))) 
        #print(x.shape)
        x = self.avgPool1(F.relu(self.conv3(x)))
        #print(x.shape)
        x = self.drop1(x)
        """
        
        return x
    
# ------------------------------------------------------------------------------------------------------
# Performing weighted selection without replacement for location tags    
def WeightedSelectionWithoutReplacement(weights, m):
    # https://stackoverflow.com/questions/352670/weighted-random-selection-with-and-without-replacement
    elt = [(math.log(random.random()) / weights[i], i) for i in range(len(weights))]
    return [x[1] for x in heapq.nlargest(m, elt)]   
    
# ------------------------------------------------------------------------------------------------------
#%% Main function: Including training and testing model

if __name__ == '__main__':
    
    file_name=r"E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\3-D imaging results Guoyi\TrueModel2019\ResultsFeaturesDict\feat_label_setup1_v6.mat"
    
    data = scipy.io.loadmat(file_name)
    # Data arranged as [batch_size,n_feature,feature_length]
    # In Pytorch images are represented as [batch_siz,channels,height,width]
    # This can be compared to images, with 4 channels and height of feature_length, width of 1
    dataName = 'featVec_2'
    dataX_alltag = data[dataName] # Can change this to featVec, featVec_1, featVec_2, featVec_3
    dataX_alltag = dataX_alltag[:,[0,4,12,13],:]   # Features: Using only RSSI and Phase with object 
    dataY = data['labelVec']    # Labels: #Objects {0,5}
    dataLocTag = data['LocTag'] # Location Tag: Based on pre-defined scheme nData x 5
    dataPosTag = data['PosTag'] # Posture Tag: 1=stand, 2=sit, for three people looks like 211
    datanStand = data['nStand'] # Number of standing people
    datanSit = data['nSit']     # Number of sitting people
    nData = dataX_alltag.shape[0]
    nFeat = dataX_alltag.shape[1]
    feat_length = dataX_alltag.shape[2]
    
    nTagOriginal = 80
    nRecv = 4
    nFreq = 50
    nTagSkip = 1 # Skip 4 tags
    nTag = int(nTagOriginal/nTagSkip) 
    dataX = np.zeros((nData,nFeat,nTag*nRecv*nFreq))
    
    # ------------------------------------------------------------------------------------------------------
    # Now we are skipping some tags
    for i in range(nData):
        for j in range(nFeat):
            temp = dataX_alltag[i,j,:]
            if dataName == 'featVec':
                # featVec arranged as [All tags... all receivers... all freq]
                temp = temp.reshape((nTagOriginal,nRecv,nFreq))
                temp = temp[::nTagSkip,:,:]            
            elif dataName == 'featVec_1':
                #featVec_1 arranged as [all frequencies... all tags... all receivers]
                temp = temp.reshape((nFreq,nTagOriginal,nRecv))
                temp = temp[:,::nTagSkip,:]
            elif dataName == 'featVec_2':
                #featVec_2 arranged as [all receivers... all tags... all frequencies]
                temp = temp.reshape((nRecv,nTagOriginal,nFreq))
                temp = temp[:,::nTagSkip,:]
            elif dataName == 'featVec_3':
                #featVec_3 arranged as [all tags... all frequencies... all receivers]
                temp = temp.reshape((nTagOriginal,nFreq,nRecv))
                temp = temp[::nTagSkip,:,:]
            else:
                print('\n****Enter Correct dataName****\n')

            dataX[i,j,:] = temp.reshape(1,-1)
    
    data_normalize = False               # Normalize data?

    print(len(dataY))
    
    # ------------------------------------------------------------------------------------------------------
    # Now making the code location-specific training
    locTag = np.concatenate((np.arange(1,28), np.arange(32,34)),axis=0) # Possible Location Tags
    wtLocTag = np.ones(len(locTag)) # Weights for each LocTag
    wtLocTag[12] = 5 #Increase the weight for locTag=13, center location
        
    fractionTrainLocTag = "na" # If not using this input "na"
    #numTrainLocTag = math.ceil(fractionTrainLocTag*len(locTag))
    #pickLocTag = locTag[np.sort(WeightedSelectionWithoutReplacement(wtLocTag,numTrainLocTag))] 
    pickLocTag = [11,12,13,16,17,21,32]
    trainIdx = [] #Indices for training
    idxZeroObj = [] #Indices with no object during test
    for i in range(dataLocTag.shape[0]):
        if dataY[i]!= 0:
            if np.sum(np.isin(dataLocTag[i],pickLocTag))>=math.ceil(dataY[i]/2):# once seeing any one of these locations
                trainIdx.append(i) # This currently does not contain any no object data
        if np.sum(dataLocTag[i]) == 0:
            idxZeroObj.append(i) 
            
    # Add 50% (pre-selected) of the 0 object cases randomply from all 0 object
    idxZeroObjTrain = np.sort(random.sample(idxZeroObj,int(len(idxZeroObj)/2)))        
    trainIdx = np.sort(trainIdx + (np.ndarray.tolist(idxZeroObjTrain)))    
    
    # Now the remaining is the testdata, getting its indices
    testIdx = []
    for i in range(dataLocTag.shape[0]):
        if not(np.isin(i,trainIdx)):
            testIdx.append(i)
    

    # ------------------------------------------------------------------------------------------------------
    # Now performing normalization based on the training data
    # Finding mean of all the features in dataset.X[trainIdx] of shape (#train_data,#feat,#samples_per_feature)
    # This is memory consuming, optimize this so as to not create xreshaped
    if data_normalize == True:
        print('Performing feature normalization')
        #scaler = preprocessing.StandardScaler().fit(dataX[trainIdx].transpose((1,0,2)).reshape((nFeat,-1)).reshape((-1,nFeat)))
        dataXtrain = dataX[trainIdx].transpose(1,0,2)
        dataXtrain = dataXtrain.reshape(nFeat,-1)
        dataMean = np.mean(dataXtrain,axis=1).reshape(nFeat,1)
        dataStd = np.std(dataXtrain,axis=1).reshape(nFeat,1)
        dataX = dataX.transpose(1,0,2)
        dataX = dataX.reshape(nFeat,-1)

        dataX = (dataX-dataMean)/dataStd
        dataX = dataX.reshape(nFeat,nData,feat_length)
        dataX = dataX.transpose(1,0,2)

    
    X_train = dataX[trainIdx]
    X_train = torch.from_numpy(X_train).float().to(device)
    # X_train.is_cuda #True
    Y_train = dataY[trainIdx]
    Y_train = torch.from_numpy(Y_train).long().to(device)
    X_test = dataX[testIdx]
    Y_test = dataY[testIdx]
    X_test = torch.from_numpy(X_test).float().to(device)
    Y_test = torch.from_numpy(Y_test).long().to(device)  
    X = torch.from_numpy(dataX).float().to(device)
    Y = torch.from_numpy(dataY).long().to(device)
    """
    ***************************************************************************
    Whenever DataLoader is used, if shuffle is True, make sure you do not use X 
    and Y from totaldata_loader.dataset - as those are not shuffled.
    Unable to find shuffled indices, so keeping shuffle = False here.
    ***************************************************************************
    """    
    batchsize_train = 50
    trainShuffleIdx = (np.arange(len(X_train)))
    np.random.shuffle(trainShuffleIdx) # Shuffling trainShuffleIdx
    numBatches = math.ceil(len(trainIdx)/batchsize_train)
    # ------------------------------------------------------------------------------------------------------
    # %% Start Training Block
    #hyperparameter definition    
    model = occupancyCounting()
    model.to(device)
    learning_rate = 0.002   # 0.02, 0.008
    momentum = 0.9 # 0.1
    random_seed=1
    torch.backends.cudnn.enabled = False
    torch.manual_seed(random_seed)
    # Original no weight decay
    optimizer = optim.SGD(model.parameters(), lr=learning_rate, momentum=momentum, weight_decay = 1e-3)

    EPOCHS = 300
    train_loss_epoch = []      # Training Loss every epoch
    train_acc_epoch = []
    train_time_epoch = []      # Training time accumulating from the 0th epoch 
    
    test_loss_n_epoch = []     # Test Loss every nth epoch
    test_acc_n_epoch =[]       # Test accuracy every nth epoch
    total_acc_n_epoch = []     # Total accuracy every nth epoch
    count_n_epoch = []         # Counting every nth epoch
    
    plt.figure()
    
    start_time = time.time()
    for epoch in range(EPOCHS):
        # ------------------------------------------------------------------------------------------------------
        train_loss = 0
        correct_train = []
        trainIdxStart = 0
        while trainIdxStart < len(X_train):
            #print(len(Y_train))
            if trainIdxStart + batchsize_train < len(X_train):
                trainIdxStop = trainIdxStart + batchsize_train
            else:
                trainIdxStop = len(X_train) 
            trainIdx_batch = trainShuffleIdx[trainIdxStart:trainIdxStop]
            X_trainbatch = X_train[trainIdx_batch]
            Y_trainbatch = Y_train[trainIdx_batch]
            trainIdxStart = trainIdxStop
            model.train()
            output_train = model(X_trainbatch)                 
            loss = F.cross_entropy(output_train,Y_trainbatch)         # Computing total loss of this batch
            optimizer.zero_grad()                    # Clean-up step for PyTorch
            loss.backward()                          #calculate the gradient decent
            optimizer.step()                         #update the weight
                               
            train_loss = train_loss + (loss*len(trainIdx_batch)/len(trainIdx)) # Loss for each epoch
            
            with torch.set_grad_enabled(False):
                # This is not TRUE accuracy, because there is drop-off layer. 
                model.eval()
                output_train = model(X_trainbatch)
                symbol_train = output_train.data
                symbol_train = symbol_train.max(dim=1).indices
                correct_train.append(Y_trainbatch.eq(symbol_train))
            
            
        train_time_epoch.append(time.time()-start_time)     
           
        correct_train = [ item for sublist in correct_train for item in sublist]   
        train_accuracy = sum(correct_train).to(dtype=torch.float)/len(trainIdx)
        train_acc_epoch.append(float(train_accuracy))
        train_loss_epoch.append(float(train_loss))   # Appending Training loss for the epoch
        #print('Epoch: %i' %(epoch+1),', Train loss: %0.2f' %(float(train_loss)),', Train Accuracy: %0.2f' %train_accuracy)

        """
        if (epoch+1)%10 == 0:
            # Update Plot with some random case every 5 epochs
            plt.clf()
            plt.plot(Y_trainbatch[0], 'bo')
            plt.plot(symbol_train[0], 'rx')
            plt.title('epoch: %i , Training loss is ' %(epoch+1) + '%f'%float(train_loss))
            plt.show()
            plt.pause(0.1)
        """
        #print('Epoch: %i' %(epoch+1),', Train loss: %0.2f' %(float(train_loss)),', Train Accuracy: %0.2f' %train_accuracy)
          
        if (epoch+1)%1 == 0:
            count_n_epoch.append(epoch+1)              # Current Epoch
            # ------------------------------------------------------------------------------------------------------
            #evaluation: # Update test loss & accuracy every n epochs
            #print('Epoch: %i' %(epoch+1),', Train loss: %0.2f' %(float(train_loss)),', Train Accuracy: %0.2f' %train_accuracy)
            
            with torch.set_grad_enabled(False):    
                correct_test = []
                test_loss = 0
                
                #print("Hello")
 
                model.eval()
                output_test = model(X_test)
                test_loss = test_loss + (F.cross_entropy(output_test,Y_test))
                symbol_test = output_test.data
                symbol_test = symbol_test.max(dim=1).indices
                correct_test.append(Y_test.eq(symbol_test))
                correct_test = [ item for sublist in correct_test for item in sublist]   
                test_accuracy = sum(correct_test).to(dtype=torch.float)/len(testIdx)
                print('Epoch: %i' %(epoch+1),', Train loss: %0.2f' %(float(train_loss)),', Train Accuracy: %0.2f' %train_accuracy)
                print('Epoch: %i' %(epoch+1),', Test loss: %0.2f' %(float(test_loss)),', Test Accuracy: %0.2f' %test_accuracy)
                test_loss_n_epoch.append(float(test_loss))
                test_acc_n_epoch.append(test_accuracy)
               
    # ------------------------------------------------------------------------------------------------------
    # Some sanity checks by looking at the data
    y = Y_test.data.cpu().numpy()              # Test data: As taken from above, same shuffled as y_pred.
    y_pred = symbol_test.data.cpu().numpy()    # Predicted symbol
    a = np.concatenate((y,y_pred),axis=1)
    # Only looking at non zero points
    idxYnoZero = y!=0
    a1 = np.concatenate((np.reshape(y[idxYnoZero],(-1,1)),np.reshape(y_pred[idxYnoZero],(-1,1))),axis=1)
    
    # ------------------------------------------------------------------------------------------------------
    # Printing train and test results for each epoch
    # print(model.state_dict())
    print('total training time is',train_time_epoch[-1])
    
    plt.figure()
    plt.plot(np.arange(EPOCHS)+1,train_loss_epoch,'b')
    plt.plot(count_n_epoch,test_loss_n_epoch,'r')
    plt.title("Training & Testing loss")
    plt.xlabel("Epoch")
    plt.show()

    plt.figure()
    plt.plot(np.arange(EPOCHS)+1,train_acc_epoch,'b')
    plt.plot(count_n_epoch,test_acc_n_epoch,'r')
    plt.title('Training & Testing accuracy')
    plt.xlabel("Epoch")
    plt.show()
        
    confusion_mat = sklearn.metrics.confusion_matrix(y,y_pred)
    # Make sure y and y_pred are shuffled similarly
    fig, ax = plt.subplots()
    sns.set(font_scale=1.2) # For label size
    sns.heatmap(confusion_mat, annot=True, cmap="YlGnBu") #Font size
    ax.set_ylim([6,0])
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()
    
    # ------------------------------------------------------------------------------------------------------
    # %% Final evaluation: Checking results for entire dataset
    
    correct2 = []
    totaldata_loss = 0
    with torch.no_grad():
        model.eval()      # Should already be in this state...
        output = model(X)
        totaldata_loss = totaldata_loss + (F.cross_entropy(output,Y) )
        symbol = output.detach()
        symbol = symbol.max(dim=1).indices
        correct2.append(Y.eq(symbol))
        correct2 = [item for sublist in correct2 for item in sublist]   
        totaldata_accuracy = sum(correct2).to(dtype=torch.float)/len(Y)

    # ------------------------------------------------------------------------------------------------------
    # Sanity checks by looking at the data
    y_totaldata = Y.data.cpu().numpy() # This is shuffled if totaldata_loader Dataset has shuffle, so correct implementation
    y_pred_totaldata = symbol.data.cpu().numpy()
    a2 = np.concatenate((y_totaldata,y_pred_totaldata),axis=1)
    # Only looking at non zero points
    idxYnoZero = y_totaldata!=0
    a3 = np.concatenate((np.reshape(y_totaldata[idxYnoZero],(-1,1)),np.reshape(y_pred_totaldata[idxYnoZero],(-1,1))),axis=1)
    
    # ------------------------------------------------------------------------------------------------------
    # Displaying overall statistics
    print('Total data loss is: ', totaldata_loss)
    print('Total data accuracy is: ', totaldata_accuracy)
    
    # Plot confusion matrix
    confusion_mat2 = sklearn.metrics.confusion_matrix(y_totaldata,y_pred_totaldata)
    fig, ax = plt.subplots()
    sns.set(font_scale=1.2) # For label size
    sns.heatmap(confusion_mat2, annot=True, cmap="YlGnBu") #Font size
    ax.set_ylim([6,0])
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()
    
    # ------------------------------------------------------------------------------------------------------
    # %% This is data summary section
    print("************* Printing data summary *************")
    print('Fraction of location tags for training  = ',fractionTrainLocTag)
    print('Number of Location tags for training = ',len(pickLocTag)) #numTrainLocTag
    print('Picked Tag locations = ',pickLocTag)
    print('Train: #0 = ',len(idxZeroObjTrain))
    print('Train: #1 = ',sum(Y_train == 1))
    print('Train: #2 = ',sum(Y_train == 2))
    print('Train: #3 = ',sum(Y_train == 3))
    print('Train: #4 = ',sum(Y_train == 4))
    print('Train: #5 = ',sum(Y_train == 5)) 
    print('******************************')    
    print('Test: #0 = ',sum(Y_test == 0))
    print('Test: #1 = ',sum(Y_test == 1))
    print('Test: #2 = ',sum(Y_test == 2))
    print('Test: #3 = ',sum(Y_test == 3))
    print('Test: #4 = ',sum(Y_test == 4))
    print('Test: #5 = ',sum(Y_test == 5))     
        
    # This is summary, observe what locations are giving wrong result, etc.
    # Folllowing are only valid when totaldata_loader has 'shuffle=False'
    wrongDataIdxTotal = np.where(correct2==False) #Indices when error occurs in entire dataset
    wrongDataLocTag = dataLocTag[wrongDataIdxTotal,:]
    wrongDataLocTag = wrongDataLocTag.reshape((wrongDataLocTag.shape[1],wrongDataLocTag.shape[2]))
    wrongDataPosTag = dataPosTag[wrongDataIdxTotal,:]
    wrongDataPosTag = wrongDataPosTag.reshape((wrongDataPosTag.shape[1],wrongDataPosTag.shape[2]))
    wrongData1Obj = np.where(dataY[wrongDataIdxTotal] == 1)
    wrongDataNumObj = dataY[wrongDataIdxTotal]
    wrongDataNumPosLocTag = np.concatenate((wrongDataNumObj,wrongDataPosTag,wrongDataLocTag),axis=1)
    #wrongDataName = data['dataNames'][wrongDataIdxTotal]
    
    # Emptying cache
    torch.cuda.empty_cache()