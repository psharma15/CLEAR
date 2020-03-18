% This algorithm uses a Linear RTI model, similar to Patwari: Radio 
% Tomographic Imaging with Wireless Networks
% The variations of this include Xu (2019): Compressive Sensing Based Radio
% Tomography Imaging with spatial diversity
% Basically - different regularization parameters and other factors - all
% combinations can be easity adapted - as well as FISTA, OMP can be
% implemented to solve this linear inverse problem as well.

% The code firsts develops the inverse problem for RTI

% -------------------------------------------------------------------------
% Pragya Sharma
% ps847@cornell.edu
% 25 Feb 2020
% -------------------------------------------------------------------------

% function [imgMFrti] = rti_lin(yObj,yNoObj,tagPosition,rxPosition,...
%     freq,roomSize,voxelSize,opts)
% Input
% yObj: data with object (complex), in a column
% yNoObj: Data without object (complex), in a column
% Above data are arranged as nTag*nRx*nFreq
% tagPosition: column containing tag coordinates in meters, row is tag Num
% rxPosition: column containing rx cooordinates in meters, row is rx Num
% freq: column containing freq in Hz, row is freq Num
% roomSize: Arranged as [xSt, xEnd; ySt, yEnd; zSt, zEnd] in meters
% voxelSize: [x,y,z] in meters
% opts: TBC
% ---------------------------------------------------------------------
% codePath5 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\RTI';
% addpath(codePath5);

% Reading the saved data from Impinj RFID reader and processing that to generate image by 
% solving y=Ax linear model developed from Scatter model. 

% close all;
clear all;clc
mm = '12'; dd = '05'; yyyy = '2019';
codePath = ['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\',...
    '3-D imaging results Guoyi\TrueModel2019\TrueScaleExperimentData_',mm,'_',dd,'_',yyyy,'\'];
dataPath = ['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\',...
    '3-D imaging results Guoyi\TrueModel2019\TrueScaleExperimentData_',mm,'_',dd,'_',yyyy,'\'];
opts.savePath = ['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\',...
    'Algorithms\ExptData\Guoyi_exptTrueSize_NovDec2019\',mm,dd,yyyy,'\'];
codePath1 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\MP';
codePath2 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\Fourier';
codePath3 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\Cluster';
codePath4 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\LeastSquare';
codePath5 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\RTI';
codePath5 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\RTI';
addpath(codePath1);addpath(codePath2);addpath(codePath3);addpath(codePath4); addpath(codePath5);
opts.saveFig = 0; fSz=14;

addpath(codePath);addpath(dataPath);
fileName = 'Guoyi_x6y6';
% Load the data.
load([dataPath,'data_wo_',fileName,'_',yyyy,mm,dd,'.mat']); % Without object data
data_wo = [chindlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];
load([dataPath,'data_w_',fileName,'_',yyyy,mm,dd,'.mat']); % With object data
data_w = [chindlist, tagindexlist, antennalist, rssiimpinjlist, rssiimpinjlist_d, phasedeglist];
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist msgfreqlist

% Read tag and receiver positions and frequency from the pre-processed
% data corresponding to data arranged in data_w and data_wo
[tagPosition, rxPosition, freq] = tag_antenna_positions3D_func();
freq = freq(:); % Make sure this is a column vector
if size(tagPosition,2) ~= 3 || size(rxPosition,2) ~=3
    fprintf('Sizes of tag and rx positions are not in expected format. Returning...\n');
end

    
nFreq = length(freq);
nTag = size(tagPosition,1);
nRx = size(rxPosition,1);
nData = nTag*nRx*nFreq;
% Room size in meters, each row indicates x, y and z limits
roomSize = [0,3.6;0,3.6;0,2.4]  ; % [0,3.8;0,3.8;0,2.4] [0,3.6;0,3.6;0,2.4] %1.2, 2.0; 1.2, 2.0; 0, 2, [-1,4.0;-1,4.0;-1,3.0]
% voxel size in meters, each row indicates x, y and z direction values
voxelSize = [0.12;0.12;0.5]; % [0.12;0.12;0.32]
xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
nVoxel = [length(xVoxel), length(yVoxel), length(zVoxel)];

%% Select frequencies/ tags/ receivers
%{
freqIdx = 6:8:46;
sParamObj = sParamObj(:,:,freqIdx);
sParamNoObj = sParamNoObj(:,:,freqIdx);
phStdDevNoObj = phStdDevNoObj(:,:,freqIdx);
phStdDevObj = phStdDevObj(:,:,freqIdx);
rssiStdDevObj = rssiStdDevObj(:,:,freqIdx);
rssiStdDevNoObj = rssiStdDevNoObj(:,:,freqIdx);
freq = freq(freqIdx);
nFreq = length(freq);
nData = nTag*nRx*nFreq;
% nTag = ceil(nTag*0.8); % Randomly keep only 80% of the tags
% tagIdx = randperm(size(tagPosition,1),nTag);
% G_w_s = G_w_s(tagIdx,:,:);
% G_wo_s = G_wo_s(tagIdx,:,:);
% tagPosition = tagPosition(tagIdx,:);
%}

%%

% Pre-processing included making both sParamObj and sParamObj 0 when not
% read.
idxEmpt = find(sParamObj == 0); % SparamObj = 0 if not read
numEmpt = length(idxEmpt);
fprintf('Percentage Tx-Rx-Freq pair lost: %3.2f%% \n',numEmpt*100/numel(sParamObj));

figure
plot(1:nData,rssiStdDevNoObj(:)); hold on;
plot(1:nData,rssiStdDevObj(:));
plot(idxEmpt,ones(length(idxEmpt),1),'*');
ylabel('RSSI Std Dev');
xlabel('Data number');
legend('No Object','With Object','No data read');

figure
plot(1:nData,phStdDevNoObj(:)); hold on;
plot(1:nData,phStdDevObj(:));
plot(idxEmpt,10.*ones(length(idxEmpt),1),'*');
ylabel('Phase Std Dev (\circ)');
xlabel('Data number');
legend('No Object','With Object','No data read');

%%
% Applying a phase standard deviation threshold and/or RSSI std dev thresh
opts.tagStdDevCorr = 1; % Tag Phase and RSSI response has variation over time, correcting that

if opts.tagStdDevCorr == 1
    opts.phStdTh = 12; % In degrees. A very high value means no rejection. 
    opts.rssiStdTh  = 1.2*max(max(rssiStdDevNoObj(:)),max(rssiStdDevObj(:)));
    fprintf('Using RSSI Std Dev threshold: %3.2f,\nPhase Std Dev threshold: %3d\260 \n', opts.rssiStdTh,opts.phStdTh);
    idxPhRssiHighStd = (phStdDevNoObj > opts.phStdTh) | (phStdDevObj > opts.phStdTh)| ...
        (rssiStdDevNoObj > opts.rssiStdTh) | (rssiStdDevObj > opts.rssiStdTh);   
    sParamObj(idxPhRssiHighStd) = 0;
    sParamNoObj(idxPhRssiHighStd) = 0;
    idxEmpt = find(sParamObj == 0); % SparamObj = 0 if not read
    fprintf('Percentage Tx-Rx-Freq pair lost: %3.2f%% \n',length(idxEmpt)*100/numel(sParamObj));
end

yObj = sParamObj;
yNoObj = sParamNoObj;

yObj = abs(yObj).*1e-3; % Converting from mW to W, with only amplitude
yNoObj = abs(yNoObj).*1e-3; 

nFreq = length(freq);
nTag = size(tagPosition,1);
nRx = size(rxPosition,1);

% -----------------------------------------------------------------------
% Example plots of a few cases, looking at different frequencies
%{
nFig = [2,2];
idx = 1;
a1 = reshape(yNoObj,nTag*nRx,[]);
b1 = reshape(yObj,nTag*nRx,[]);
figure('Position',[200,200,800,550]);
subplot(nFig(1),nFig(2),1);
plot(a1(idx,:),'o-'); hold on; plot(b1(idx,:),'*-');legend('No Obj','Obj');
subplot(nFig(1),nFig(2),2);
plot(a1(idx+1,:),'o-'); hold on; plot(b1(idx+1,:),'*-');legend('No Obj','Obj');
subplot(nFig(1),nFig(2),3);
plot(a1(idx+2,:),'o-'); hold on; plot(b1(idx+2,:),'*-');legend('No Obj','Obj');
subplot(nFig(1),nFig(2),4);
plot(a1(idx+4,:),'o-'); hold on; plot(b1(idx+4,:),'*-');legend('No Obj','Obj');
%}
opts.model = 'Ellipse';
opts.lambda = 0.15; % 15 cm weight 
W = genW(tagPosition,rxPosition,freq, roomSize,voxelSize,opts);

% -----------------------------------------------------------------------
% Using frequency information, multiple options:
% Take a simple mean, across all frequencies
% other possible variations like Xu with spatial diversity 
% or just use one frequency - try the variations and see what works best.
% Use all frequencies
opts.freqSel = 'mean';
switch(opts.freqSel)
    case 'mean'
        yObjMean = zeros(nTag,nRx);
        yNoObjMean = zeros(nTag,nRx);
        for i = 1:nTag
            for j = 1:nRx
                idxNonZeroObj = yObj(i,j,:) ~=0;
                if sum(idxNonZeroObj)
                    yObjMean(i,j) = mean(yObj(i,j,idxNonZeroObj),3);
                else
                    yObjMean(i,j) = 0;
                end
                idxNonZeroNoObj = yNoObj(i,j,:)~=0;
                if sum(idxNonZeroNoObj)
                    yNoObjMean(i,j) = mean(yNoObj(i,j,idxNonZeroNoObj),3);
                else
                    yNoObjMean(i,j) = 0;
                end
            end
        end

        yObjTagRx = yObjMean(:);
        yNoObjTagRx = yNoObjMean(:);
    case 'fixed'
        % Use one fixed frequency
        freqDesired = 925e6;
        freqIdx = find(freq<=freqDesired,1,'last');
        if isempty(freqIdx)
            warning('Frequency outside available range')
            return
        end
        fprintf('Frequency selected: %f\n',freq(freqIdx));
        yObjTagRx = yObj(:,:,freqIdx);
        yNoObjTagRx = yNoObj(:,:,freqIdx);
        yObjTagRx = yObjTagRx(:);
        yNoObjTagRx = yNoObjTagRx(:);
        
    case 'all'
        % Use all the frequencies
        fprintf('RTI: Using all frequencies individually.\n');
        yObjTagRx = yObj(:);
        yNoObjTagRx = yNoObj(:);
        W = repmat(W,length(freq),1);
    otherwise
        fprintf('Enter correct method for using frequency information.\n');
end

idxNonZero = (yObjTagRx~=0) & (yNoObjTagRx~=0);
y_rssi = 10.*(log10(yObjTagRx(idxNonZero)) - log10(yNoObjTagRx(idxNonZero))); % Taking only the RSS in dB

W = W(idxNonZero,:);
% figure
% imagesc(W)
% title('W'); xlabel('Voxel');ylabel('Signal Link');
% Now different methods to solve.

%% MF
fSz1=16;
[~,imgMFrti,tCompMF] = matchFilt(W,y_rssi,roomSize,voxelSize);
title(['MF RTI, T(ms): ',num2str(round(tCompMF*1000),'%.0f'),' Frequency: ',opts.freqSel],'FontSize',fSz1)

[maxImgMF,maxIdx] = (max(abs(imgMFrti)));
pkVoxel = xyzVoxelCoord(maxIdx,:);
fprintf('******* Max MF image: %f *********\n',maxImgMF);

threshRatio = 0.88;
thresh = maxImgMF*threshRatio;
imgMFabs = reshape(abs(imgMFrti),nVoxel);
imgMFabs(imgMFabs<thresh) = 0;
imgMFabs(imgMFabs>=thresh) = 1;

% opts.saveFig = 1;
imgMFrti(abs(imgMFrti)<thresh) = 0;
visImg(imgMFrti,roomSize,voxelSize);
title(['MF RTI, T(ms): ',num2str(round(tCompMF*1000),'%.0f'),' Frequency: ',opts.freqSel,', Th: ',num2str(threshRatio)],'FontSize',fSz)
% view(90,0)
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_mf_calib',num2str(opts.calibType),'.fig']);
end

[~, clusters,a1] = i4block_components(imgMFabs, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.3; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.2;
opts.elemNumTh = 0.4; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end

%% OMP
opts.RRfig = 0;
[~,imgOMPcomp,tCompOMP] = voxelOMP(W,y_rssi,opts,roomSize,voxelSize,opts.savePath);


%% FISTA
opts2.ifWavelet = 0; opts2.wavName = 'Daubechies';opts.par =4;
opts2.nVoxel = nVoxel;
if opts2.ifWavelet
    addpath(genpath('C:\Program Files\MATLAB\R2018a\toolbox\Wavelab850'));
    WavePath();
end

% ***************************************************************************
% FISTA is written for complex, make sure it works for this.
[~,imgFista,tCompFista] = lsFISTA_l1(W,y_rssi,roomSize,voxelSize,opts2);
title(['FISTA RTI, T(ms): ',num2str(round(tCompFista*1000),'%.0f')],'FontSize',fSz1)

% end

