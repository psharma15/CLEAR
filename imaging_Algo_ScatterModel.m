
% -----------------------------------------------------------------------------------------
% RF Scatter model is developed using paper: Y. Ma and E. C. Kan, “Ubiquitous tagless object locating 
% with ambient harmonic tags,” in IEEE International Conference on Computer Communications 
% (INFOCOM), 2016.
% -----------------------------------------------------------------------------------------

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
% This is non reversible step, can't go back and rerun this section
% fprintf('************Discarding RSSI values...******************\n');
% 
% sParamObjAbs = abs(sParamObj);sParamObjAbs(sParamObjAbs==0) = 1;
% sParamNoObjAbs = abs(sParamNoObj); sParamNoObjAbs(sParamNoObjAbs==0) = 1;
% sParamObj = sParamObj./sParamObjAbs;
% 
%% Add noise and average it.
% snr = 20;
% nAvg = 10000;
% sParamObjAvg = sParamObj;
% sParamNoObjAvg = sParamNoObj;
% for iter = 1:nAvg-1
%     [sParamNoisyObj,sParamNoisyNoObj,calcSNRObj,calcSNRNoObj] = addNoise(sParamObj,...
%     sParamNoObj,snr);
%     sParamObjAvg = sParamObjAvg+sParamNoisyObj;
%     sParamNoObjAvg = sParamNoObjAvg+sParamNoisyNoObj;
% end
% sParamNoObj = sParamNoObjAvg/nAvg;
% sParamObj = sParamObjAvg/nAvg;

%% K-space view and processing
%{
posScat = [1.8,1.8,1.2]; % Considering middle of the capture volume
opts3 = [];
[K,~] = kspaceXYZ(tagPosition, rxPosition, freq, posScat,opts3);
%}

%% Inverse object reflectivity reconstruction for each voxel
opts.genA = 1; opts.loadA = 0;
opts.saveA = 0;
opts.calibType = 1; %4.2
opts.RssiThresh = 0.7; %
opts.phStdDevCorr = 0;
opts.phStdDev = phStdDevObj;

[A,b1] = ...
    genExptAb(sParamObj,sParamNoObj,tagPosition,rxPosition,freq,...
    roomSize,voxelSize,opts);
opts.genA=0;opts.loadA=0;
opts.calibType = 4.2;
[~,b2] = ...
    genExptAb(sParamObj,sParamNoObj,tagPosition,rxPosition,freq,...
    roomSize,voxelSize,opts);
opts.calibType = 0;
[~,b3] = ...
    genExptAb(sParamObj,sParamNoObj,tagPosition,rxPosition,freq,...
    roomSize,voxelSize,opts);
opts.ensembCalib = 0;
fprintf('*****Check FINAL calibration type.*****\n');
if opts.ensembCalib == 1
    fprintf('Final: Ensemble of calibration 1 and 4.2.\n');
    b = (b1/norm(b1))+(b2/norm(b2));
    opts.calibType = 'En';
else
    b = b2; opts.calibType = 4.2; fprintf('Final: Calibration 4.2.\n');
%     b = b1; opts.calibType = 1; fprintf('Final: Calibration 1.\n');
%     b = b3; opts.calibType = 0; fprintf('Final: Calibration 0.\n');

end
opts.noCalib = 0; % 1 if No calibration is desired
if opts.noCalib
    fprintf('******* DISCARDING any Calibration. *******\n');
    b = sParamObj(:);
end

clearvars xVoxel yVoxel zVoxel
[xyzVoxelCoord,~,~,~] = genXYZ(roomSize,voxelSize);

for i = 1:size(A,1)
    if b(i) == 0
        A(i,:) = 0; % Just rejecting data with 0 b signal
    end
end

close all;

%% Phase standard deviation correction: less weigtage if high std dev with 
% or without object. 
% When no phase threshold 
is used earlier.
% opts.phStdDevCorr = 0;
if opts.phStdDevCorr 

    phStdDev = phStdDevObj(:);
    phStdDevNoObj = phStdDevNoObj(:);

    phStdDev(phStdDev == 0) = 1; % Replace 0 std dev by 1
    % phStdDevNoObj(phStdDevNoObj == 0) = 1;

    fprintf('Scaling A with std dev. \n');

    for i =  1:size(A,1)
    %     A(i,:) = A(i,:)./(phStdDevObj(i).*phStdDevNoObj(i));
        A(i,:) = A(i,:)./(phStdDev(i));
    end
    % b = b./(phStdDevObj.*phStdDevNoObj);
%     b = b./(phStdDev);
    
end
 
%% ----------------------------------------------------------------------------------------------------
% Matched filtering (MF) algorithm

[imgMFnorm,imgMFcomp,tCompMF] = matchFilt(A,b,roomSize,voxelSize);
title(['MF, T(ms): ',num2str(round(tCompMF*1000),'%.0f'),', Calib: ',num2str(opts.calibType)],'FontSize',12)

[maxImgMF,maxIdx] = (max(abs(imgMFcomp)));
pkVoxel = xyzVoxelCoord(maxIdx,:);
fprintf('******* Max MF image: %f *********\n',maxImgMF);

threshRatio = 0.8;
thresh = maxImgMF*threshRatio;
imgMFabs = reshape(abs(imgMFcomp),nVoxel);
imgMFabs(imgMFabs<thresh) = 0;
imgMFabs(imgMFabs>=thresh) = 1;

% opts.saveFig = 1;
imgMFcomp(abs(imgMFcomp)<thresh) = 0;
visImg(imgMFcomp,roomSize,voxelSize);
title(['MF, T(ms): ',num2str(round(tCompMF*1000),'%.0f'),', Cal: ',num2str(opts.calibType),', Th: ',num2str(threshRatio)],'FontSize',fSz)
% view(90,0)
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_mf_calib',num2str(opts.calibType),'.fig']);
end

[~, clusters,a] = i4block_components(imgMFabs, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.3; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.61; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0.6; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end

%% ----------------------------------------------------------------------------------------------------
% Voxel based Orthogonal Matching Pursuit (OMP)
opts.RRfig = 0;
% Asing = single(A); % bsing = single(b);
[~,imgOMPcomp,tCompOMP] = voxelOMP(A,b,opts,roomSize,voxelSize,opts.savePath);
threshRatioOMP = 0;
title(['OMP, T(ms): ',num2str(round(tCompOMP*1000),'%.0f'),', Cal: ',num2str(opts.calibType),', Th: ',num2str(threshRatioOMP)],'FontSize',fSz)
% view(90,0)
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_omp_calib',num2str(opts.calibType),'.fig']);
end

[maxImgOMP,maxIdxOMP] = (max(abs(imgOMPcomp)));
fprintf('******* Max OMP image: %f *********\n',maxImgOMP);

pkVoxelOMP = xyzVoxelCoord(maxIdxOMP,:);

imgVoxOMPabs = reshape(abs(imgOMPcomp),nVoxel);
imgVoxOMPabs(imgVoxOMPabs>0) = 1;

[ c, clusters,~ ] = i4block_components(imgVoxOMPabs,roomSize,voxelSize );
fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.3; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.61; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0.6; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

for i = 1:size(clusterOut.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusterOut.centroid(i,:),clusterOut.elemNum(i));
end

%% ----------------------------------------------------------------------------------------------------
% Dictionary based Orthogonal Matching Pursuit (OMP)
% There are 3 dictionary files: dict1 (Object signal, no calibration),
% dictCalib1(Calibration 1), dictCalib42 (Calibration 4.2)
load([opts.dictFile,'.mat']);

if length(b) ~= size(dict1,1)
    warning('Error in dictionary size. Stopping DictOMP...');
else
    if opts.calibType == 1
        dict = dictCalib1; fprintf('Dict: Calib 1.\n');
    elseif opts.calibType == 4.2
        dict = dictCalib42; fprintf('Dict: Calib 4.2.\n');
    elseif opts.calibType == 0
        dict = dictCalib0; fprintf('Dict: Calib 0.\n');
    end
    if opts.ensembCalib == 1
        dict = dictCalib1 + dictCalib42; fprintf('Dict: Calib Ensemble 1, 4.2.\n');
    end
    if opts.noCalib == 1
        dict = dict1; % If no calibration, discard above steps, use no calibrated dict
        fprintf('***** Dict: No Calibration. *****\n');
    end
    opts.dict1Obj = 1;
    if opts.dict1Obj
        fprintf('******** Using dictionary of ONLY 1 object. ********\n');
        dict = dict(:,1:6);
    end
    
    opts.RRfig = 1;
    % Stopping k of dictionary based OMP, has to be LESS THAN SIZE(dict,2)
    opts.sc = size(dict,2);
%     opts.sc = 1;
    opts.kVal = 2; % Number of objects + 1

    [objDictOMP,xDictOMP,normResDictOMP] = dictOMP(dict,b,opts);
end

%% ----------------------------------------------------------------------------------------------------
% Using Conjugate Gradient Least Squares
% imgBrightnessCG = ...
%     lsCGexpt(A,b,roomSize,voxelSize,opts.savePath,fileName);

%% ----------------------------------------------------------------------------------------------------
% Using L1 norm based reqularization and solving with Fast Iterative Shrinkage Thresholding 
% Algorithm (FISTA). See 2009 paper by Beck and Teboulle, "A Fast Iterative Shrinkage-Thresholding 
% Algorithm," for the algorithm details.

for Linear Inverse Problems∗
bNonZeroIdx = b~=0;
ARed = A(bNonZeroIdx,:);
bRed = b(bNonZeroIdx);
opts4.ifWavelet = 0; opts4.wavName = 'Daubechies';opts.par =4;
opts4.nVoxel = nVoxel;
if opts4.ifWavelet
    addpath(genpath('C:\Program Files\MATLAB\R2018a\toolbox\Wavelab850'));
    WavePath();
end

% ***************************************************************************
[imgFista,imgFistaComp,tCompFista] = lsFISTA_l1(ARed,bRed,roomSize,voxelSize,opts4);
% xticks([0.2 0.3 0.4 0.5 0.6 0.7]);
% xticklabels({'0.2','0.3','0.4','0.5','0.6','0.7'});
% savefig([dataPath,'\data_16tx8rx_',dataName,'_fista_freq',num2str(nFreq)]);
% view([-70,14])

[maxImgFista,maxIdxFista] = (max(abs(imgFistaComp)));
fprintf('******* Max FISTA image: %f *********\n',maxImgFista);
pkVoxelFista = xyzVoxelCoord(maxIdxFista,:);

imgFistaAbs = reshape(abs(imgFistaComp),nVoxel);
threshRatioFista = 0.75;
thresh = maxImgFista*threshRatioFista;

imgFistaComp(abs(imgFistaComp)<thresh)=0;
imgBrightness = visImg(imgFistaComp,roomSize,voxelSize);
title(['Fista, T(ms): ',num2str(round(tCompFista*1000),'%.0f'),', Cal: ',num2str(opts.calibType),', Th: ',num2str(threshRatioFista)],'FontSize',fSz)
% view(90,0)
% opts.saveFig = 1;
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_fista_calib',num2str(opts.calibType),'.fig']);
    fprintf('Saving FISTA figure... \n');
end

imgFistaAbs(imgFistaAbs>=thresh) = 1;
imgFistaAbs(imgFistaAbs<thresh) = 0;

[ c, clusters,~ ] = i4block_components(imgFistaAbs,roomSize,voxelSize );
fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.3; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.5; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0.5; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

for i = 1:size(clusterOut.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusterOut.centroid(i,:),clusterOut.elemNum(i));
end

%% ----------------------------------------------------------------------------------------------------
% Ensemble image of MF, OMP and FISTA
nAlgo = 3; % MF, OMP, FISTA
% wt = (1/nAlgo).*ones(nAlgo,1);
wt = [0.25;0.5;0.5];
% imgMFcomp,imgOMPcomp,imgFistaComp
imgMFnorm = wt(1).*((abs(imgMFcomp)-min(abs(imgMFcomp)))./(max(abs(imgMFcomp))-min(abs(imgMFcomp))));
imgOMPnorm = wt(2).*((abs(imgOMPcomp)-min(abs(imgOMPcomp)))./(max(abs(imgOMPcomp))-min(abs(imgOMPcomp))));
imgFISTAnorm = wt(3).*((abs(imgFistaComp)-min(abs(imgFistaComp)))./(max(abs(imgFistaComp))-min(abs(imgFistaComp))));

imEns = imgMFnorm+imgOMPnorm+imgFISTAnorm;
imgBrightness = visImg(imEns,roomSize,voxelSize);
title(['Ensemble Image, calib',num2str(opts.calibType)],'FontSize',fSz);
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_Ensemble_calib',num2str(opts.calibType),'.fig']);
    fprintf('Saving Ensemble figure... \n');
end
% rmpath(codePath1);rmpath(codePath2);
%% 
% clearvars -except A

