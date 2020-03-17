% Combining RSSI and phase based methods - ie - RTI and original setup.
% The original equations are y_rssi = Wx, and b = Ax, combining using 
% unknown Lagrangian multiplier lambda.

% codePath5 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\RTI';
% addpath(codePath5);

%% -----------------------------------------------------------------------------------

% Method 1: Using linear combination of both the system of equations:
% y1 + a.y2 = (A1 + aA2).x

yObj = sParamObj;
yNoObj = sParamNoObj;

yObj = abs(yObj).*1e-3; % Converting from mW to W, with only amplitude
yNoObj = abs(yNoObj).*1e-3; 

nFreq = length(freq);
nTag = size(tagPosition,1);
nRx = size(rxPosition,1);

opts.model = 'Ellipse';
opts.lambda = 0.15; % 15 cm weight 
W = genW(tagPosition,rxPosition,freq, roomSize,voxelSize,opts);

% -----------------------------------------------------------------------
% Using frequency information, multiple options:
% Take a simple mean, across all frequencies
% other possible variations like Xu with spatial diversity 
% or just use one frequency - try the variations and see what works best.
% Use all frequencies
opts.freqSel = 'all';
switch(opts.freqSel)
    
    case 'mean'
        fprintf('RTI: Using mean of all the freqeuncies for each tag-rx.\n');
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
        fprintf('RTI: Using one fixed frequency. \n');
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

% Finding lambda by optimizing the condition number of A1
%{
A_norm = A./norm(A); % To make both same size

W = W./norm(W);
% A is from phase based and W is from rti
A1_func = @(lambda) A_norm(idxNonZero,:) + lambda.*W;

lambda = -0.2:0.1:0.2;
nPts = length(lambda);
maxOneInfNormVec = zeros(nPts,1);
twoNormCondVec = zeros(nPts,1);
% oneNormCondVec = zeros(nPts,1);
for i = 1:nPts
    A1 = A1_func(lambda(i));
%     infNormVec(i) = cond(A1,Inf); 
    maxOneInfNormVec(i) = max(max(sum(abs(A1),2)),max(sum(abs(A1),1)));
    twoNormCondVec(i) = cond(A1,2);
%     oneNormCondVec(i) = cond(A1,1);
end
%
figure
plot(lambda,log10(maxOneInfNormVec),'*-')
hold on
plot(lambda,log10(twoNormCondVec),'o-')
plot(lambda,log10(oneNormCondVec),'s-')
legend('Max(Infinity norm, One Norm)','2 Norm Condition Number')
ylabel('log10(Condition Number)')
xlabel('lambda')
title('Condition Number of A_{ph} + lambda*A_{rssi}')
%}

% Setup requires y1 and y2 to be of the same size, does freqSel can only be
% 'all'
if ~strcmp(opts.freqSel,'all')
    warning('Correct the frequency selection');
end
lambda = 0.8;
A1 = (1-lambda).*A(idxNonZero,:)+lambda.*W;

% b is coming from phase after some calibration
% y_rssi coming from rti algorithm
y_ph_rssi = (1-lambda).*b(idxNonZero)./norm(b(idxNonZero)) + lambda.*(y_rssi)./norm(y_rssi);

% -------------------------------------------------------------------------
% Solution is similar as earlier, use inverse problem solutions to estimate
% x. Using MF

fSz1=12;
[~,imgMF_comb,tCompMF_comb] = matchFilt(A1,y_ph_rssi,roomSize,voxelSize);
title(['MF Comb 1 Ph-RTI, T(ms): ',num2str(round(tCompMF_comb*1000),'%.0f'),' Frequency: ',opts.freqSel],'FontSize',fSz1)

[maxImgMF,maxIdx] = (max(abs(imgMF_comb)));
pkVoxel = xyzVoxelCoord(maxIdx,:);
fprintf('******* Max MF image: %f *********\n',maxImgMF);

threshRatio = 0.85;
thresh = maxImgMF*threshRatio;
imgMFabs = reshape(abs(imgMF_comb),nVoxel);
imgMFabs(imgMFabs<thresh) = 0;
imgMFabs(imgMFabs>=thresh) = 1;

% opts.saveFig = 1;
imgMF_comb(abs(imgMF_comb)<thresh) = 0;
visImg(imgMF_comb,roomSize,voxelSize);
title(['MF Comb 1 Ph-RTI, T(ms): ',num2str(round(tCompMF_comb*1000),'%.0f'),' Frequency: ',opts.freqSel,', Th: ',num2str(threshRatio)],'FontSize',fSz)
% view(90,0)
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_mf_comb_ph_rti',num2str(opts.calibType),'.fig']);
end

[~, clusters,a1] = i4block_components(imgMFabs, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.4; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.61; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end

%{
% OMP
opts.RRfig = 0;
[~,imgOMP_comb,tCompOMP_comb] = voxelOMP(A1,y_ph_rssi,opts,roomSize,voxelSize,opts.savePath);

% FISTA
opts2.ifWavelet = 0; opts2.wavName = 'Daubechies';opts.par =4;
opts2.nVoxel = nVoxel;
if opts2.ifWavelet
    addpath(genpath('C:\Program Files\MATLAB\R2018a\toolbox\Wavelab850'));
    WavePath();
end

% ***************************************************************************
% FISTA is written for complex, make sure it works for this.
[~,imgFista_comb,tCompFista_comb] = lsFISTA_l1(A1,y_ph_rssi,roomSize,voxelSize,opts2);
title(['FISTA RTI, T(ms): ',num2str(round(tCompFista_comb*1000),'%.0f')],'FontSize',fSz1)

% end
%}
%}

%% Method 2: Another solution by augmenting the linear system of equations 
% to [y1; y2] = [A1; A2] x

yObj = sParamObj;
yNoObj = sParamNoObj;

yObj = abs(yObj).*1e-3; % Converting from mW to W, with only amplitude
yNoObj = abs(yNoObj).*1e-3; 

nFreq = length(freq);
nTag = size(tagPosition,1);
nRx = size(rxPosition,1);

opts.model = 'Ellipse';
opts.lambda = 0.15; % 15 cm weight 
W = genW(tagPosition,rxPosition,freq, roomSize,voxelSize,opts);

% -----------------------------------------------------------------------
% Using frequency information, multiple options:
% Take a simple mean, across all frequencies
% other possible variations like Xu with spatial diversity 
% or just use one frequency - try the variations and see what works best.
% Use all frequencies
stdThresh = 0.
rssiValThresh = 0.8;
rssiNfreqThresh = 0.8;
opts.freqSel = 'mean';

figure; hold on;
count = 1;
stdRSSIfreq = zeros(nTag,nRx); % Std dev of RSSI ratio obj/noObj over different frequencies.
for i = 1:nTag
    for j = 1:nRx
        idxNonZeroObj = yObj(i,j,:)~=0;
        if sum(idxNonZeroObj)
            stdRSSIfreq(i,j) = std(yObj(i,j,idxNonZeroObj)./yNoObj(i,j,idxNonZeroObj));
        else
            stdRSSIfreq(i,j) = 0;
        end
        
        scatter(count,stdRSSIfreq(i,j));
        count = count+1;
    end
end
xlabel('Data Number'); ylabel('RSSI std deviation');
% Looking at the mean of std dev, not rejecting 0 values - that will bias
% if more in number - assuming few cases for now -ie tag not read entirely
% at any frequency.
fprintf('Mean RSSI deviation over different frequencies: %3.2f\n',mean(stdRSSIfreq(:)));
stdRSSIfreqThresh = 1.25*mean(stdRSSIfreq(:));
opts.stdRSSIreject = 0;

if opts.stdRSSIreject == 1
    numRejectRSSIfreq = sum(stdRSSIfreq(:)>stdRSSIfreqThresh);
    fprintf('Rejecting heavy multipatch channels: STD RSSI Freq.\n');
    fprintf('# Reject channels = %d\n',numRejectRSSIfreq);    
end

switch(opts.freqSel)
    case 'mean'
        fprintf('RTI: Using mean of all the freqeuncies for each tag-rx.\n');
        yObjMean = zeros(nTag,nRx);
        yNoObjMean = zeros(nTag,nRx);
        for i = 1:nTag
            for j = 1:nRx
                idxNonZeroObj = yObj(i,j,:) ~=0;
                if opts.stdRSSIreject 
                    % Further rejection based on deteriorated channels
                    if stdRSSIfreq(i,j) > stdRSSIfreqThresh
                        % Reject this tag-rx pair entirely
                        continue;
                    end
                end                   
                    
                % Is there any way that when we take mean, only take mean
                % over frequency channels with some information?
%                 if sum(yObj(i,j,:)./yNoObj(i,j,:) <= rssiValThresh) < (nFreq*rssiNfreqThresh)
%                     % ie if at least fixed percentage of freqeuncies does
%                     % show a decrese, then good -those probably have a LoS
%                     % blockage to partial extent
%                     if yObj
%                 
%                     
%                 end
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
        fprintf('RTI: Using one fixed frequency. \n');
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

fSz = 12;
A2 = [A(idxNonZero,:); W];

% b is coming from phase after some calibration
% y_rssi coming from rti algorithm
y_ph_rssi_2 = [b(idxNonZero)./norm(b(idxNonZero));  (y_rssi)./norm(y_rssi)];

% -------------------------------------------------------------------------
% Solution is similar as earlier, use inverse problem solutions to estimate
% x. Using MF

fSz1=12;
[~,imgMF_comb2,tCompMF_comb] = matchFilt(A2,y_ph_rssi_2,roomSize,voxelSize);
title(['MF Comb 2 Ph-RTI, T(ms): ',num2str(round(tCompMF_comb*1000),'%.0f'),' Frequency: ',opts.freqSel],'FontSize',fSz1)

[maxImgMF,maxIdx] = (max(abs(imgMF_comb2)));
pkVoxel = xyzVoxelCoord(maxIdx,:);
fprintf('******* Max MF image: %f *********\n',maxImgMF);

threshRatio = 0.85;
thresh = maxImgMF*threshRatio;
imgMFabs = reshape(abs(imgMF_comb2),nVoxel);
imgMFabs(imgMFabs<thresh) = 0;
imgMFabs(imgMFabs>=thresh) = 1;

% opts.saveFig = 1;
imgMF_comb2(abs(imgMF_comb2)<thresh) = 0;
visImg(imgMF_comb2,roomSize,voxelSize);
title(['MF Comb 2 Ph-RTI, T(ms): ',num2str(round(tCompMF_comb*1000),'%.0f'),' Frequency: ',opts.freqSel,', Th: ',num2str(threshRatio)],'FontSize',fSz)
% view(90,0)
if opts.saveFig == 1
    savefig([opts.savePath,fileName,'_mf_comb_ph_rti',num2str(opts.calibType),'.fig']);
end

[~, clusters,a1] = i4block_components(imgMFabs, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.4; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.61; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end


% OMP
opts.RRfig = 0;
[~,imgOMP_comb,tCompOMP_comb] = voxelOMP(A2,y_ph_rssi_2,opts,roomSize,voxelSize,opts.savePath);

% FISTA
opts2.ifWavelet = 0; opts2.wavName = 'Daubechies';opts.par =4;
opts2.nVoxel = nVoxel;
if opts2.ifWavelet
    addpath(genpath('C:\Program Files\MATLAB\R2018a\toolbox\Wavelab850'));
    WavePath();
end

% ***************************************************************************
% FISTA is written for complex, make sure it works for this.
[~,imgFista_comb,tCompFista_comb] = lsFISTA_l1(A2,y_ph_rssi_2,roomSize,voxelSize,opts2);
title(['FISTA RTI, T(ms): ',num2str(round(tCompFista_comb*1000),'%.0f')],'FontSize',fSz1)

%}
