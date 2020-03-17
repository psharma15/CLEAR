

function [expConst,sParamCalib,xyzVoxelCoord] = ...
    genExptAb(sParamObj,sParamNoObj,tagPosition,rxPosition,freq,...
    dataName,roomSize,voxelSize,savePath)
%%
c = 3e8;                                 % Speed of light
nFreq = length(freq);                 % Number of frequencies that we use for reconstruction using IFT
[nTag, ~] = size(tagPosition);                % Number of tags
[nRecv, ~] = size(rxPosition);              % Number of receivers

figure;
scatter3(tagPosition(:,1),tagPosition(:,2),tagPosition(:,3),'ro');
hold on
scatter3(rxPosition(:,1),rxPosition(:,2),rxPosition(:,3),'k*');
axis 'equal'; axis 'tight';
legend('Tags','Receivers');

%% Extract signal at each receiver corresponding to each tag and frequency
% and perform calibration

sParamCalib = zeros(nTag,nRecv,nFreq);

for freqNum = 1:length(freq)
    lamda = c/freq(freqNum);
    for recvNum = 1:nRecv
        for tagNum = 1:nTag
            distTagRx = norm(tagPosition(tagNum,:)-rxPosition(recvNum,:));
            for recvPairNum = 1:nRecv
                if recvPairNum ~= recvNum
                    sParamObjRatio = sParamObj(tagNum,recvNum,freqNum)...
                        /sParamObj(tagNum,recvPairNum,freqNum);
                    sParamNoObjRatio = sParamNoObj(tagNum,recvNum,freqNum)...
                        /sParamNoObj(tagNum,recvPairNum,freqNum);
                    sParamCalib(tagNum,recvNum,freqNum) = ...
                        sParamCalib(tagNum,recvNum,freqNum) + ...
                        (sParamObjRatio/sParamNoObjRatio-1)*...
                        exp(-1j*(2*pi/lamda)*distTagRx);
                end
            end
        end
    end
end


%% Calculating A and b
expConst = genA(roomSize,voxelSize,tagPosition,rxPosition,freq);
sParamCalib = single(sParamCalib(:));
save([savePath,'\matAb',dataName,'.mat'],'expConst','sParamCalib',...
    'roomSize','voxelSize', '-v7.3')

end
