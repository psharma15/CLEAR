% This is object detection code using data from CST in 3D 
% Setup is same as: "Ubiquitous tagless object locating with ambient 
% harmonic tags," Yunfei Ma and Edwin C. Kan, The 35th Annual IEEE 
% International Conference on Computer Communications, San Francisco, CA, 
% USA, Apr. 2016.
% But instead of Fourier, Least squares has been used to obtain the image,
% using Conjugate Gradient Least Squares.
% Pragya Sharma (ps847@cornell.edu)
% October 14, 2018

function imgBrightness = ...
    lsCG(expConst,sParamCalib,roomSize,voxelSize,savePath,dataName)

%% Image Inverse Reconstruction - Complex 
% Image: brightness as a function of each voxel using CGLS
% A = expConst, b = sParamCalib, regConst = shift
regConst = 0.3;
maxIter = 2;

% load([savePath,'\lsMats',dataName,'.mat'])
tic;
% Using CGLS algorithm from following:
%   01 Sep 1999: First version.
%                Per Christian Hansen (DTU) and Michael Saunders (visiting
%                DTU).

[imgReconstructComplex,flag,resNE,iter] = cgls(expConst,sParamCalib,regConst,[],maxIter,[],[]);
tComp = toc;

fprintf('Time: %3.2f s, Flag: %d, Rel. residual: %f, #Iter: %d\n',tComp,flag,resNE,iter);

save([savePath,'\lsCGimg',dataName,'.mat'],'imgReconstructComplex', '-v7.3')
% clear expConst 

%% Visualizing reconstructed image
imgBrightness = visImg(imgReconstructComplex,roomSize,voxelSize);

%% 
imgThresh = 120;
a = alphamap('rampup',256);
a(1:imgThresh)=0;
alphamap(a); 

title(['CGLS, #iter = ',num2str(iter),', reg = ',...
    num2str(regConst),', imgTh = ',num2str(imgThresh),', time = ',...
    num2str(tComp),' s'],'FontSize',14)

end


