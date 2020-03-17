% This is object detection code using data from CST in 3D 
% Algorithm source is: "Ubiquitous tagless object locating with ambient 
% harmonic tags," Yunfei Ma and Edwin C. Kan, The 35th Annual IEEE 
% International Conference on Computer Communications, San Francisco, CA, 
% USA, Apr. 2016.
% Pragya Sharma (ps847@cornell.edu)
% 15 Nov 2018

function [imgBrightness,imgComplex,tComp] = matchFilt(A,b,roomSize,...
    voxelSize)
%% 

% Starting matched filtering approximately - not same as above paper. 
tic;
imgComplex = A'*b; % Simplified implementation
tComp = toc;

fprintf('Time for matched filtering = %4.3f s.\n',tComp);

maxImgInt = (max(abs(imgComplex)));
fprintf('Maximum intensity value is = %g. \n',maxImgInt);
% save([savePath,'\mfilt_img',dataName,'.mat'],'imgComplex', '-v7.3')

%% Visualize image

imgBrightness = visImg(imgComplex,roomSize,voxelSize);
% [imgBrightness,g] = visImg_grad(imgComplex,roomSize,voxelSize);

a = alphamap('rampup',256);
imgThresh = 180;
a(1:imgThresh)=0;
alphamap(a); 

title(['Matched filtering, imgTh = ',...
    num2str(imgThresh),', time = ',...
    num2str(tComp),' s'],'FontSize',12)




end
