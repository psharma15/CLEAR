%% 
% Solves b = Ax using orthogonal matching pursuit (OMP), where x is the
% image in the voxel domain. The current stopping criteria is number of
% voxels being occupied with the object.
% -------------------------------------------------------------------------
% Last modified on: 04 Jan 2018
% Pragya Sharma, ps847@cornell.edu
% -------------------------------------------------------------------------

function [imgBrightness,imgComplex,tComp] = voxelOMP(A,b,opts,roomSize,...
    voxelSize,savePath)
%% 
% Defining options
opts.sc = 0.002; % Fractional sparsity (percent value is opts.sc*100
opts.scPercent = 1;
opts.norm = 0;
opts.regConst = 1e5; %5e-2

% Starting OMP
tic;
[x,idxSet,~,iter,res] = omp(A,b,opts);
tComp = toc;

idxSet(idxSet~=0) = x;
imgComplex = idxSet;

resNorm = norm(res(end));
fprintf('Voxel OMP stopped at iter = %d,\n with residual norm = %g,\n in %3.2f s.\n',iter,...
    resNorm,tComp);

% save([savePath,'\voxOMP_img',dataName,'.mat'],'imgComplex', '-v7.3')

%% Visualize image

imgBrightness = visImg(imgComplex,roomSize,voxelSize);

a = alphamap('rampup',256);
imgThresh = 0;
a(1:imgThresh)=0;
alphamap(a); 
% 
% title(['Voxel OMP, #iter = ',num2str(iter),', SC = ',...
%     num2str(opts.sc),', resNorm = ',num2str(resNorm),', imgTh = ',...
%     num2str(imgThresh),', time = ',...
%     num2str(tComp),' s'],'FontSize',12)

title(['Voxel OMP, #iter = ',num2str(iter),', resNorm = ',num2str(resNorm),...
    ', time = ', num2str(tComp),' s'],'FontSize',12)
end

