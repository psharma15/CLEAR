%% 
% Solve J(f) = ||b-Ax|| + a((||x||^1)_1) 
% This code is using FISTA algorithm to solves above LASSO problem.
% -------------------------------------------------------------------------
% Created on: 13 Nov 2018
% Last modified on: 
% Pragya Sharma, ps847@cornell.edu
% -------------------------------------------------------------------------

function [imgBrightness,imgComplex,tComp] = lsFISTA_l1(A,b,roomSize,...
    voxelSize,opts)
%% 
% Defining options
opts.maxIter = 10; %20
% Gamma: g(x)=gamma*||x||_1, ie, weightage of the sparsity term
% Order of 100 for simulation.. 5e-2 for expt true scale, 10,0.05
opts.gamma = max(abs(A'*b))*0.9; %0.9
% opts.gamma = 5;
opts.gammaScale =0.9;% 0.9, Each iteration gamma is scaled by gammaScale
opts.verbose = 1;
opts.tol = 1e-8; %15e-8
opts.backtracking = 0;
% L0: Lipschitz Constant (initialization for backtracking, fixed otherwise)
% Keep L greater than max(eig(ATA)) so that step size is between (0,1/L]
% 7.7e4, 1.6e5 max(eig(A'*A)) - check! 8.051e5
%4e5 for all the experiments and sim with 6 freq, Above 6e5 for 50 freq (may change with ph thresh)
opts.L0 = 9e5; % With wavelet, it can be low 1.2e5

% opts.xinit = abs((A'*b)./norm(A)); % Poor starting. Did not converge with large iterations.


maxFistaAttempt = 5;
numFistaAttempt = 0;
imgComplex = 0;

while (numFistaAttempt<=maxFistaAttempt) && (sum(imgComplex)==0)
    % Starting FISTA
    fprintf('-------------- CHECK THE FISTA VERSION. --------------\n');

    tic;
    [imgComplex,iter,relNorm,L] = fista_lasso_v3(A,b,opts);
    tComp = toc;
    
    if sum(imgComplex)==0
        fprintf('Fista convergence failed, restarting with different intial point...\n');
        if numFistaAttempt == 2
            opts.xinit = abs((A'*b)./norm(A));
        else
            opts.xinit = randn(size(A,2),1); % Poor starting. Did not converge with large iterations.
    
        end
    end
    numFistaAttempt = numFistaAttempt + 1;
end

   
fprintf('FISTA stopped at iter = %d with relNorm = %g in %3.2f s.\n',iter,...
    relNorm,tComp);

% save([savePath,'\lsFISTA_l1_img',dataName,'.mat'],'imgComplex', '-v7.3')

%% Visualize image
fprintf('Isnan result: %d\n',sum(isnan((imgComplex(:)))));
% imgBrightness = [];
imgBrightness = visImg(imgComplex,roomSize,voxelSize);

a = alphamap('rampup',256);
imgThresh = 40;
a(1:imgThresh)=00;
alphamap(a); 

title(['FISTA l1, #iter = ',num2str(iter),', reg = ',...
    num2str(opts.gamma),', L = ',num2str(L),', imgTh = ',...
    num2str(imgThresh),', time = ',...
    num2str(tComp),' s'],'FontSize',12)
% 
end

