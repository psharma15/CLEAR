%% 
% Solves b = Ax using orthogonal matching pursuit (OMP), where A is the
% dictionary of fixed set of objects at particular positions. This library
% or dictionary is generated using genLibraryBP.m code. 
% The current stopping criteria is number of objects, needs to be updated,
% as this is unknown.
% There is NO IMAGE GENERATED in this case.
% -------------------------------------------------------------------------
% Last modified on: 04 Jan 2018
% Pragya Sharma, ps847@cornell.edu
% -------------------------------------------------------------------------

function [objIdx,x,normRes] = dictOMP(dict,b,opts)
%% 
% Defining options
opts.norm = 1;
opts.regConst = 5e-2;
opts.scPercent = 0;

% Starting OMP
tic;
[x,~,objIdx,iter,normRes] = omp(dict,b,opts);
tComp = toc;

fprintf('Dictionary OMP stopped at iter = %d, \nwith res norm = %g, \nin %4.3f s.\n',iter,normRes(end),tComp);
fprintf(' The detected case is %d.\n',objIdx);
% save([savePath,'\voxOMP_img',dataName,'.mat'],'imgComplex', '-v7.3')

end

