%% 
% Solve J(f) = ||g-Tf|| + a((||f||^k)_k) + b((||D|f|||^k)_k)
% This code is implementing approach of the following paper:
% M. Çetin (Cetin), W. C. Karl, "Feature-enhanced synthetic aperture radar 
% image formation based on nonquadratic regularization", IEEE Trans. Image 
% Process., vol. 10, no. 4, pp. 623-631, Apr. 2001.
% -------------------------------------------------------------------------
% savePath = ['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\',...
%     'Algorithms\LeastSquare\ProcessedData'];
% dataName = 'RxXYZ6';
% [expConst,sParamCalib,xyzVoxelCoord] = ...
%     genExptAb(sParamObj,sParamNoObj,tagPosition,rxPosition,freq,...
%     dataName,roomSize,voxelSize,savePath);

% load([savePath,'\lsMats',dataName,'.mat'])
% Room size in meters, each row indicates x, y and z limits
% roomSize = [-0.5, 0.5; -0.5, 0.5; -0.25, 0.25]; 
% voxel size in meters, each row indicates x, y and z direction values
% voxelSize = [0.06;0.06;0.06]; 

function fNew = lsCetin(expConst,sParamCalib,roomSize,...
    voxelSize,savePath,dataName)
%%
xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 

nx = length(xVoxel);
ny = length(yVoxel);
nz = length(zVoxel);

nVoxel = nx*ny*nz;

%%
eps = 1e-6; % Epsilon
k = 0.5; % k-Norm
const1 = 5; % lambda1 const 
const2 = 10; % lambda2 const 
stepSz = 1; % step size
del = 5e-3; % Set this delta threshold stopping condition
delCG = 1e-6; % Stopping for solving conjugate gradient 
maxIter = 10;
maxIterCG = 4;

%% Defining 3-D derivative operator

% Defining dx
tempI = [-eye(ny),eye(ny),zeros(ny,(ny*nx) - (2*ny))];
temp1 = zeros((nx-1)*ny,nx*ny);
countRow = 1;
for i = 1:(nx-1)
    temp1(countRow:countRow+ny-1,:) = circshift(tempI,countRow-1,2);
    countRow = countRow+ny;
end

dx = kron(eye(nz), temp1); 
clearvars tempI temp1 xVoxel yVoxel zVoxel;

% Defining dy
d1 = -eye(ny,ny) + diag(ones(ny-1,1),1);
d1 = d1(1:end-1,:);
dy = kron(eye(nx*nz,nx*nz),d1);
clearvars d1;

% Defining dz
dz = zeros((nz-1)*nx*ny, nVoxel);
dz(1,:) = [-1, zeros(1,nx*ny-1), 1, zeros(1,nVoxel - (nx*ny+1))];
for i = 2:(nx*ny*(nz-1))
    dz(i,:) = circshift(dz(i-1,:),1);    
end

d = [dx;dy;dz]; %3D derivative operator
clearvars dx dy dz

%%
% f = Complex image
%
tic;
fOld = zeros(nVoxel,1);
fNew = ones(nVoxel,1); % To keep high norm
% a = (3*nx*ny*nz - nx*ny - ny*nz - nz*nx);

res = @(x,xOld)norm(x-xOld)/norm(xOld); % There will be some problem if norm(xOld) is 0
iter = 1;
resNorm = res(fNew,fOld);
while (iter <= maxIter) && (resNorm > del)
    L1 = diag((1./(((abs(fOld).^2)+eps).^(1-(k/2)))));
    L2 = diag((1./((((d*abs(fOld)).^2)+eps).^(1-(k/2)))));
    phi = diag(exp(-1j.*angle(fOld)));
    H = 2*(expConst'*expConst) + (k*const1^2).*L1 + ...
        (k*const2^2).*(phi'*((d'*L2)*(d*phi)));
    b = (1-stepSz).*(H*fOld) + (stepSz*2).*(expConst'*sParamCalib);
    
    % Solving Eqn. (15) by CGLS, set maxIter and tolerance later
    [fNew,flag,resNE,~] = cgls(H,b,[],[],maxIterCG,[],fOld); 
    fprintf('In iter = %d, resNorm = %f, flag = %d, resCGLS = %f\n',iter,resNorm,flag,resNE);
    resNorm = res(fNew,fOld);
    fOld = fNew;
    iter = iter + 1;
end

totIter = iter-1;
tComp = toc;

fprintf('Time: %3.2f s\n',tComp);

clearvars L1 L2 phi H b 
save([savePath,'\lsCetinImg',dataName,'.mat'],'fNew', '-v7.3')

%% Visualize image

imgBrightness = visImg(fNew,roomSize,voxelSize);

a = alphamap('rampup',256);
imgThresh = 120;
a(1:imgThresh)=0;
alphamap(a); 

title(['Cetin, k = ',num2str(k),' l1 = ',num2str(const1),', l2 = ',...
    num2str(const2), ', #iter = ',num2str(totIter),' imgTh = ',...
    num2str(imgThresh),' time = ',num2str(tComp),' s'],'FontSize',12)


end

