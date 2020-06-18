% Wavelet transformation for Orthogonal Matching Pursuit (OMP).

% Input:
function [out] = wavFilt3d(opts)
% wavFilt generated PHI = A*PSI, where PSI here is the matrix with wavelet
% filter coefficients. Works only for 3D data at multilevel
% Input: 
% opts: Options, default values are defined below
% Output: 
% Either of the following, depends on opts.reconstruct
% phi = A*psi, where psi is the wavelet filter matrix
% x = psi*Image, where Image is the wavelet transformed result of x.

%%
if ~isfield(opts,'imSz')
    opts.imSz = [0.1, 0.72; 0.1, 0.72; 0, 0.3]; % meters
end
if ~isfield(opts,'voxSz')
    opts.voxSz = [0.02;0.02;0.02]; % meters
end
if ~isfield(opts,'freq')
    opts.freq = [1.7006;1.7591;1.8185;1.8779;1.9364;1.9958;2.0543;2.1137;2.1722;2.2316;2.2901]*1e9;
end
if ~isfield(opts,'posRxTx')
    opts.posRxTx = 1;
end
if ~isfield(opts,'A')
    [~,~,rxPosition,tagPosition] = posRxTx(opts.posRxTx);
    opts.A = genA(opts.imSz,opts.voxSz,tagPosition,rxPosition,opts.freq);
    
end
if ~isfield(opts,'wName')
    opts.wName = 'db8'; % db8 wavelets by default
end
if ~isfield(opts,'wLevel')
    opts.wLevel = 2; % Level 2 by default
end
if ~isfield(opts,'reconstruct')
    % By default this code generates Phi, keep this value 1 to reconstruct 
    % x, ie, inverse wavelet transform of the obtained result.
    opts.reconstruct = 0;
end

%% 
[m,n] = size(opts.A);
[~,~,~,nVoxel] = genXYZ(opts.imSz,opts.voxSz);
opts.A = reshape(opts.A,[m,nVoxel(1),nVoxel(2),nVoxel(3)]);
nPowOf2 = nextpow2(nVoxel);

for i = 1:3
    if 2^nPowOf2(i) ~= nVoxel(i)
        fprintf('Voxel number should be in powers of 2.\n');
        return;
    end
end

[Lo,Ho] = wfilters(opts.wName,'r');
Lo = Lo(:)'; Ho = Ho(:)'; % Ensure filters are in row format
lFilt = length(Lo); % Length of filter
fprintf(' Using %s wavelet at level %d.\n',opts.wName,opts.wLevel);

% Not calculating psi directly: [n x n]
phi = zeros([m,nVoxel(1),nVoxel(2),nVoxel(3)]);

if opts.reconstruct == 0
    fprintf('Generating Phi...\n');
    for i = 1:m
        for k = 1:opts.wLevel
            wfilt = 
    % Write Phi generation here
    
    phi(i,:) = ;
    end
    out = phi;
else
    clear opts.A
    x = zeros(nNextPow,1);
    if ~isfield(opts,'imgWav')
        fprintf('Input image in wavelet domain for reconstruction.\n');
        return
    else
        opts.imgWav = [opts.imgWav(:); zeros(nNextPow-n,1)];
        fprintf('Reconstructing image from wavelet domain.\n');
    end
    for colNum = 1:2:n-1
        if mod(colNum+1,5000) == 0
            fprintf('Column Number %d out of %d\n',colNum, n);
        end
        
        % Write Inverse transform here
        
        
        x = x + psi_1.*opts.imgWav(colNum) + psi_2.*opts.imgWav(colNum+1);
    end
    x = x(1:n);
    out = x;
end
clear psi_1 psi_2
        
    
