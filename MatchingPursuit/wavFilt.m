% Wavelet transformation for Orthogonal Matching Pursuit (OMP).

% Input:
function [out] = wavFilt(opts)
% wavFilt generated PHI = A*PSI, where PSI here is the matrix with wavelet
% filter coefficients. Works only for 1 level, 1 D data.
% Input: 
% opts: Options, default values are defined below
% Output: 
% Either of the following, depends on opts.reconstruct
% phi = A*psi, where psi is the wavelet filter matrix
% x = psi*Image, where Image is the wavelet transformed result of x.

%%
if ~isfield(opts,'imSz')
    opts.imSz = [0.1, 0.8; 0.1, 0.8; 0, 0.3]; % meters
end
if ~isfield(opts,'voxSz')
    opts.voxSz = [0.015;0.015;0.015]; % meters
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
    opts.wName = 'db4'; % db2 wavelets by default
end
if ~isfield(opts,'reconstruct')
    % By default this code generates Phi, keep this value 1 to reconstruct 
    % x, ie, inverse wavelet transform of the obtained result.
    opts.reconstruct = 0;
end

%% 
[m,n] = size(opts.A);
nNextPow = 2^nextpow2(n);
opts.A = [opts.A, zeros(m,nNextPow-n)];

[Lo,Ho] = wfilters(opts.wName,'r');
Lo = Lo(:); Ho = Ho(:); % Converting to column format
lFilt = length(Lo); % Length of filter

% Not calculating psi directly: [n x n]
phi = zeros(m,n);

if opts.reconstruct == 0
    fprintf('Generating Phi...\n');
    tic
    for colNum = 1:2:n-1
        
        
        if colNum <= n-2
            psi_1 = zeros(nNextPow,1);
            psi_2 = zeros(nNextPow,1);
            psi_1(colNum:(colNum+lFilt-1)) = Lo;
            psi_2(colNum:(colNum+lFilt-1)) = Ho;
        else
            % Handling last 2 columns
            psi_1(colNum:(colNum+(lFilt/2)-1)) = Lo(1:(lFilt/2));
            psi_1(1:(lFilt/2)) = Lo(((lFilt/2)+1):end);
            psi_2(colNum:(colNum+(lFilt/2)-1)) = Ho(1:(lFilt/2));
            psi_2(1:(lFilt/2)) = Ho(((lFilt/2)+1):end);
        end
        phi(:,colNum) = opts.A*psi_1;
        phi(:,colNum+1) = opts.A*psi_2;
        if (colNum+1) == 5000
            timed5000 = toc;
        end
        if mod(colNum+1,5000) == 0
            fprintf('Column Number %d out of %d\n',colNum,n);
            timeRemain = ((ceil(n/5000) - ((colNum+1)/5000))*timed5000)/60;
            fprintf('Approx. remaining time = %f min\n',timeRemain);
        end

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
        if colNum <= n-2
            psi_1 = zeros(nNextPow,1);
            psi_2 = zeros(nNextPow,1);
            psi_1(colNum:(colNum+lFilt-1)) = Lo;
            psi_2(colNum:(colNum+lFilt-1)) = Ho;
        else
            % Handling last 2 columns
            psi_1(colNum:(colNum+(lFilt/2)-1)) = Lo(1:(lFilt/2));
            psi_1(1:(lFilt/2)) = Lo(((lFilt/2)+1):end);
            psi_2(colNum:(colNum+(lFilt/2)-1)) = Ho(1:(lFilt/2));
            psi_2(1:(lFilt/2)) = Ho(((lFilt/2)+1):end);
        end
        x = x + psi_1.*opts.imgWav(colNum) + psi_2.*opts.imgWav(colNum+1);
    end
    x = x(1:n);
    out = x;
end
clear psi_1 psi_2
        
    
