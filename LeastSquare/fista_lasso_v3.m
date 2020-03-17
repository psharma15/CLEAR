% This function implements Fast Iterative Shrinkage-Thresholding Algorithm
% (FISTA) for solving linear least squares (LS) problem with l1 (lasso)
% regularization. 
% Optimizes the problem of form: min_X 0.5*||Ax - b||^2 + gamma*||x||_1
% -------------------------------------------------------------------------
% Beck, Amir, and Marc Teboulle. "A fast iterative shrinkage-thresholding 
% algorithm for linear inverse problems." SIAM journal on imaging sciences 
% 2.1 (2009): 183-202.
% -------------------------------------------------------------------------
% Created on: 13 Nov 2018
% Last modified on: 15 Nov 2018
% Pragya Sharma, ps847@cornell.edu
% -------------------------------------------------------------------------
% This is using wavelets with objective function like
% min_X {||TX-b||+gamma*||X||_1}
% Where X: wavelet transform of x, T = A.Phi
% NOTE: FISTA cannot be used to solve with X only in the L1 norm term as 
% the proximal operator part derivation is not valid. It will change.

function [x,iter,relNorm,L] = fista_lasso_v3(A,b,opts)
%% Checking input values and defining default constants
if size(b,2) ~= 1
    fprintf('b is not a column vector. Please define correctly...\n');
    return;
end
if ~isfield(opts,'xinit')
    %opts.xinit = single(zeros(size(A,2),1)+(1e-8 + 1j*1e-8)); % Making a column vector
    opts.xinit = (zeros(size(A,2),1)+(1e-8 + 1j*1e-8));
elseif (length(opts.xinit) ~= size(A,2)) || (size(opts.xinit,2) > 1)
        fprintf('The size of xinit is wrong. Re-defining \n');
        opts.xinit = zeros(size(A,2),1)+(1e-8 + 1j*1e-8);
end
if ~isfield(opts,'maxIter') % Maximum number of iterations
    opts.maxIter = 10;
end
if ~isfield(opts,'tol') % Tolerance value
    opts.tol = (1e-4);
end
if ~isfield(opts,'eta') % Multiplying factor for Lipschitz constant eta^i_k
    opts.eta = 1.25;
end
if ~isfield(opts,'L0') % Initial Lipschitz constant for backtracking
    opts.L0 = 1; % max(eig(AHA)) ~ 5.6e4
end
if ~isfield(opts,'gamma') % Regularization parameter for ||x||_1
    opts.gamma = (1e-2);
end
if ~isfield(opts,'backtracking') % Backtracking on by default
    opts.backtracking = 1; 
end
if ~isfield(opts,'ifWavelet')
    opts.ifWavelet = 0;
end
if ~isfield(opts,'verbose') % Display 
    opts.verbose = 0; 
end

if opts.ifWavelet == 1
    fprintf('********************* FISTA Wavelet Version ***********************\n');
    fprintf('Update Lipscitz constant if needed.\n');
    if ~isfield(opts,'wavName')
        opts.wavName = 'Daubechies';
    end
    if ~isfield(opts,'par')
        opts.par = 4;
    end
    if ~isfield(opts,'wavL')
        opts.wavL = 1;
    end
    fprintf(['Using Wavelet: ',opts.wavName,' Par: ',num2str(opts.par),'\n']);
    % Getting the filter
    qmf = MakeONFilter(opts.wavName,opts.par);
    figure;plot(qmf)
    nx = opts.nVoxel(1);
    ny = opts.nVoxel(2);
    nz = opts.nVoxel(3);
    n = nx*ny*nz;
    % Assuming nx,ny,nz are powers of 2
    if mod(log2(n),1)~=0 
        warning('Input voxel numbers in powers of 2')
        return;
    end
    if opts.calcT 
        T = zeros(size(A));
        for ii = 1:n
            ek = zeros(n,1);
            ek(ii) = 1;
%             psi = IWT_PO(ek,opts.wavL,qmf); if ii == 1, fprintf('\nUsing DCT\n\n'); end
            psi = idct(ek); if ii == 1, fprintf('\nUsing DCT\n\n'); end
            T(:,ii) = A*psi;
        end
        save(['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\'...
            'Algorithms\Wavelet\','T_',opts.wavName(1:4),num2str(opts.par), num2str(opts.wavL)],'T');
    else
        load(['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\'...
            'Algorithms\Wavelet\','T_',opts.wavName(1:4),num2str(opts.par), num2str(opts.wavL)],'T');
    end
    A = T;
    
%     opts.xinit = pinv(A)*b;
end

%% Some useful functions

function res = gradfx(A,b,x)
    % Gradient of f(x) = AH(Ax-b)
%     (size(A))
%     size(x)
%     size(b)
    res = A'*(A*x - b);
end

function res = proxOpL1(x,gamma,L)
    % Proximal operator of L1-norm term g(x) with complex x. 
    % Preserving the phase information
    
    % ************Update this for wavelet
    res = max(abs(x) - gamma/L, 0).*exp(1j*angle(x));
end

function res = calc_fx(A,b,x)
    % Calculating f(x) = 0.5 ||Ax-b||^2
    res = 0.5*norm(A*x - b)^2;
end

function res = calc_gx(x,gamma)
    % Calculating g(x) =  gamma*||x||_1
    
    if numel(gamma) == 1 
        res = gamma*sum(abs(x));
    elseif (size(gamma,2) == 1) && (numel(gamma) == numel(x))
        res = sum(gamma(:).*abs(x(:)));
    end
end

function res = calcFx(A,b,x,gamma)
    % Calculating F(X) = f(x) + g(x)
    res = calc_fx(A,b,x) + calc_gx(x,gamma);
end

function res = calcQL(A,b,x,y,L,gamma)
    % Q_L(x,y) = f(y) + <x-y,grad(f(y))> + (L/2)||x-y||^2 + g(x)
    res = calc_fx(A,b,y)+((x-y)'*gradfx(A,b,y))+...
        ((L/2)*norm(x-y)^2)+calc_gx(x,gamma);
end

%% Carrying out the iterations

iter = 0; % k
xOld = opts.xinit; % x_(k-1)
yOld = opts.xinit; % y_k
tOld = 1; % t_k
L = opts.L0; % L_k, here i_k = 0 which is smallest nonnegative integer.

% AHA = A'*A; % Check if this needs to be saved or recalculated
% AHb = A'*b;
% x = proxOpL1(yOld,opts.gamma,L);

fprintf('Carrying out FISTA iterations. \n');
while iter < opts.maxIter
    iter = iter + 1;
    while opts.backtracking
        % Finding i_k
        pL = proxOpL1(yOld - (1/L)*gradfx(A,b,yOld),opts.gamma,L); % pL(y_k)
        FpL = calcFx(A,b,pL,opts.gamma); % F(pL(y_k))
        QLpL = calcQL(A,b,pL,yOld,L,opts.gamma);
        
        if opts.verbose
            fprintf('iter = %d, L = %g, FpL = %g, QLpL = %g\n',...
                iter,L,FpL,abs(QLpL));
        end   
        
        if abs(FpL) <= abs(QLpL) % Originally ONLY LESS THAN CONDITION: I'm testing equals to
            % When they are equal, do not do backtracking, it diverges.
            % Check this condition for the complex case
%             fprintf('Stopping: abs(FpL)<abs(QLpL)\n');
            break
        end
        
        L = opts.eta * L; % L_k
    end
    
    x = proxOpL1((yOld - (1/L)*gradfx(A,b,yOld)),opts.gamma,L); % x_k
    costFxOld = calcFx(A,b,xOld,opts.gamma);
    costFx = calcFx(A,b,x,opts.gamma);
    if costFx > costFxOld
        x = xOld;
        fprintf("Solution is diverging, stopping ...\n");
        break;
    end
    relNorm = norm(x - xOld)/length(x); % This needs to be updated. Not correct.
    relNorm2 = norm(x-xOld)/norm(x);
    if opts.verbose
        fprintf('iter = %d, L = %g, costFx = %g, relNorm = %g\n',iter,...
            L,costFx,relNorm2);
    end   
    if relNorm2 < opts.tol
        break;
    end
    
    t = (1 + sqrt(1 + 4*tOld^2))/2; % t_(k+1)
    y = x + ((tOld - 1)/t)*(x - xOld); % y_(k+1)    
    yOld = y;  
    tOld = t;
    xOld = x;
    
    opts.gamma = opts.gamma*opts.gammaScale;
end

if opts.ifWavelet 
    fprintf('Converting x back from wavelet.\n');
    X = x;
    x = zeros(n,1);
    for ii = 1:n
        ek = zeros(n,1);
        ek(ii) = 1;
%         psi = IWT_PO(ek,opts.wavL,qmf); % IWT_PO
        psi = idct(ek);
        x = x + psi*X(ii);
    end
end
end

