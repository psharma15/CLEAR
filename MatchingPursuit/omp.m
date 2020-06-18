% Orthogonal Matching Pursuit algorithm (OMP), solves b = Ax, where b is
% the available output data, A is the dictionary and x is the unknown input
% value
% Pragya Sharma, ps847@cornell.edu
% Last update on: 1/4/19

function [x,idxSet,objIdx,iter,normRes] = omp(A,b,opts)
%% OMP 
% Input: 
% A: This is the dictionary with columns as basis vectors, size M x N.
% b: This is the output data we have, size M.
% opts: Following options can be defined
% opts.sc: Stopping criteria based on sparsity in percentage.
% Output:
% x: This is the unknown input vector values for non-zero positions, size<N.
% idxSet: This is the index values that are non-sparse and takes the values
% given by x.
% iter: Number of OMP iterations, related to stopping criteria.
% res: This is the residual vector.
% -------------------------------------------------------------------------
%% 
% Initializations
b = b(:); % b needs to be a column vector
[~,n] = size(A);
res = b;
x = zeros(n,1);
iter = 0;
tol = 0;
idxSet = zeros(n,1); % Binary vector with 1s corresponding to selected atom 
atomNumLib = []; %Library of atom numbers selected
if ~isfield(opts,'scPercent')
    opts.scPercent = 0; % Stopping sparcity criteria is not % by default
end
if opts.scPercent == 1
    opts.sc = ceil(opts.sc*n);
end
    
% Normalizing each column of A if option is selected
if opts.norm == 1
    fprintf('Normalizing columns of A for OMP \n');
    for i = 1:n
        A(:,i) = A(:,i)/norm(A(:,i));
    end
end

normRes = zeros(opts.sc+1,1);
normRes(1) = norm(res);

% Stopping rule could be changed, currently 'tol' is basically just number
% of iterations as the stopping condition is maxIter.
while  tol < opts.sc
    innerProd = abs(A'*res);
    corrMax = 0;
    for i = 1:n
        j = 1;
        idxColRemove = 0; % Remove this index from consideration if 1
        while j <= length(atomNumLib)
            if i == atomNumLib(j)
                idxColRemove = 1;
            end
            j = j + 1;
        end
        if (idxColRemove == 0) && (innerProd(i) > corrMax)
            corrMax = innerProd(i);
            atomMax = i;
        end
    end
    
    idxSet(atomMax) = 1;
    atomNumLib = [atomNumLib;atomMax]; % This is the set of indices of important columns of A
    atomNumLib = sort(atomNumLib);
    x = ((A(:,atomNumLib)'*A(:,atomNumLib)+opts.regConst*eye(length(atomNumLib)))\A(:,atomNumLib)')*b; % LS solution of x 
    res = b - (A(:,atomNumLib)*x); % Check res-Anew*x or b-Anew*x
    iter = iter + 1;
    tol = iter; % Update stopping rule
    normRes(tol+1,1) = norm(res);
    fprintf('Residual in Iter %d = %.2f.\n',iter,normRes(tol,1));
    
end

objIdx = find(idxSet~=0);
% objIdx
% x

tk = zeros(opts.sc,1);
for i = 2:opts.sc+1
    % The indices are slightly different, as normRes(0) needs to be defined
    % at index 1.
    tk(i-1) = normRes(i)/normRes(i-1); 
end

if opts.RRfig == 1
    figure('Position',[700,500,500,250])
    plot(1:opts.sc,tk(1:end));
    xlabel('k');ylabel('RR(k)');
    hold on
    
    plot(opts.kVal,tk(opts.kVal),'*','markersize',12);
%     ylim([0.9,1.01])
    title(['RR curve, k = ',num2str(opts.kVal)]);
    
    figure('Position',[700,500,500,250])
    plot(objIdx,abs(x));
    title(['Stopping k = ',num2str(opts.sc)]);
    xlabel('objIdx','FontSize',12)
    ylabel('abs(x)','FontSize',12)
    title('abs(x) for each library element');
end
end
    

