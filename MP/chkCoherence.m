% Checking mutual coherence of columns of input matrix A

function [condTrue,M,lowBoundM] = chkCoherence(A,k0)
% CHKMUTUALCOHERENCE
% -------------------------------------------------------------------------
% The columns of A are ai, i = 1,..., n. Number of rows = m.
% M = max|ai^H * aj|, i!=j, and ai^H*ai = 1, where H: Hermitian operator.
% The lower bound is M >= sqrt((n-m)/m(n-1))
% Sparsity condition is M < 1/(2k0 - 1)
% -------------------------------------------------------------------------
% Input: 
% A is the matrix for which the mutual coherence is defined. 
% Output:
% condTrue: 1 if the sparsity condition is satisfied, 0 ow.
% M: Structure, M.MC = Max Mutual Coherence, M.i and M.j are column numbers
% for which maximum coherence is achieved.
% -------------------------------------------------------------------------
% Pragya Sharma, ps847@cornell.edu
% 01/18/2019
% -------------------------------------------------------------------------

[m,n] = size(A);

for i = 1:n
    A(:,i) = A(:,i)/norm(A(:,i)); % Making  ai^H*ai = 1
end

M.MC = 0;
for colNum1 = 1:n
    for colNum2 = 1:n
        if colNum1 ~= colNum2
            corr = abs(A(:,colNum1)'*A(:,colNum2));
            if corr>M.MC
                M.MC = corr;
                M.i = colNum1;
                M.j = colNum2;
            end
        end
    end
end

lowBoundM = sqrt((n-m)/(m*(n-1)));
sparseCond = 1/(2*k0-1);
if M.MC < sparseCond
    condTrue = 1;
else
    condTrue = 0;
end
end