function [Q, R]= qr_h3(A,tol)
% f is a chebmatrix
% tic
[numRows,numCols] = size(A);
if ( numCols == 1 )
    % Trivial case: If A has only one column we simply scale it.
    R = sqrt(A'*A);
    if ( R ~= 0 )
        Q = A./R;
    else
        error('A=0');
    end
else
    % Legendre-Vandermonde matrix:
    L = legpoly(0:numCols-1, [0 numRows], 'norm', 1);
    %convert to chebmatrix
    Lm=cell(numRows,numCols);
    nele=numRows*numCols;
    for k=1:nele
        [ro,co]=ind2sub([numRows,numCols],k);
        Lm{k}=newDomain(restrict(L(:,co),[ro-1,ro]),[0,1]);        
    end
    Lm=chebmatrix(Lm);
    fprintf('Auxiliary orthonormal chebmatrix constructed.');
    % Call abstract QR:
    if ( nargin ==1 )
        tol = eps;
    end
    [Q, R] = absQR(A, Lm, @innerProdCbm, @norm,tol);
end
% toc
end

function out=innerProdCbm(A,B)
out=sum(sum(cellfun(@(a,b) sum(a.*b),A,B)));
end