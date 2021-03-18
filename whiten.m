function varargout = whiten(X,m)
% Return the whitend signals and the whitening matrix.
% 
% Format:
%      [Y,W] = whiten(X,m)
%      [Y,W] = whiten(X)
%          Y = whiten(X,m)
%          Y = whiten(X)
%
% Parameters
%       Y -- the whitened signals, each row representing a signal
%       W -- the whitening matrix, i.e., Y = W * X
%       X -- the input signals, each row representing a signal
%       m -- the dimention of the whitened signals,i.e.,Y is m-by-T dimention, where
%            T is the length of signals. Default value is the dimention of X.
%
% Author: Zhilin Zhang
% March 13, 2006



[n,T] = size(X); 
if nargin == 1
    m = n;        % default value
end

if m < n   % assumes white noise
 	[U,D] 	= eig((X*X')/T); 
	[puiss,k] = sort(diag(D));
 	ibl = sqrt(puiss(n-m+1:n)-mean(puiss(1:n-m)));  % subspace of signals and noise - that of noise
 	bl 	= ones(m,1) ./ ibl ;
 	W	= diag(bl)*U(1:n,k(n-m+1:n))';
 	IW 	= U(1:n,k(n-m+1:n))*diag(ibl);

    fprintf('whitening...Lost %g%% energy\n',100*sum(puiss(1:n-m))/(sum(diag(D))));
else       % assumes no noise
 	IW 	= sqrtm((X*X')/T);
 	W	= inv(IW);
end;

Y = W * X;

varargout{1} = Y;
varargout{2} = W;
