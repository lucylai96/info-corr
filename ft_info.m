function [CTFI, idx, TFI, varCTFI] = ft_info(r_s1,r_s2,s1,s2)
% output [CorrectTruncated_FisherInfo, idx, naiveTFI, varTFI ]
ds = (s2-s1)*(pi/180);
%estimating from a truncated matrix
T = size(r_s1,1);
% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./ds;

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

% singular decomposition
[E] = sort(eig(Sigma),'descend');
idx = find(E<0.001,1);
%idx = (size(r_s1,1) + size(r_s2,1))-5; %cut off eig at K=2T-5
%idx = T;
%[U,S,V] = svd(Sigma);
[EVect, EVal] = eig(Sigma);
V = flip(EVect,2);
%var_vec = flipud(diag(EVal));
%pct = cumsum(var_vec) / sum(diag(EVal));% Calculate the cumulative percentage of the variances.
% Find the 99% variance
%idx = find(pct >= 0.99, 1);
%figure(1);hold on;plot(E)

%A = V(:,1:ceil(idx-5))'; % low-d subspace
A = V(:,1:ceil(T))'; % low-d subspace

g_prime = (A*f_prime')';
V_Sigma = A*Sigma*A';

% estimating Fischer information
TFI = (g_prime/V_Sigma)*g_prime'; %naive estimate

% bias-correction (derive this!)
K = size(g_prime,2);

CTFI = TFI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*ds^2));
varCTFI =  (2*CTFI^2)/(2*T-K-5) * (1 + 4*(2*T-3)/(T*CTFI*ds^2) + 4*K*(2*T-3)/(T*CTFI*ds^2)^2); % variance of FI bias corrected estimator

%from inv degrees to radians
%CTFI = 1/(((1/CTFI)*(pi/180)^2)); 
%varCTFI= 1/(((1/varCTFI)*(pi/180)^2)); 
%TFI = 1/(((1/TFI)*(pi/180)^2)); 

end