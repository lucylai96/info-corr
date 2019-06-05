function [CFI idx FI varCFI] = f_info(r_s1,r_s2,s1,s2)
% estimating f'
ds = (s2-s1)*(pi/180);
f_prime = (mean(r_s2)-mean(r_s1))./ds;

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

[E] = sort(eig(Sigma),'descend');
idx = max(find(E<0.001,1),0);
if isempty(idx)
    idx=0;
end

% estimating Fischer information
FI = (f_prime/Sigma)*f_prime'; %f_prime*inv(Sigma)*f_prime';

% bias-correction
T = size(r_s1,1);
N = size(r_s1,2);
CFI = FI*((2*T-N-3)/(2*T-2)) - ((2*N)/(T*ds^2));
varCFI =  (2*CFI^2)/(2*T-N-5) * (1 + 4*(2*T-3)/(T*CFI*ds^2) + 4*N*(2*T-3)/(T*CFI*ds^2)^2); % variance of FI bias corrected estimator

%CFI = CFI^-2;

% CFI = 1/((((1/CFI)*(pi/180)^2)));
% varCFI = 1/((((1/varCFI)*(pi/180)^2)));
% FI = 1/((((1/FI)*(pi/180)^2)));
end