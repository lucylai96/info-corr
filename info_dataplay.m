function info_dataplay
% purpose: look at harvey lab data 

%deResp: matrix of trials x neurons structure of deconvolved spike sum proxies from calcium traces
%mvSpd: vector of movement speeds (abstract number, larger = faster), one per trial
%visCon: vector of stimulus constrasts, one per trial
%visOri: vector of stimulus drift directions, one per trial
addpath('/Users/lucy/Google Drive/Harvard/Rotations/drugo-lab/')


load('j1_171018.mat')

%% take only data from the highest contrasts for now

highConResp = deResp(visCon==max(visCon),:);
highConMov = mvSpd(visCon==max(visCon));
highVisOri = visOri(visCon==max(visCon));

%% separate out data by orientations

highConResp0 = highConResp(highVisOri==0,:);
highConResp0 = trunc(highConResp0,800);

highConResp15 = highConResp(highVisOri==15,:);
highConResp15 = trunc(highConResp15,800);

highConResp30 = highConResp(highVisOri==30,:); 
highConResp30 = trunc(highConResp30,800);

highConResp45 = highConResp(highVisOri==45,:); 
highConResp45 = trunc(highConResp45,800);

highConResp60 = highConResp(highVisOri==60,:); 
highConResp60 = trunc(highConResp60,800);

highConResp75 = highConResp(highVisOri==75,:); 
highConResp75 = trunc(highConResp75,800);

highConResp90 = highConResp(highVisOri==90,:);
highConResp90 = trunc(highConResp90,800);

% ~500 neurons and ~800 trials

%% F information
% take two diff orientations 15 deg apart

f_info(highConResp0, highConResp15, 0, 15)
f_info(highConResp15, highConResp30, 15, 30)
f_info(highConResp30, highConResp45, 30, 45)
f_info(highConResp45, highConResp60, 45, 60)
f_info(highConResp60, highConResp75, 60, 75)

f_info(highConResp15, highConResp60, 15, 60)

%[U,S,V] = svd(Sigma);


%% truncated data (how many trials can you take out before the F info differs)
% has to be less trials for both conditions 

t_highConResp0 = trunc(highConResp0, 200);
t_highConResp15 = trunc(highConResp15, 200);
t_highConResp30 = trunc(highConResp30, 200);

ft_info(t_highConResp0, t_highConResp15, 0, 15);

%% estimated FI from full information vs estimated FI from estimated information 


end

function CFI = f_info(r_s1,r_s2,s1,s2)
% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./(s2-s1);

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

% estimating Fischer information
FI = f_prime*inv(Sigma)*f_prime';

% bias-correction
T = size(r_s1,1);
N = size(r_s1,2);
CFI = FI*((2*T-N-3)/(2*T-2)) - ((2*N)/(T*(s2-s1)^2));
CFI = 1/((((1/CFI)*(pi/180)))^2);
end

function FI = ft_info(r_s1,r_s2,s1,s2)
%estimating from a truncated matrix

% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./(s2-s1);

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

% singular decomposition
[E] = sort(eig(Sigma),'descend');
idx = find(E<0.001,1);
%idx = size(r_s1,1) + size(r_s2,1);
[U,S,V] = svd(Sigma);

A = V(:,1:idx)'; % low-d subspace

%rt_s1 = (A*r_s1')';
%rt_s2 = (A*r_s2')';
%g_prime = (mean(rt_s2)-mean(rt_s1))./(s2-s1);

g_prime = (A*f_prime')';
V_Sigma = A*Sigma*A';

% estimating Fischer information
TFI = g_prime*inv(V_Sigma)*g_prime';

% bias-correction (derive this!)
K = size(g_prime,2);
T = size(r_s1,1);
CTFI = TFI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(s2-s1)^2));
CTFI = 1/(((1/CTFI)*(pi/180))^2); %from inv degrees to radians

end

function t_matrix = trunc(matrix, num)
% truncate the data matrix
% input: matrix is the matrix you want to truncate, num is how many trials you want to keep
% output: t_matrix output matrix

t_matrix  = matrix(1:num,:);

end
