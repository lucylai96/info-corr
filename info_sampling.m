function info_sampling
% purpose simulating and calculating F info

setPretty

%% verification that FI + bias correction works
% verifying that sampling the F info from a generated dataset with
% statistics 

load('par.mat')

T = 800;
N = 500;
k = 1;
%% do by cutting off TRIALS
%for dim = 252:4:800 %for these dimensionalities
for dim = 252:4:800 %for these dimensionalities

T = dim;
d_th = 15;

R_S1 = r_s1(1:T, :);
R_S2 = r_s2(1:T, :);

[CFI(k), CFI2(k), FI(k), FI2(k), BC(k)] = f_info(R_S1,R_S2,d_th);

k = k+1;
end
% CFI is the "true information" in that population (treating it like the
% ground truth)
% FI is the same thing but in rads
% BC is the bias-corrected estimate 
% CFI2 is the bias-corrected information from samples 
% FI2 is the same thing but in rads

% plot true

% plot bias-corrected information against the # of trials 
plot(252:4:800,CFI2) %
line([200 800],[CFI2(end) CFI2(end)],'Color','r')
%line([(N+2)/2 (N+2)/2],[min(CFI) max(CFI)],'Color','k')
xlabel('# trials')
ylabel('bias-corrected')


plot(CFI, CFI2) %CFI is 
figure;
plot(FI, FI2)
%% do by cutting off NEURONS

%% verification that FI + bias correction works
% verifying that sampling the F info from a generated dataset with
% statistics 



end 

function [CFI, CFI2, FI, FI2, BC] = f_info(r_s1,r_s2,d_th)
T = size(r_s1,1);
N = size(r_s1,2);
% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./(d_th);

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FI using original parameters
FI = (f_prime/Sigma)*f_prime'; %f_prime*inv(Sigma)*f_prime';
BC = FI*((2*T-N-3)/(2*T-2)) - ((2*N)/(T*(d_th)^2));
BC = 1/((((1/BC)*(pi/180)^2))); % go to degrees ^2
CFI = 1/((((1/FI)*(pi/180)^2))); % go to degrees ^2

% FI using unbiased sample parameters
Sig = (2.*Sigma)/(T*d_th^2);
du_dth = mvnrnd(f_prime,Sig); % 1xN
S = wishrnd(Sigma/(2*(T-1)),(2*(T-1))); %Wishart
FI2 = (du_dth/S)*du_dth'; %f_prime*inv(Sigma)*f_prime';

% bias-correction
CFI2 = FI2*((2*T-N-3)/(2*T-2)) - ((2*N)/(T*(d_th)^2));
CFI2 = 1/((((1/CFI2)*(pi/180)^2)));

end