function info_sampling
% purpose simulating and calculating F info

setPretty


%% verification that FI + bias correction works
% verifying that sampling the F info from a generated dataset with
% statistics

load('par.mat')

T = 800;

k = 1;
skp = 20;
%% do by cutting off TRIALS
%for dim = 252:4:800 %for these dimensionalities
for dim = 10:skp:800 %for these dimensionalities
    
    T = dim;
    d_th = 15;
   % if T<(N+5)/2
   %     N = 2*T-6;
   % end
    
    for i = 1:50
         N=500;
        rp1 = randperm(800);
        rp2 = randperm(500);
        if T>(N+5)/2 && T>=N % invertible and N<T (ideal
        R_S1 = r_s1(rp1(1:dim), rp2(1:N));
        R_S2 = r_s2(rp1(1:dim), rp2(1:N));
        [CFI(k,i), CFI2(k,i), FI(k,i), FI2(k,i), BC(k,i)] = f_info(R_S1,R_S2,d_th);
        
        elseif T>(N+5)/2 && T<N % invertible and N>T
         N = T;
         R_S1 = r_s1(rp1(1:dim), rp2(1:N));
         R_S2 = r_s2(rp1(1:dim), rp2(1:N));
        
         [CFI(k,i), CFI2(k,i), FI(k,i), FI2(k,i), BC(k,i)] = f_info(R_S1,R_S2,d_th);
        elseif T<(N+5)/2 % non-invertible
        R_S1 = r_s1(rp1(1:dim), rp2(1:N));
        R_S2 = r_s2(rp1(1:dim), rp2(1:N));
        [CFI(k,i), CFI2(k,i), FI(k,i), FI2(k,i), BC(k,i)] = ft_info(R_S1,R_S2,d_th);
         
        end
    end
    k = k+1;
end
% CFI is the "true information" in that population (treating it like the
% ground truth)
% FI is the same thing but in rads
% BC is the bias-corrected estimate
% CFI2 is the bias-corrected information from samples
% FI2 is the same thing but in rads

% plot bias-uncorrected vs bias corrected
figure;
hold on
shadedErrorBar([10:skp:800]',mean(BC,2),std(BC,[],2),[],1)
%plot(252:4:800,mean(BC,2)) %
xlabel('# trials')
ylabel('bias-corrected FI (rads^{-2})')
line([0 800],[BC(end) BC(end)],'Color','r')
line([250 250],[-20 140],'Color','b')
line([500 500],[-20 140],'Color','b')
prettyplot

% plot bias-corrected information against the # of trials

% estimated vs true FI
figure;
%shadedErrorBar(mean(CFI(100:end,:),2),mean(CFI2(100:end,:),2),std(CFI2(100:end,:),[],2),[],1)
errorbar(mean(CFI(14:end,:),2),mean(CFI2(14:end,:),2),std(CFI2(14:end,:),[],2),std(CFI2(14:end,:),[],2),[],[],'ro')
errorbar(mean(CFI(1:13,:),2),mean(CFI2(1:13,:),2),std(CFI2(1:13,:),[],2),std(CFI2(1:13,:),[],2),[],[],'bo')
xlabel('true FI (rads^{-2})')
ylabel('bias-corrected sampled FI (rads^{-2})')
axis equal
dline
prettyplot

%plot(CFI,CFI2,'o') %plotting all predicted vs actual


% std of bias-corrected estimator vs # trials
figure;
plot([10:skp:800],std(BC,[],2))
xlabel('# trials')
ylabel('std of bias-corrected FI (rads^{-2})')
prettyplot

plot(CFI, CFI2) %CFI is
figure;
plot(FI, FI2)
%% do by cutting off NEURONS
k = 1;
T = 800;
for dim = 200:4:500 %for these dimensionalities
    
    N = dim;
    d_th = 15;
    
    for i = 1:30
        rp2 = randperm(500);
        
        R_S1 = r_s1(1:T,rp2(1:dim));
        R_S2 = r_s2(1:T,rp2(1:dim));
        
        [CFI_N(k,i), CFI2_N(k,i), FI_N(k,i), FI2_N(k,i), BC_N(k,i)] = f_info(R_S1,R_S2,d_th);
    end
    
    k = k+1;
end


% plot bias-uncorrected vs bias corrected
figure;
hold on
shadedErrorBar([200:4:500]',mean(BC_N,2),std(BC_N,[],2),[],1)
%plot(252:4:800,mean(BC,2)) %
xlabel('# neurons')
ylabel('bias-corrected FI (rads^{-2})')
line([200 500],[BC_N(end) BC_N(end)],'Color','r')
prettyplot

% estimated vs true FI
figure;
%shadedErrorBar(mean(CFI(100:end,:),2),mean(CFI2(100:end,:),2),std(CFI2(100:end,:),[],2),[],1)
errorbar(mean(CFI_N,2),mean(CFI2_N,2),std(CFI2_N,[],2),std(CFI2_N,[],2),[],[],'o')
xlabel('true FI (rads^{-2})')
ylabel('bias-corrected sampled FI (rads^{-2})')
axis equal
dline
prettyplot

% std of bias-corrected estimator vs # trials
figure
plot([200:4:500],std(CFI2_N,[],2))
xlabel('# neurons')
ylabel('std of bias-corrected FI (rads^{-2})')
prettyplot

% std of bias-corrected estimator vs # trials
figure
plot([200:4:500],std(CFI2_N,[],2))
xlabel('# neurons')
ylabel('std of bias-corrected FI (rads^{-2})')
prettyplot

%% verification that FI + bias correction (on truncation) works
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


function [CFI, CFI2, FI, FI2, BC] = ft_info(r_s1,r_s2,d_th)
%estimating from a truncated matrix
T = size(r_s1,1);
% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./(d_th);

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FI using original parameters
%FI = (f_prime/Sigma)*f_prime'; %f_prime*inv(Sigma)*f_prime';
%CFI = 1/((((1/FI)*(pi/180)^2))); % go to degrees ^2


% low-d space
[E] = sort(eig(Sigma),'descend'); % singular decomposition
idx = find(E<0.001,1)-5;
%idx = (size(r_s1,1) + size(r_s2,1))-5; %cut off eig at K=2T-5
idx = T;
K = idx;

[EVect, EVal] = eig(Sigma);
V = flip(EVect,2);
A = V(:,1:idx)'; % low-d subspace

g_prime = (A*f_prime')';
V_Sigma = A*Sigma*A';

% FI using original parameters
FI = (g_prime/V_Sigma)*g_prime'; %f_prime*inv(Sigma)*f_prime';
BC = FI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(d_th)^2));
BC = 1/((((1/BC)*(pi/180)^2))); % go to degrees ^2
CFI = 1/((((1/FI)*(pi/180)^2))); % go to degrees ^2


% FI using unbiased sample parameters
V_Sig = (2.*V_Sigma)/(T*d_th^2);
dx_dth = mvnrnd(g_prime,V_Sig); % 1xN
S = wishrnd(V_Sigma/(2*(T-1)),(2*(T-1))); %Wishart
FI2 = (dx_dth/S)*dx_dth'; %f_prime*inv(Sigma)*f_prime';

% bias-correction
CFI2 = FI2*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(d_th)^2));
CFI2 = 1/((((1/CFI2)*(pi/180)^2)));
 
end