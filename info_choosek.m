function info_choosek
%choosek(0)
choosek(4)
end
function choosek(ck)
% choose the numbder of eigenvectors you want to keep
load ModelData
FI_TRUE
% resp = 50 neurons, 2 orientations, 100000 trials
%resp = resp(:,:,1:300); %reconfigure into trials by neurons by conditions
new_resp(:,:,1) = squeeze(resp(:,1,:))';
new_resp(:,:,2) = squeeze(resp(:,2,:))';
highConRespOri = new_resp;

%ck = 0; just tak away neurons
%ck = 1;% just take away neurons all the way down after N<T
%ck = 2;% choose K=T
%ck = 3;% choose K = max 99% of the variance
%ck = 4;% choose K based on the one that max(FI)

load('j1_171018_estFI.mat')
highConRespOri = new_resp;
ORI = [-7 0]*pi/180; %rad
DORI = [1 2];
or_corr=[1 2];
ds = diff(ORI(or_corr));

skp= 5; % step size

b = 50; %bootstraps
sz = length(10:skp:300); % size of the calculated vectors
bs_TFI_N = zeros(sz,1);
err_TFI_N = zeros(sz,1);
bs_naiveTFI_N = zeros(sz,1);
err_naiveTFI_N = zeros(sz,1);
var_TFI_N = zeros(sz,1);

bs_TFI_T = zeros(sz,1);
err_TFI_T = zeros(sz,1);
bs_naiveTFI_T = zeros(sz,1);
err_naiveTFI_T = zeros(sz,1);
var_TFI_T = zeros(sz,1);

TFI_N = zeros(b,1);
naiveTFI_N = zeros(b,1);
varTFI_N = zeros(b,1);
TFI_T = zeros(b,1);
naiveTFI_T = zeros(b,1);
varTFI_T = zeros(b,1);
if ck ==0
    k = 1;
    %% take away only neurons
    
    for dim = 5:skp:50 %for these dimensionalities
        %% do by cutting off NEURONS
        t = 100000; %use all the trials
        n = dim; %only use "dim" # of neurons
        for b  = 1:50 %bootstrap 100 times
            rp1=randperm(100000);
            rp2=randperm(50);
            
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
            t_highConRespOri_bsN = zeros(t,n,length(ori));
            for i  = 1:2
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsN(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i);
                
            end
            
            for i  = 1 %comparisons are # of orientations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                % if  t>(n+2)/2 %if you have more trials than neurons
                [TFI_N(b,i) idx(b,i) naiveTFI_N(b,i) varTFI_N(b,i)] = f_inforad(t_highConRespOri_bsN(:,:,i), t_highConRespOri_bsN(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                % end
            end
        end
        %each row is a dimension (# neurons)
        bs_TFI_N(k,:) = mean(TFI_N(:,1)); % here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_N(k,:) = std(TFI_N(:,1));
        var_TFI_N(k,:) = mean(varTFI_N(:,1)); %any variance on this reflects the variance of the 6 orientation difference conditions and bootstrapping
        
        idx_TFI_N(k,:) = mean(idx(:,1),1);
        
        bs_naiveTFI_N(k,:) = mean(naiveTFI_N(:,1));
        err_naiveTFI_N(k,:) = std(naiveTFI_N(:,1));
        k = k+1;
    end
    
    
    figure;hold on;
    shadedErrorBar([5:skp:dim]',bs_TFI_N(1:length([5:skp:dim])),err_TFI_N(1:length([5:skp:dim])),{'g','LineWidth',1},1)
    line([0 100],[FI_TRUE FI_TRUE],'Color','r')
    prettyplot
    
elseif ck ==1
    %% just take away neurons all the way down
    k = 1;
    for dim = 5:skp:100 %for these dimensionalities
        %% here only use 1 orientation difference
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:100 %bootstrap 100 times
            n = 50; % keep all neurons
            rp1=randperm(100000);
            rp2=randperm(50);
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),2);
            t_highConRespOri_bsT = zeros(t,n,2);
            for i  = 1:2
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i); %t x n
            end
            
            for i  = 1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 50;
                if t>(n+5)/2 && t>=n %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_inforad(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                elseif t<n
                    n = t;
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_inforad(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                    
                end
                %                 elseif t>(n+5)/2 && t<n
                %                     n = t;
                %
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
                %                 elseif t<(n+5)/2 %less trials than neurons
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
            end
            
        end
        
        %each row is a dimension (# trials)
        bs_TFI_T(k,:) = mean(TFI_T(:,1)); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_T(k,:) = std(TFI_T(:,1));
        var_TFI_T(k,:) = mean(varTFI_T(:,1)); %any variance on this reflects the variance of bootstrapping
        
        idx_TFI_T(k,:) = mean(idx,1);
        
        bs_naiveTFI_T(k,:) = mean(naiveTFI_T(:,1));
        err_naiveTFI_T(k,:) = std(naiveTFI_T(:,1));
        
        k = k+1;
    end% dimensions
    
    figure;hold on;
    shadedErrorBar([5:skp:dim]',bs_TFI_T(1:length([5:skp:dim])),err_TFI_T(1:length([5:skp:dim])),[],1)
    line([0 100],[FI_TRUE FI_TRUE],'Color','r')
    line([27.5 27.5],[-300 400],'Color','b')
    prettyplot
    
    
    
    
elseif ck ==2  %% choose K=T
    
    
    k = 1;
    for dim = 5:skp:100 %for these dimensionalities
        %% here only use 1 orientation difference
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:100 %bootstrap 100 times
            n = 50; % keep all neurons
            rp1=randperm(100000);
            rp2=randperm(50);
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),2);
            t_highConRespOri_bsT = zeros(t,n,2);
            for i  = 1:2
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i); %t x n
            end
            
            for i  = 1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 50;
                if t>(n+5)/2 %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_inforad(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                else
                 
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info1(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                    
                end
                %                 elseif t>(n+5)/2 && t<n
                %                     n = t;
                %
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
                %                 elseif t<(n+5)/2 %less trials than neurons
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
            end
            
        end
        
        %each row is a dimension (# trials)
        bs_TFI_T(k,:) = mean(TFI_T(:,1)); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_T(k,:) = std(TFI_T(:,1));
        var_TFI_T(k,:) = mean(varTFI_T(:,1)); %any variance on this reflects the variance of bootstrapping
        
        idx_TFI_T(k,:) = mean(idx,1);
        
        bs_naiveTFI_T(k,:) = mean(naiveTFI_T(:,1));
        err_naiveTFI_T(k,:) = std(naiveTFI_T(:,1));
        
        k = k+1;
    end% dimensions
    figure;hold on;
    shadedErrorBar([5:skp:dim]',bs_TFI_T(1:length([5:skp:dim])),err_TFI_T(1:length([5:skp:dim])),[],1)
    line([0 100],[FI_TRUE FI_TRUE],'Color','r')
    line([27.5 27.5],[-300 400],'Color','b')
 
    
    load('ck_neurons.mat')
    shadedErrorBar([5:skp:dim]',bs_TFI_N(1:length([5:skp:dim])),err_TFI_N(1:length([5:skp:dim])),{'g','LineWidth',1},1)
       prettyplot
    title('choose k = t')
    
elseif ck ==3 %% choose k = 99% variance
    
    k = 1;
    for dim = 5:skp:100 %for these dimensionalities
        %% here only use 1 orientation difference
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:100 %bootstrap 100 times
            n = 50; % keep all neurons
            rp1=randperm(100000);
            rp2=randperm(50);
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),2);
            t_highConRespOri_bsT = zeros(t,n,2);
            for i  = 1:2
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i); %t x n
            end
            
            for i  = 1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 50;
                if t>(n+5)/2 %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_inforad(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                else
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info2(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                    
                end
                %                 elseif t>(n+5)/2 && t<n
                %                     n = t;
                %
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
                %                 elseif t<(n+5)/2 %less trials than neurons
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
            end
            
        end
        
        %each row is a dimension (# trials)
        bs_TFI_T(k,:) = mean(TFI_T(:,1)); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_T(k,:) = std(TFI_T(:,1));
        var_TFI_T(k,:) = mean(varTFI_T(:,1)); %any variance on this reflects the variance of bootstrapping
        
        idx_TFI_T(k,:) = mean(idx,1);
        
        bs_naiveTFI_T(k,:) = mean(naiveTFI_T(:,1));
        err_naiveTFI_T(k,:) = std(naiveTFI_T(:,1));
        
        k = k+1;
    end% dimensions
    figure;hold on;
    shadedErrorBar([5:skp:dim]',bs_TFI_T(1:length([5:skp:dim])),err_TFI_T(1:length([5:skp:dim])),[],1)
    line([0 100],[FI_TRUE FI_TRUE],'Color','r')
    line([27.5 27.5],[-300 400],'Color','b')
    prettyplot
    
     load('ck_neurons.mat')
    shadedErrorBar([5:skp:dim]',bs_TFI_N(1:length([5:skp:dim])),err_TFI_N(1:length([5:skp:dim])),{'g','LineWidth',1},1)
       prettyplot
           title('choose k = 99% variance')
    
elseif ck ==4 %%  choose K based on the one that max(FI)
    
    
    k = 1;
    for dim = 5:skp:100 %for these dimensionalities
        
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:100 %bootstrap 100 times
            n = 50; % keep all neurons
            rp1=randperm(100000);
            rp2=randperm(50);
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),2);
            t_highConRespOri_bsT = zeros(t,n,2);
            for i  = 1:2
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i); %t x n
            end
            
            for i  = 1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 50;
                if t>(n+5)/2 %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_inforad(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                else
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info3(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                    
                end
                %                 elseif t>(n+5)/2 && t<n
                %                     n = t;
                %
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
                %                 elseif t<(n+5)/2 %less trials than neurons
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds);%this is the bias-corrected FI in units of rad^-2
                %
            end
            
        end
        
        %each row is a dimension (# trials)
        bs_TFI_T(k,:) = mean(TFI_T(:,1)); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_T(k,:) = std(TFI_T(:,1));
        var_TFI_T(k,:) = mean(varTFI_T(:,1)); %any variance on this reflects the variance of bootstrapping
        
        idx_TFI_T(k,:) = mean(idx,1);
        
        bs_naiveTFI_T(k,:) = mean(naiveTFI_T(:,1));
        err_naiveTFI_T(k,:) = std(naiveTFI_T(:,1));
        
        k = k+1;
    end% dimensions
    
    figure;hold on;
    shadedErrorBar([5:skp:dim]',bs_TFI_T(1:length([5:skp:dim])),err_TFI_T(1:length([5:skp:dim])),[],1)
    line([0 100],[FI_TRUE FI_TRUE],'Color','r')
    line([27.5 27.5],[-300 400],'Color','b')
      load('ck_neurons.mat')
    shadedErrorBar([5:skp:dim]',bs_TFI_N(1:length([5:skp:dim])),err_TFI_N(1:length([5:skp:dim])),{'g','LineWidth',1},1)
       prettyplot
       axis([10 100 -300 300])
       title('choose k that max(F)')
    
    
    
end





end

function [CFI idx FI varCFI] = f_inforad(r_s1,r_s2,ds)
% estimating f'

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

function [CTFI, idx, TFI, varCTFI] = ft_info1(r_s1,r_s2,ds)
%% choose K = T
% output [CorrectTruncated_FisherInfo, idx, naiveTFI, varTFI ]

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

function [CTFI, idx, TFI, varCTFI] = ft_info2(r_s1,r_s2,ds)
%% choose 99% variance
% output [CorrectTruncated_FisherInfo, idx, naiveTFI, varTFI ]

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
var_vec = flipud(diag(EVal));
pct = cumsum(var_vec) / sum(diag(EVal));% Calculate the cumulative percentage of the variances.
 %Find the 99% variance
idx = find(pct >= 0.99, 1);
%figure(1);hold on;plot(E)

%A = V(:,1:ceil(idx-5))'; % low-d subspace
A = V(:,1:idx)'; % low-d subspace

g_prime = (A*f_prime')';
V_Sigma = A*Sigma*A';

% estimating Fischer information
TFI = (g_prime/V_Sigma)*g_prime'; %naive estimate

% bias-correction (derive this!)
K = idx;

CTFI = TFI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*ds^2));
varCTFI =  (2*CTFI^2)/(2*T-K-5) * (1 + 4*(2*T-3)/(T*CTFI*ds^2) + 4*K*(2*T-3)/(T*CTFI*ds^2)^2); % variance of FI bias corrected estimator

%from inv degrees to radians
%CTFI = 1/(((1/CTFI)*(pi/180)^2));
%varCTFI= 1/(((1/varCTFI)*(pi/180)^2));
%TFI = 1/(((1/TFI)*(pi/180)^2));

end


function [CTFI, idx, TFI, varCTFI] = ft_info3(r_s1,r_s2,ds)
%% choose K that max FI


% output [CorrectTruncated_FisherInfo, idx, naiveTFI, varTFI ]

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

for idx = 1:(2*size(r_s1,1))
    
    A = V(:,1:idx)'; % low-d subspace

    g_prime = (A*f_prime')';
    V_Sigma = A*Sigma*A';

    
% estimating Fischer information
TFI(idx) = (g_prime/V_Sigma)*g_prime'; %naive estimate

% bias-correction (derive this!)
K = idx;

CTFI(idx) = TFI(idx)*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*ds^2));
varCTFI(idx) =  (2*CTFI(idx)^2)/(2*T-K-5) * (1 + 4*(2*T-3)/(T*CTFI(idx)*ds^2) + 4*K*(2*T-3)/(T*CTFI(idx)*ds^2)^2); % variance of FI bias corrected estimator

    
    
end

%[r,c]=find(CTFI==median(CTFI));
[minValue,c] = min(abs(CTFI-median(CTFI)));

%line([c c],[min(CTFI_final)-10 max(CTFI_final)+10],'Color','r')

CTFI = CTFI(c);
idx = c;
TFI = TFI(c);
varCTFI = varCTFI(c);
end