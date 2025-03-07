function info_datacomp
clear all
comp = 0;
%% This code loads estimates Fisher information using Harvey Lab dataset with
%% the direct bias corrected estimator, CV estimator with stopping, then makes the main figures of the paper
addpath('/Users/lucy/Google Drive/Harvard/Rotations/drugo-lab/')
addpath('/Users/lucy/Google Drive/Harvard/MatTools/')
setPretty
processed =1;
if processed ==0
    load('j1_171018_estFI_matchRank.mat')
    
    %%% DECODER WITH EARLY STOPPING
    
    skp= 60; % step size
    k = 1;
    ds = 15*(pi/180) ;
    boot = 10; %bootstraps
    sz = length(10:skp:400); % size of the calculated vectors
    
    NR=20; % how many cross-validation splits for early stopping
    
    fracTR=1/3; %fraction training
    fracTE=1/3; %fraction test
    
    FIVAL_ES_N = zeros(sz,6);
    FIVALerr_ES_N = zeros(sz,6);
    FITR_ES_N = zeros(sz,6);
    FITRerr_ES_N = zeros(sz,6);
    
    FIVAL_ES_T = zeros(sz,6);
    FIVALerr_ES_T = zeros(sz,6);
    FITR_ES_T = zeros(sz,6);
    FITRerr_ES_T = zeros(sz,6);
    
    FIVAL_T = zeros(boot,6);
    FIVAL_N = zeros(boot,6);
    FITR_T = zeros(boot,6);
    FITR_N = zeros(boot,6);
    
    
    for dim = 10:skp:800 %for these dimensionalities
        %% do by cutting off NEURONS
        %     t = 500; %use all the trials
        %     n = dim; %only use "dim" # of neurons
        %     for b  = 1:50 %bootstrap 100 times
        %         rp1=randperm(800);
        %         rp2=randperm(500);
        %
        %         t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
        %         t_highConRespOri_bsN = zeros(t,n,length(ori));
        %         for i  = 1:length(ori)
        %             % bootstrapping to get a reduced matrix of t trials
        %             t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);
        %             % bootstrapping to get a reduced matrix of n neurons
        %             t_highConRespOri_bsN(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i);
        %         end
        %
        %         for i  = 1:length(ori)-1 %comparisons are # of orientations - 1
        %             %here the columns are the diff oriation difference
        %             %conditions and the rows are the bootstraps
        %             % if  t>(n+2)/2 %if you have more trials than neurons
        %             [FIVAL_N(b,i), FITR_N(b,i)] = EarlyStopping(t_highConRespOri_bsN(:,:,i), t_highConRespOri_bsN(:,:,i+1),ds,fracTR,fracTE,NR,1);
        %             % end
        %         end
        %     end
        %     %each row is a dimension (# neurons)
        %     FIVAL_ES_N(k,:) = mean(FIVAL_N,1); % here is taking mean across all conditions (should also do for just one condition at a time)
        %     FIVALerr_ES_N(k,:) = sem(FIVAL_N,1);
        %     FITR_ES_N(k,:) = mean(FITR_N,1); %training
        %     FITRerr_ES_N(k,:) = sem(FITR_N,1);
        %
        
        %figure;shadedErrorBar([10:10:500]',mean(bs_ES_N,2),std(bs_ES_N,[],2),[],1)
        
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:boot %bootstrap 100 times
            n = 500; % keep all neurons
            rp1=randperm(800);
            rp2=randperm(500);
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
            t_highConRespOri_bsT = zeros(t,n,length(ori));
            for i  = 1:length(ori)
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i); %t x n
            end
            
            for i  = 1:length(ori)-1 %comparisons are # of orientations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 500;
                [FIVAL_T(b,i), FITR_T(b,i)] = EarlyStopping(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ds,fracTR,fracTE,NR,1);
                
            end
        end
        %each row is a dimension (# trials)
        FIVAL_ES_T(k,:) = mean(FIVAL_T,1); % here is taking mean across all conditions (should also do for just one condition at a time)
        FIVALerr_ES_T(k,:) = sem(FIVAL_T,1);
        FITR_ES_T(k,:) = mean(FITR_T,1); %trainingg
        FITRerr_ES_T(k,:) = sem(FITR_T,1);
        k = k+1;
    end %dimensions # trials
    
    
elseif processed ==1
    load('CV_early_stop.mat')
    load('j1_new.mat')
    load('j1_dimreduce.mat')
    
    
    
    names = {'0-15^o'; '15-30^o'; '30-45^o'; '45-60^o'; '60-75^o';'75-90^o'};figure; hold on;
    figure; hold on;
    for i = 1:6
        subplot (2,3,i); hold on;
        plot([0 800],[mean(bs_TFI_T(end,i)) mean(bs_TFI_T(end,i))],'k--')
        line([252.5 252.5],[-50 300],'Color','k')
        axis([ 0 800 -50 300])
        title(names{i})
        prettyplot
    end
    equalabscissa(2,3)
    subplot 231
    xlabel('# trials')
    ylabel('Fisher Information (rad^{-2})')
    
    for i = 1:6
        subplot (2,3,i); hold on;
        shadedErrorBar([10:skp:dim]',bs_TFI_T(:,i),sqrt(var_TFI_T(:,i)),{'b','markerfacecolor',[0 0 1]},1)
        shadedErrorBar([310:skp:dim]',bs_naiveTFI_T(6:end,i),err_naiveTFI_T(6:end,i),{'g','markerfacecolor',[0 0 1]},1)
        prettyplot
    end
    
    for i = 1:6
        subplot (2,3,i);shadedErrorBar(10:skp:dim,FIVAL_ES_T(:,i),FIVALerr_ES_T(:,i),{'r','markerfacecolor',[1 0 0]},1)
        hold on; shadedErrorBar(10:skp:dim,FITR_ES_T(:,i),FITRerr_ES_T(:,i),{'r','markerfacecolor',[1 0 0]},1)
        % plot([0 800],[80 80],'r--');
        
    end
    
    %
    %     % use all neurons, plot of FI vs trials
    %     figure;shadedErrorBar(10:skp:dim,mean(FIVAL_ES_T,2),mean(FIVALerr_ES_T,2),{'r','markerfacecolor',[1 0 0]},1)
    %     hold on; shadedErrorBar(10:skp:dim,mean(FITR_ES_T,2),mean(FITRerr_ES_T,2),{'r','markerfacecolor',[1 0 0]},1)
    %     plot([0 800],[80 80],'r--');axis([0 800 -50 250])
    %     xlabel('# trials')
    %     ylabel('Fisher Information (rad^{-2})')
    %     prettyplot
    %
    %     xplot = [10    70   130   190   250   310   370   430   490 490  550   610   670   730   790];
    %     shadedErrorBar(xplot,mean(bs_TFI_T,2),mean(err_TFI_T,2)*sqrt(6),{'b','markerfacecolor',[0 0 1]},1)
    %     shadedErrorBar(xplot,mean(bs_naiveTFI_T,2),mean(err_naiveTFI_T,2),{'g','markerfacecolor',[0 1 0]},1)
    %     plot([0 800],[mean(bs_TFI_T(end,:),2) mean(bs_TFI_T(end,:),2)],'k--')
    %     xlabel('# trials')
    %     ylabel('Fisher Information (rad^{-2})')
    %     prettyplot
    %
    
    %% error (MSE)
    
    figure;hold on;
    for i = 1:6
        subplot (2,3,i); hold on;
        bias1 = (bs_TFI_T(end,i)-bs_TFI_T(:,i)).^2;
        bias2 = (bs_TFI_T(end,i)-bs_naiveTFI_T(:,i)).^2;
        bias3 = (bs_TFI_T(end,i)-FIVAL_ES_T(:,i)).^2;
        plot(10:skp:dim,sqrt(bias1+var_TFI_T),'b')
        plot(310:skp:dim,sqrt(bias2(6:end)+(err_naiveTFI_T(6:end,i)).^2),'g')
        plot(10:skp:dim,sqrt(bias3+ (FIVALerr_ES_T(:,i)).^2) ,'r')
        title(names{i})
        line([252.5 252.5],[0 800],'Color','k')
        prettyplot
    end
    subplot(2,3,4)
    
    xlabel('# trials')
    ylabel('RMSE')
    equalabscissa(2,3)
    
    
    %% error vs ratio of trials to neurons
    
    figure;hold on;
    colors_p = gradientCol(size(bs_TFI_NT,1),1);
    
    for n = 1:size(bs_TFI_NT,1)
        sqerr(n,:) = (NT(n,:)-mean(FI)).^2;
        relerr(n,:) = sqrt(sqerr(n,:))./mean(FI);
        tn_ratio(n,:) = [100:100:800]./(n*100);
        error_mat(:,n) =  relerr(n,:);
        plot([100:100:800]./(n*100),relerr(n,:),'o-','LineWidth',1,'Color', colors_p(n,:))
        %title(strcat("#neurons = ",num2str(n*100)))
        if n == 1
            xlabel('trials-to-neurons ratio T/N')
            ylabel('Relative error')
        end
        prettyplot
    end
    legend('N=100','N=200','N=300','N=400','N=500')
    figure;
    imagesc(error_mat)
    xlabel('# neurons (x100)')
    ylabel('trials-to-neurons ratio T/N')
    set(gca,'YDir','normal');
    
    h = colorbar;
    ylabel(h, 'Relative error')
    
end
%% parameters for information estimation

% %%% Variational Bayes Logistic Regression
% FIVBTR = NaN(NTT,NP);
% FIVBVAL = NaN(NTT,NP);
% fracTRVB = 2/3;
% for t = 1:NTT
%     NVAL = NT(t);
%     for p=1:NP
%         rng(p)
%         indtmp = randperm(N); %subsample data to desired size
%         D1 = squeeze(resp(:,or_corr(1),indtmp(1:NVAL)))';
%         D2 = squeeze(resp(:,or_corr(2),indtmp(1:NVAL)))';
%         [FIVBVAL(t,p), FIVBTR(t,p)] = VBLogReg(D1,D2,ds,fracTRVB,NR,1);
%     end
% end

end