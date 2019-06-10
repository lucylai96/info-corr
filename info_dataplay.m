function info_dataplay
% purpose: look at harvey lab data

%deResp: matrix of trials x neurons structure of deconvolved spike sum proxies from calcium traces
%mvSpd: vector of movement speeds (abstract number, larger = faster), one per trial
%visCon: vector of stimulus constrasts, one per trial
%visOri: vector of stimulus drift directions, one per trial
addpath('/Users/lucy/Google Drive/Harvard/Rotations/drugo-lab/')
addpath('/Users/lucy/Google Drive/Harvard/MatTools/')

setPretty

processed =0; %first step (truncate trials)
processed =1; %second step (truncate neurons and trials)
processed =2; %plot stuff
processed =3; %choose dimensionality and then calculate based on both trials and neurons
%processed =4; %just plotting information against trials
if processed ==0
    load('j1_171018.mat')
    
    %% take only data from the highest contrasts for now
    
    highConResp = deResp(visCon==max(visCon),:);
    highConMov = mvSpd(visCon==max(visCon));
    highVisOri = visOri(visCon==max(visCon));
    ori = unique(visOri); %[0, 15, 30, 45, 60, 75, 90];
    
    %% separate out data by oriations
    for i  = 1:length(ori)
        temp = highConResp(highVisOri==ori(i),:); % sort by oriation
        highConRespOri(:,:,i)= trunc(temp, 800); %truncate all conditions to 800 trials, this is the matrix with all high contrast responses separated by oriation
    end
    
    % ~500 neurons and 800 trials
    
    %% Fisher information
    % take two diff oriations 15 deg apart, calculate the FI
    
    for i  = 1:length(ori)-1 %comparisons are # of oriations - 1
        FI(i) = f_info(highConRespOri(:,:,i), highConRespOri(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
        
    end
    
    names = {'0-15^o'; '15-30^o'; '30-45^o'; '45-60^o'; '60-75^o';'75-90^o'};
    figure; plot(FI,'LineWidth',1.5)
    set(gca,'xtick',[1:length(FI)],'xticklabel',names)
    xlabel('\delta\theta')
    ylabel('Fisher Information (rad^{-2})')
    prettyplot
    %f_info(highConResp15, highConResp60, 15, 60)
    
    %[U,S,V] = svd(Sigma);
    
    
    %% truncated data (how many trials can you take out before the F info differs)
    % has to be less trials for both conditions
    
    k = 1;
    for t = 100:50:250
        
        for b  = 1:100 %bootstrap 100 times
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
            
            % bootstrapping to get a reduced matrix of t trials
            for i  = 1:length(ori)
                t_highConRespOri_bsTrials(:,:,i) = datasample(highConRespOri(:,:,i),t,1); % 3rd param: 1=subsampling trials, 2=subsampling neurons
            end
            
            for i  = 1:length(ori)-1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                TFI_trials(b,i) = ft_info(t_highConRespOri_bsTrials(:,:,i), t_highConRespOri_bsTrials(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                
            end
            
        end
        bs_TFI(k,:,1) = mean(TFI_trials,1); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI(k,:,1) = sem(TFI_trials,1);
        
        k = k+1;
        % bs_TFI = %bootstrapped data
    end
    
elseif processed ==1
    load('j1_171018_estFI.mat')
    bs_TFI_NT = zeros(5,6,8);
    err_TFI_NT = zeros(5,6,8);
    k = 1;
    figure;
    rp1=randperm(800);
    rp2=randperm(500);
    for t = 100:100:800 %for these 8 truncated # of trials
        
        for n = 100:100:500 %for these 5 truncated # of neurons
            TFI_NT = zeros(50,6);
            
            for b  = 1:50 %bootstrap 100 times
                
                
                t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
                t_highConRespOri_bsNT = zeros(t,n,length(ori));
                
                
                for i  = 1:length(ori)
                    % bootstrapping to get a reduced matrix of t trials
                    t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);%datasample(highConRespOri(:,:,i),t,1, 'Replace', false); % 3rd param: 1=subsampling trials, 2=subsampling neurons
                    % bootstrapping to get a reduced matrix of n neurons
                    t_highConRespOri_bsNT(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i);
                end
                
                for i  = 1:length(ori)-1 %comparisons are # of oriations - 1
                    %here the columns are the diff oriation difference
                    %conditions and the rows are the bootstraps
                    if t>(n+5)/2 %if you have more trials than neurons
                        TFI_NT(b,i) = f_info(t_highConRespOri_bsNT(:,:,i), t_highConRespOri_bsNT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    else %less trials than neurons
                        TFI_NT(b,i) = ft_info(t_highConRespOri_bsNT(:,:,i), t_highConRespOri_bsNT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    end
                end
                
            end
            %each row here is the condition of # of neurons (should have 5 rows in
            %the end 100:100:500. each "page" is the truncated # of trials (should
            %be 8 in the end for 100:100:800)
            bs_TFI_NT(n*.01,:,k) = mean(TFI_NT,1); %here is taking mean across all conditions (should also do for just one condition at a time)
            err_TFI_NT(n*.01,:,k) = sem(TFI_NT,1);
            
            
        end
        
        colors_p = gradientCol(size(bs_TFI_NT(:,:,k),2),3);
        subplot(2,4,k); hold on;
        set(gca, 'ColorOrder',  colors_p)
        plot(repmat([100:100:500]',1,6),bs_TFI_NT(:,:,k), 'LineWidth',2)
        title(strcat("#trials = ",num2str(t)))
        if k==5
            xlabel('# neurons')
            ylabel('estimated Fisher Information (rad^{-2}')
            legend('0-15^o', '15-30^o', '30-45^o', '45-60^o','60-75^o','75-90^o')
        end
        prettyplot
        
        k = k+1;
        
    end
    
    
    for n = 1:size(bs_TFI_NT,1)
        errorbar(FI,bs_TFI_NT(n,:), err_TFI_NT(n,:),'.','MarkerSize',20,'LineWidth',1,'Color', colors_p(k,:))
        k = k+1;
        % bs_TFI = %bootstrapped data
    end
    
elseif processed==2
    load('j1_171018_estFI.mat')
    
    %% estimated FI from full information vs estimated FI from estimated information
    colors_p = gradientCol(size(bs_TFI_NT,3),2);
    figure; hold on;
    
    % n =5;
    plt = 1;
    for k = 1:size(bs_TFI_NT,3)
        for n = 1:size(bs_TFI_NT,1)
            subplot(1,5,n); hold on;
            errorbar(FI,bs_TFI_NT(n,:,k),err_TFI_NT(n,:,k),'.','MarkerSize',15,'LineWidth',1,'Color', colors_p(k,:))
            prettyplot
            
            % axis square
            axis([40 120 -30 150])
            
            title(strcat("#neurons = ",num2str(n*100)))
            
        end
        
    end
    equalabscissa(1,5)
    
    for n = 1:size(bs_TFI_NT,1)
        subplot(1,5,n); hold on;
        dline
    end
    
    xlabel('Fisher Information from full data (rad^{-2})')
    ylabel('Fisher Information from truncated dataset(rad^{-2})')
    legend('T=100','T=200','T=300','T=400','T=500','T=600','T=700','T=800','unity')
    
    
    %% FI vs trials
    figure;hold on;
    NT = reshape(mean(bs_TFI_NT,2),5,8);
    NT_err = reshape(mean(err_TFI_NT,2),5,8); %sem from bootstrapping for each neuron/trial condition
    %collapse over all 6 degree differences
    for n = 1:size(bs_TFI_NT,1)
        subplot(1,5,n); hold on;
        shadedErrorBar(100:100:800,NT(n,:),NT_err(n,:),[],1)
        title(strcat("#neurons = ",num2str(n*100)))
        plot([0 800],[NT(5,8) NT(5,8)],'r--')
        if n == 1
            xlabel('# trials')
            ylabel('Fisher Information (rad^{-2})')
        end
        prettyplot
    end
    
    equalabscissa(1,5)
    
    %% FI vs neurons
    
    figure;hold on;
    % NT = reshape(mean(bs_TFI_NT,2),5,8);
    % NT_err = reshape(mean(err_TFI_NT,2),5,8);
    %collapse over all 6 degree differences
    for t = 1:size(bs_TFI_NT,3)
        subplot(2,4,t); hold on;
        shadedErrorBar(100:100:500,NT(:,t),NT_err(:,t),[],1)
        title(strcat("#trials = ",num2str(t*100)))
        if t == 5
            xlabel('# neurons')
            ylabel('Fisher Information (rad^{-2})')
        end
        prettyplot
    end
    
    equalabscissa(2,4)
    
    
    %% neurons and trials with error bars on both sides (doesn't look good)
    
    %      NT = reshape(mean(bs_TFI_NT,2),5,8);
    %      NT_err = reshape(mean(err_TFI_NT,2),5,8);
    %      mean(NT_err,1)
    %
    %       colors_p = gradientCol(size(bs_TFI_NT,3),2);
    %     figure; hold on;
    %
    %     % n =5;
    %     plt = 1;
    %     for k = 1:size(bs_TFI_NT,3)
    %         for n = 1:size(bs_TFI_NT,1)
    %             subplot(1,5,n); hold on;
    %             errorbar(mean(FI),NT(n,k),NT_err(n,k),NT_err(n,k),NT_err(n,k),NT_err(n,k),'.','MarkerSize',15,'LineWidth',1,'Color', colors_p(k,:))
    %             prettyplot
    %
    %             % axis square
    %             axis([40 120 -30 150])
    %
    %             title(strcat("#neurons = ",num2str(n*100)))
    %
    %         end
    %
    %     end
    %     equalabscissa(1,5)
    %
    %     for n = 1:size(bs_TFI_NT,1)
    %         subplot(1,5,n); hold on;
    %         dline
    %     end
    %
    %     xlabel('Fisher Information from full data (rad^{-2})')
    %     ylabel('Fisher Information from truncated dataset(rad^{-2})')
    %     legend('T=100','T=200','T=300','T=400','T=500','T=600','T=700','T=800','unity')
    %
    
    
    %% error vs number of trials
    
    figure;hold on;
    
    for n = 1:size(bs_TFI_NT,1)
        subplot(1,5,n); hold on;
        sqerr(n,:) = (NT(n,:)-mean(FI)).^2;
        plot(100:100:800,sqerr(n,:))
        title(strcat("#neurons = ",num2str(n*100)))
        if n == 1
            xlabel('# trials')
            ylabel('Squared error')
        end
        prettyplot
    end
    
    equalabscissa(1,5)
    
    
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
    
elseif processed==3 % choose a dimensionality, and reduce either neurons or trials and look at FI.
    
    load('j1_171018_estFI.mat')
    skp= 60; % step size
    k = 1;
    b = 50; %bootstraps
    sz = length(10:skp:400); % size of the calculated vectors
    bs_TFI_N = zeros(sz,6);
    err_TFI_N = zeros(sz,6);
    bs_naiveTFI_N = zeros(sz,6);
    err_naiveTFI_N = zeros(sz,6);
    var_TFI_N = zeros(sz,6);
    
    bs_TFI_T = zeros(sz,6);
    err_TFI_T = zeros(sz,6);
    bs_naiveTFI_T = zeros(sz,6);
    err_naiveTFI_T = zeros(sz,6);
    var_TFI_T = zeros(sz,6);
    
    TFI_N = zeros(b,6);
    naiveTFI_N = zeros(b,6);
    varTFI_N = zeros(b,6);
    TFI_T = zeros(b,6);
    naiveTFI_T = zeros(b,6);
    varTFI_T = zeros(b,6);
    
    neuron =0;
    for dim = 10:skp:800 %for these dimensionalities
        %% do by cutting off NEURONS
        if neuron == 1
        t = 500; %use all the trials
        n = dim; %only use "dim" # of neurons
        for b  = 1:50 %bootstrap 100 times
            rp1=randperm(800);
            rp2=randperm(500);
            
            t_highConRespOri_bsTrials = zeros(t,size(highConRespOri,2),length(ori));
            t_highConRespOri_bsN = zeros(t,n,length(ori));
            for i  = 1:length(ori)
                % bootstrapping to get a reduced matrix of t trials
                t_highConRespOri_bsTrials(:,:,i) = highConRespOri(rp1(1:t),:,i);
                % bootstrapping to get a reduced matrix of n neurons
                t_highConRespOri_bsN(:,:,i)= t_highConRespOri_bsTrials(:,rp2(1:n),i);
                
            end
            
            for i  = 1:length(ori)-1 %comparisons are # of orientations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                % if  t>(n+2)/2 %if you have more trials than neurons
                [TFI_N(b,i) idx(b,i) naiveTFI_N(b,i) varTFI_N(b,i)] = f_info(t_highConRespOri_bsN(:,:,i), t_highConRespOri_bsN(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                % end
            end
        end
        %each row is a dimension (# neurons)
        bs_TFI_N(k,:) = mean(TFI_N,1); % here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_N(k,:) = std(TFI_N,[],1);
        var_TFI_N(k,:) = mean(varTFI_N,1); %any variance on this reflects the variance of the 6 orientation difference conditions and bootstrapping
        
        idx_TFI_N(k,:) = mean(idx,1);
        
        bs_naiveTFI_N(k,:) = mean(naiveTFI_N,1);
        err_naiveTFI_N(k,:) = std(naiveTFI_N,[],1);
        
        end
        %figure;shadedErrorBar([10:10:500]',mean(bs_TFI_N,2),std(bs_TFI_N,[],2),[],1)
        
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:50 %bootstrap 100 times
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
            
            for i  = 1:length(ori)-1 %comparisons are # of oriations - 1
                %here the columns are the diff oriation difference
                %conditions and the rows are the bootstraps
                n = 500;
               % if t>(n+5)/2 && t>=n %if you have more trials than neurons
                if t>(n+5)/2 %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    %elseif t>(n+5)/2 && t<n
                    %elseif t<n %removing neurons all the way down
                    %    n = t;
                    
                    %     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    
                elseif t<(n+5)/2 %less trials than neurons
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    
                end
            end
        end
        %each row is a dimension (# neurons)
        bs_TFI_T(k,:) = mean(TFI_T,1); %here is taking mean across all conditions (should also do for just one condition at a time)
        err_TFI_T(k,:) = std(TFI_T,[],1);
        var_TFI_T(k,:) = mean(varTFI_T,1); %any variance on this reflects the variance of the 6 orientation difference conditions and bootstrapping
        
        idx_TFI_T(k,:) = mean(idx,1);
        
        bs_naiveTFI_T(k,:) = mean(naiveTFI_T,1);
        err_naiveTFI_T(k,:) = std(naiveTFI_T,[],1);
        
        k = k+1;
    end% dimensions
    %save('j1_new.mat')
   
    
    % FI found with neurons vs FI found with trials
    figure;plot(mean(bs_TFI_N,2),mean(bs_TFI_T,2),'o')
    errorbar(mean(bs_TFI_N,2),mean(bs_TFI_T,2),mean(err_TFI_T,2),mean(err_TFI_T,2),mean(err_TFI_N,2),mean(err_TFI_N,2),'o')
    axis equal; dline;
    xlabel('I_F (reduced dim using neurons');ylabel('I_F reduced dim using trials')
    prettyplot
    
    % FI vs neurons and vs trials
    figure;
    subplot 121
    shadedErrorBar([10:skp:dim]',mean(bs_TFI_N,2),std(bs_TFI_N,[],2),[],1)
    xlabel('# neurons');ylabel('I_F (rad^{-2})')
    title('# trials = 500')
    prettyplot
    subplot 122;
    shadedErrorBar([10:skp:dim]',mean(bs_TFI_T,2),std(bs_TFI_T,[],2),[],1)
    hold on;line([250 250],[-40 120],'Color','r')
    xlabel('# trials');
    title('# neurons = dim')
    prettyplot
    equalabscissa(1,2)
    suptitle('reducing dimensionality by reducing neurons vs. reducing trials')
    
    
    % above plot overlayed one on top of another
    % plot each condition separately
    names = {'0-15^o'; '15-30^o'; '30-45^o'; '45-60^o'; '60-75^o';'75-90^o'};
    figure;hold on;
    for i = 1:6
        subplot (2,3,i)
        shadedErrorBar([10:skp:dim]',bs_TFI_N(:,i),err_TFI_N(:,i),[],1)
        
        hold on;
        shadedErrorBar([10:skp:dim]',bs_TFI_T(:,i),err_TFI_T(:,i),{'ro-','markerfacecolor',[1 0 0]},1)
        hold on;line([250 250],[-40 120],'Color','r')
        title(names{i})
        prettyplot
        
    end
    subplot 234;hold on;
    xlabel('rank');ylabel('I_F (rad^{-2})')
    equalabscissa(2,3)
    
    suptitle('reducing dimensionality by reducing neurons vs. reducing trials')
    prettyplot
    
    % necess dimensionality
    figure;
    subplot 121
    plot([10:skp:dim]',mean(idx_TFI_N,2))
    xlabel('# neurons');ylabel('necessary dimensionality (\lambda_k > 0)')
    title('# trials = 500'); axis equal; axis([0 dim 0 dim]);dline;
    prettyplot
    subplot 122;
    plot([10:skp:dim]',mean(idx_TFI_T,2)); axis equal; axis([0 dim 0 dim]);dline;
    hold on;line([250 250],[0 dim],'Color','r')
    xlabel('# trials');
    title('# neurons = 500')
    prettyplot; equalabscissa(1,2)
elseif processed ==4
    % load('j1_new.mat')
    
    k = 1;
    for dim = 10:skp:800 %for these dimensionalities
        %% here only use 1 orientation difference
        %% do by cutting off TRIALS
        t = dim;  %has to be t<(n+2)/2, new dim is K = 2T-2, only use "dim" # of trials
        % n = 400;  % keep all neurons
        for b  = 1:50 %bootstrap 100 times
            n = 500; % keep all neurons
            rp1=randperm(800);
            rp2=randperm(500);
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
                n = 500;
                if t>(n+5)/2 && t>=n %if you have more trials than neurons
                    
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                elseif t<n
                    n = t;
                    [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                    
                end
                %                 elseif t>(n+5)/2 && t<n
                %                     n = t;
                %
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = f_info(t_highConRespOri_bsT(:,1:n,i), t_highConRespOri_bsT(:,1:n,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
                %
                %                 elseif t<(n+5)/2 %less trials than neurons
                %                     [TFI_T(b,i) idx(b,i) naiveTFI_T(b,i) varTFI_T(b,i)] = ft_info(t_highConRespOri_bsT(:,:,i), t_highConRespOri_bsT(:,:,i+1),ori(i),ori(i+1));%this is the bias-corrected FI in units of rad^-2
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
end



%% comparison with naiveFI
xplot = [10    70   130   190   250   310   370   430  490  550   610   670   730   790];
% use all neurons, plot of FI vs trials
figure; hold on;
shadedErrorBar(xplot,bs_TFI_T(2:end,1),err_TFI_T(2:end,1),{'b','markerfacecolor',[0 0 1]},1)
shadedErrorBar(xplot,mean(bs_naiveTFI_T,2),mean(err_naiveTFI_T,2),{'g','markerfacecolor',[0 1 0]},1)
plot([0 800],[mean(bs_TFI_T(end,:),2) mean(bs_TFI_T(end,:),2)],'k--')
xlabel('# trials')
ylabel('Fisher Information (rad^{-2})')
prettyplot

figure; hold on;
%variance doesn't match?
shadedErrorBar(xplot,mean(bs_TFI_T,2),mean(err_TFI_T,2)*sqrt(6),[],1) %std across bootstrapped samples
% shadedErrorBar(xplot,mean(bs_TFI_T,2),std(bs_TFI_T,[],2),[],1) %std across the 6 different orientation differences
shadedErrorBar(xplot,mean(bs_TFI_T,2),sqrt(mean(var_TFI_T,2)),{'r','markerfacecolor',[1 0 0]},1)

xlabel('# trials')
ylabel('Fisher Information (rad^{-2})')

prettyplot

figure;
shadedErrorBar([10:skp:dim]',mean(bs_naiveTFI_T,2),std(bs_naiveTFI_T,[],2),[],1)
hold on;shadedErrorBar([10:skp:dim]',mean(bs_TFI_T,2),std(bs_TFI_T,[],2),{'g','markerfacecolor',[0 1 0]},1)

end



%
% function [CFI idx] = f_info(r_s1,r_s2,s1,s2)
% % estimating f'
% f_prime = (mean(r_s2)-mean(r_s1))./(s2-s1);
%
% % estimating sigma
% Sigma = 0.5*(cov(r_s1) + cov(r_s2));
%
% [E] = sort(eig(Sigma),'descend');
% idx = max(find(E<0.001,1),0);
% if isempty(idx)
%     idx=0;
% end
%
% % estimating Fischer information
% FI = (f_prime/Sigma)*f_prime'; %f_prime*inv(Sigma)*f_prime';
%
% % bias-correction
% T = size(r_s1,1);
% N = size(r_s1,2);
% CFI = FI*((2*T-N-3)/(2*T-2)) - ((2*N)/(T*(s2-s1)^2));
% CFI = 1/((((1/CFI)*(pi/180)^2)));
% end
%
% function [CTFI, idx] = ft_info(r_s1,r_s2,s1,s2)
% %estimating from a truncated matrix
% T = size(r_s1,1);
% % estimating f'
% f_prime = (mean(r_s2)-mean(r_s1))./(s2-s1);
%
% % estimating sigma
% Sigma = 0.5*(cov(r_s1) + cov(r_s2));
%
% % singular decomposition
% [E] = sort(eig(Sigma),'descend');
% idx = find(E<0.001,1);
% %idx = (size(r_s1,1) + size(r_s2,1))-5; %cut off eig at K=2T-5
% %idx = T;
% %[U,S,V] = svd(Sigma);
% [EVect, EVal] = eig(Sigma);
% V = flip(EVect,2);
% %var_vec = flipud(diag(EVal));
% %pct = cumsum(var_vec) / sum(diag(EVal));% Calculate the cumulative percentage of the variances.
% % Find the 99% variance
% %idx = find(pct >= 0.99, 1);
% %figure(1);hold on;plot(E)
%
% %A = V(:,1:ceil(idx-5))'; % low-d subspace
% A = V(:,1:ceil(T))'; % low-d subspace
%
% g_prime = (A*f_prime')';
% V_Sigma = A*Sigma*A';
%
% % estimating Fischer information
% TFI = (g_prime/V_Sigma)*g_prime';
%
% % bias-correction (derive this!)
% K = size(g_prime,2);
%
% CTFI = TFI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(s2-s1)^2));
% CTFI = 1/(((1/CTFI)*(pi/180)^2)); %from inv degrees to radians
%
%
% end

function [CTFI,c] = ft_info_test(r_s1,r_s2,s1,s2)
%estimating from a truncated matrix
T = size(r_s1,1);
% estimating f'
f_prime = (mean(r_s2)-mean(r_s1))./(s2-s1);

% estimating sigma
Sigma = 0.5*(cov(r_s1) + cov(r_s2));

% singular decomposition
[E] = sort(eig(Sigma),'descend');
%idx = find(E<0.001,1);
idx = (size(r_s1,1) + size(r_s2,1))-5; %cut off eig at K=2T-5

%[U,S,V] = svd(Sigma);
[EVect, EVal] = eig(Sigma);
V = flip(EVect,2);
%var_vec = flipud(diag(EVal));
%pct = cumsum(var_vec) / sum(diag(EVal));% Calculate the cumulative percentage of the variances.
% Find the 99% variance
%idx = find(pct >= 0.99, 1);

% figure
for idx = 1:(2*size(r_s1,1))-2
    
    A = V(:,1:idx)'; % low-d subspace
    %
    %
    g_prime = (A*f_prime')';
    V_Sigma = A*Sigma*A';
    %
    % estimating Fischer information
    TFI(idx) = (g_prime/V_Sigma)*g_prime';
    %   TFI_rad(idx) = 1/((((1/TFI(idx))*(pi/180)^2)));
    
    % bias-correction (derive this!)
    K = size(g_prime,2);
    
    
    CTFI(idx) = TFI(idx)*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(s2-s1)^2));
    CTFI_final(idx) = 1/((((1/CTFI(idx))*(pi/180)^2)));
    
    %CTFI = TFI*((2*T-K-3)/(2*T-2)) - ((2*K)/(T*(s2-s1)^2));
    %CTFI = 1/(((1/CTFI)*(pi/180)^2)); %from inv degrees to radians
    
    
end
%plot(CTFI_final)

%xlabel('k')
%ylabel('Fisher Information')
%prettyplot
%
[r,c]=find(CTFI_final==max(CTFI_final));
%line([c c],[min(CTFI_final)-10 max(CTFI_final)+10],'Color','r')

CTFI = CTFI_final(c);
%CTFI_final(2*T-5);

end
