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

highConResp15 = highConResp(highVisOri==15,:);

highConResp30 = highConResp(highVisOri==30,:); 

highConResp45 = highConResp(highVisOri==45,:); 

highConResp60 = highConResp(highVisOri==60,:); 

highConResp75 = highConResp(highVisOri==75,:); 

highConResp90 = highConResp(highVisOri==90,:);


end
