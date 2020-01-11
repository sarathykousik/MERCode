

for loop = 1: length(MER_data)
    % FT-ize data
    disp(['##########  Depth: ', num2str(loop), '/',num2str(length(MER_data))])
    rawData             = [];
    rawData.fsample     = MER_data(loop).metaData(1).SampFreq;
    rawData.label       = MER_data(loop).label;
    rawData.trial       = {MER_data(loop).data};

    len                 = length(rawData.trial{1});
    tot_time            = len/rawData.fsample;
    time_vec            = 0:1/rawData.fsample:(len-1)/rawData.fsample;
    rawData.time        = {time_vec};

    % Parallelize ft_preprocessing
    cfg                         = [];
    cfg.lpfilter                = 'yes';
    cfg.lpfreq                  =  4000;
    cfg.hpfilter                = 'yes';
    cfg.hpfreq                  =  400;
    preprocBG                   = ft_preprocessing(cfg, rawData);

    %
    cfg                         = [];
    cfg.bpfilter                = 'yes';
    cfg.bpfreq                  =  [2 400];
    cfg.bpinstabilityfix        = 'reduce';
    cfg.dftfilter               = 'yes';
    cfg.dftfreq                 = [50 100 150 200 250 300 350];
    cfg.detrend                 = 'yes';
    preprocLFP                  = ft_preprocessing(cfg, rawData);

    %

    time_tot    = floor(preprocLFP.time{1}(end));
    fs          = preprocLFP.fsample;
    sample_tot  = time_tot*fs;
    trl         = [];
    trl         = [[[0:fs/2:sample_tot-(fs/2)]+1]' [fs/2:fs/2:sample_tot]'];
    trl(:,3)    = [fs.*size(trl,1)/2]';

    % find the interesting epochs of data
    cfg                 = [];
    cfg.trl             = trl;
    dataEpochLFP        = ft_redefinetrial(cfg, preprocLFP);
    dataEpochMER        = ft_redefinetrial(cfg, preprocBG);


    % FT Visclean routine 

    cfg                    = [];
    visCleanLFP{loop}      = ft_rejectvisual(cfg,dataEpochLFP);
    visCleanMER{loop}      = ft_rejectvisual(cfg,dataEpochMER);

end

%%

cfg              = [];
% cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
% cfg.method       = 'wavelet';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';

% cfg.tapsmofrq    = 0.4 *cfg.foi;
% cfg.t_ftimwin    = 5./cfg.foi;   
cfg.pad          = 2;
% cfg.width        = 5; 
% cfg.toi          = 0:0.02:0.5;            
for loop = 1:length(visCleanMER)
    
    disp(['##########  Depth: ', num2str(loop), '/',num2str(length(MER_data))])
    cfg.foi                 = [400:25:4000];                          
    MER_power{loop}         = ft_freqanalysis(cfg, visCleanMER{loop});
    
    MER_SpcMed(loop, :)     = median(...
                                MER_power{loop}.powspctrm(:,find(MER_power{loop}.freq>500 ...
                                    & MER_power{loop}.freq<2000)),2);
    
    for trlLoop = 1:length(visCleanMER{loop}.trial)
        medVal(trlLoop,:) = median(visCleanMER{loop}.trial{trlLoop},2);
        RMSVal(trlLoop,:) = sqrt(mean(visCleanMER{loop}.trial{trlLoop}.^2, 2));
    end
                                
    MER_MedVal(loop, :)     = median(medVal,1);
    MER_RMSVal(loop, :)     = sqrt(mean(RMSVal.^2,1));

        
    cfg.foi          = [1:2:400];                          
    LFP_power{loop}  = ft_freqanalysis(cfg, visCleanLFP{loop});    
end

%%
for dLoop = 1:length(MER_power)
    MERim(dLoop,:,:) = MER_power{dLoop}.powspctrm;
    LFPim(dLoop,:,:) = LFP_power{dLoop}.powspctrm;
end

figure(1)
for chan = 1:5
    subplot(1,5,chan)
    imagesc(MER_power{1}.freq,depths,squeeze(MERim(:,chan,:)))
    title(MER_power{1}.label{chan})
    xlabel('Frequency in Hz')
    ylabel('Depth in mm')
end
suptitle('MER')
%%
foi = 60:200;
figure(2)
for chan = 1:5
    subplot(1,5,chan)
    imagesc(LFP_power{1}.freq(foi),depths,squeeze(LFPim(:,chan,foi)))
    title(LFP_power{1}.label{chan})
    xlabel('Frequency in Hz')
    ylabel('Depth in mm')
end
suptitle('LFP')

%% Norm channels
valMat = MER_SpcMed;
normVal.SpcMed = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);
valMat = MER_RMSVal;
normVal.RMSVal = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);
valMat = MER_MedVal;
normVal.MedVal = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);


% Plot MER through depths
disp([' ' ])
disp('###  Plots   ###')
disp([' ' ])
for  loop = 1:length(MER_data)
   
    depths(loop) = MER_data(loop).depth;
    
end

h1 = figure(1);
subplot 321
plot(depths,MER_SpcMed, '-o', 'LineWidth',2,...
    'MarkerSize',6)
title('Median of spectral power')
legend(visCleanMER{1}.label, 'Location', 'northwest')

subplot 323
plot(depths,MER_MedVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('Median of MER')

subplot 325
plot(depths,MER_RMSVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('RMS of MER')

% h1 = figure(2);
subplot 322
plot(depths,normVal.SpcMed, '-o', 'LineWidth',2,...
    'MarkerSize',6)
title('Median of spectral power')
legend(visCleanMER{1}.label, 'Location', 'northwest')

subplot 324
plot(depths,normVal.MedVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('Median of MER')

subplot 326
plot(depths,normVal.RMSVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('RMS of MER')

suptitle('No-norm   < - - >    Normalized plots')





