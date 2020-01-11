


%% FT-ize data

rawData             = [];
rawData.fsample     = MER_data(1).metaData(1).SampFreq;
rawData.label       = MER_data(1).label;
rawData.trial       = {MER_data(1).data};

len                 = length(rawData.trial{1});
tot_time            = len/rawData.fsample;
time_vec            = 0:1/rawData.fsample:(len-1)/rawData.fsample;
rawData.time        = {time_vec};

%% Parallelize ft_preprocessing
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

%%

time_tot    = floor(preprocLFP.time{1}(end));
fs          = preprocLFP.fsample;
sample_tot  = time_tot*fs;
trl         = [];
trl         = [[[0:fs/2:sample_tot-(fs/2)]+1]' [fs/2:fs/2:sample_tot]'];
trl(:,3)    = [fs.*size(trl,1)/2]';

% find the interesting epochs of data
cfg = [];
cfg.trl = trl;
dataEpochLFP = ft_redefinetrial(cfg, preprocLFP);
dataEpochMER = ft_redefinetrial(cfg, preprocBG);


%% FT Visclean routine 

cfg                 = [];
visClean            = ft_rejectvisual(cfg,dataEpochMER);

%%
cfg              = [];
% cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
% cfg.method       = 'wavelet';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foi          = [400:25:4000];                          
cfg.tapsmofrq    = 0.4 *cfg.foi;
cfg.t_ftimwin    = 5./cfg.foi;   
cfg.pad          = 2;
% cfg.width        = 5; 
% cfg.toi          = 0:0.02:0.5;            
LFP_power          = ft_freqanalysis(cfg, visClean);

%%
cfg                         = [];
cfg.zlim                    = 'maxmin';
cfg.xlim                    = [0 0.1]; 
cfg.ylim                    = [1 400];
cfg.maskstyle               = 'saturation';	
cfg.channel                 = {'Central'};
cfg.shading                 = 'interp';
cfg.masknans                = 'yes';
cfg.shading                 = 'interp';
ft_singleplotER(cfg,LFP_power);

%% Plot data

plot(visClean.time{1},visClean.trial{2}(1,:))

semilogy(LFP_power.freq, LFP_power.powspctrm);

%%
A=fileread('protokoll.txt');
B = strfind(A, 'channel configuration');
what = {A(B(1)+23:B(1)+82)};
what_tok = tokenize(char(what),':*,');
channel_config = {what_tok{2},what_tok{4},what_tok{6},what_tok{8},what_tok{10}}

%%
A=fileread('protokoll.txt');
B = strfind(A, 'New site no');
stim = strfind(A, 'Stimulation');
B = B(find(B<stim(1)));
depth=[];
for loop = 1:length(B)
    char_tok = tokenize(A(B(loop)+28:B(loop)+40), ': ');
    if ~isempty(char_tok{3})
        depthEx(loop) = str2num(char_tok{3});
    else
        depthEx(loop) = -20;
    end
end

depth = [min(depthEx(find(depthEx~=-20))):1:-6, -5:0.5:max(depthEx)];

%%
what_tok = tokenize(char(what),':*,');
channel_config = {what_tok{2},what_tok{4},what_tok{6},what_tok{8},what_tok{10}}

%%

A = fileread('protokoll.txt');

for loop = 1:length(MER_data)
   
    B = strfind(A, {MERData(1).filenames(1)});

end

%%


