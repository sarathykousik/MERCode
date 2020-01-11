clear all;   clc;    close all;
ft_defaults
directory = 'J:\Surgery\Busk Gunner\391946732';
[MER_data] = MER_data_extract(directory);

%% Filter

noSites = length(MER_data);
MER_filt = MER_data;

for loop = 1:noSites
   disp(['********** Depth: ', num2str(loop)]);
   Fs = MER_data(loop).metaData(1).SampFreq;
   Fbp  = [300, 3000];    N = 2;      type = 'but';    dir = 'twopass';
   MER_filt(loop).data = ft_preproc_bandpassfilter(MER_data(loop).data, Fs, Fbp, N, type, dir);
   

   
   Fbp  = [1, 300]; 
   LFP_filt(loop).data = ft_preproc_bandpassfilter(MER_data(loop).data, Fs, Fbp, N, type, dir);
   
end

%% 
for loop = 1:length(MER_filt)
    t=0:(1/Fs):((length(MER_filt(loop).data )-1)/Fs);
    subplot(5,5,loop)
    plot(t,MER_filt(loop).data')
    title(num2str(loop))
    axis([2 50 -0.2 0.2])
end
suptitle('MER data site-wise')

%% Choose the middle 10 secs
parfor loop  = 1:length(LFP_filt)
    
   half_t = length(LFP_filt(loop).data)/20000*.5 ;
%    clean_t = [half_t-5 half_t+5]; 
   clean_samples = ceil(20000*[half_t-5 half_t+5]);
   
   
   
end



%% Artefact rejection



%% Calculate Beta power



%% Calculate BG activity




















