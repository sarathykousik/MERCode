function [MER_data] = MER_preproc_v1(directory)
%% Preprocessing of MER data for analysis - to be converted into a function
% Dependacy: Uses 'MER_extract' - modified from function provided by INOMED
%                                 MER_bigFile   
% ************************************************************************* 
% All data ouput per hemisphere
% Import - data, metaData, details of each channel, dataset, site (depth)
% Arrange data in a struct array: struct(site)   
%                                 struct(site).trajPos  = [C L A M P]
%                                 struct(site).MERdata  = [channels x samples]      
%                                 struct(site).Fs       = sampling_freq, in Hz
%                                 struct(site).filename = filename      
%                                 struct(site).folder   = parent_folder

%%
cd(directory)

noMERFiles                  = length(dir('*MER*'));
noChan                      = length(dir(['01*MER*']));
noSites                     = noMERFiles/noChan;

MER_data = [];
MER_data(noSites,noChan).data   = [];


for loop = 1:noSites
    disp(['#################### Depth: ', num2str(loop)])
    MERchanFiles            = dir(['*',num2str(loop,'%02d'),'*MER*']);
    
    for chanLoop = 1:length(MERchanFiles)
        disp(['*** MER channel: ', num2str(chanLoop),' ', MERchanFiles(chanLoop).name]);
        
        [MER_data(loop, chanLoop).data, SampFreq, MER_data(loop, chanLoop).metaData] = ...
            MER_dat(MERchanFiles(chanLoop).name,directory);
        MER_data(loop, chanLoop).filename = MERchanFiles(chanLoop).name;
        MER_data(loop, chanLoop).directory = directory;
    end
    
end

return
%%












