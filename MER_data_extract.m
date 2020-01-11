function [MERData, MERlength] = MER_data_extract(directory)
% ************************************************************************* 
% Preprocessing of MER data for analysis
% Dependacy: Uses 'MER_extract' - modified from function provided by INOMED
%                                 MER_bigFile   
% ************************************************************************* 
% [MER_data] = MER_data_extract(directory)
% Ouput per hemisphere
% Import - data, metaData, details of each channel, dataset, site (depth)
% Arrangement in struct: 
%   struct(site)   
%   struct(site).data                     = [chan x samples]      
% 
%   struct(site).metaData(chan).SiteNr    =  Site        (extracted from dat file)   
%   struct(site).metaData(chan).KanalNr)  =  Channel no. (extracted from dat file)   
%   struct(site).metaData(chan).op_id)    =  Unique per hemisphere(extracted from dat file)    
%   struct(site).metaData(chan).SampFreq) =  Sampling frequency(extracted from dat file)        
% 
%   struct(site).filenames                =  .dat file names
%   struct(site).directory                =  Parent directory 
% Written ~kss~ on 11-04-2014               

%%
cd(directory)
MER_filenames               = dir('*MER*');
noMERFiles                  = length(dir('*MER*'));

noChan                      = length(dir([MER_filenames(1).name(1:2),'*MER*']));
noSites                     = noMERFiles/noChan;

MERData = [];
MERData(noSites,1).data   = [];

EMGData = [];
EMGData(noSites,1).data   = [];
%
ft_progress('init', 'gui',     'Importing MER data...');
for loop = 1:noSites
 
    ft_progress(loop/noSites, 'Processing depths %d from %d', loop, noSites);
    
    MERchanFiles            = dir(['*',num2str(loop,'%02d'),'*MER*']);
    EMGchanFiles            = dir(['*',num2str(loop,'%02d'),'*EMG*']);
    
    for chanLoop = 1:length(MERchanFiles)
        
        [MERData(loop).data(chanLoop,:), SampFreq, MERData(loop).metaData(chanLoop)] = ...
                    ExtractMERdat(MERchanFiles(chanLoop).name,directory);
                
        MERData(loop).filenames                = {MERchanFiles(:).name};
        MERData(loop).directory                = directory;
        
        MERlength(loop)                              = length(MERData(loop).data)/SampFreq;

    end
    
end
     ft_progress('close')
%
return