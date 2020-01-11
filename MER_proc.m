function [MER_summary] = MER_proc(MER_data, toiChoice)
% ************************************************************************* 
% Preprocessing of MER data processing for extracting summary values
% ************************************************************************* 
% Input:
% MER_data  - structure from MER_data_extract function
% Output:
% Computes three summary values
% 1. RMS of filtered (300 - 3000 Hz) MER data
% 2. Median of filtered (300 - 3000 Hz) MER data
% 3. Median of pwelch estimates (500 - 2000 Hz) MER data

disp([' ' ])
disp(['Calculating power spectra of MER...  ' ])
disp([' ' ])
startTime = cputime;
Fs        = MER_data(1).metaData(1).SampFreq;
fmin_MER  = 300;    
fmax_MER  = 2500;
fmin_LFP  = 4;    
fmax_LFP  = 400;
% options = optimset('Display','off', 'MaxIter', 300 , 'MaxFunEvals', 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 , 'TolFun', 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001, 'TolX', 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
% xrange  = [10:1000 2000:3000];

[fb_MER,fa_MER]   =  ellip(4,0.1,40,[fmin_MER fmax_MER]*2/Fs);
% [fb_LFP,fa_LFP]   =  ellip(2,0.1,40,[fmin_LFP fmax_LFP]*2/Fs);
toi = ceil(toiChoice*Fs);
ft_progress('init', 'gui',     'Processing MER data...');
for dLoop = 1:length(MER_data)
    ft_progress(dLoop/length(MER_data), ...
        'Processing depths %d from %d', dLoop, length(MER_data));

    for chan = 1:size(MER_data(1).data,1)
%        [dLoop chan]
        % MER
        MER_filt =  filtfilt(fb_MER,fa_MER,MER_data(dLoop).data(chan,toi(1):toi(2)));
        
        [Pxx,freq] = pwelch(MER_filt,...
                        hanning(ceil(Fs)),ceil(Fs*0.2),ceil(Fs), Fs);
        MER_SpcMed(dLoop, chan)          = median(Pxx(500:2000));
        MER_MedVal(dLoop, chan)          = median(MER_filt);
        MER_RMSVal(dLoop, chan)          = sqrt(mean(MER_filt.^2));
        MER_Spc(dLoop,chan,:)            = Pxx;
        MER_Freq(dLoop,chan,:)           = freq;

        
%         values=Pxx(xrange);
%         pos=xrange;
% %         params(1)=0.1;
%         Estimates=fminsearch(@errorfunc1dx,0.1,options,pos,values);
%         for j=1:3
%             params=Estimates;
%             Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
%         end
%         MER_FIT(dLoop, chan)=Estimates(1);
                
%         LFP_filt =  filtfilt(fb_LFP,fa_LFP,double(MER_data(dLoop).data(chan,:)));
%         [Pxx,www] = pwelch(LFP_filt,...
%                         hanning(ceil(Fs)),ceil(Fs*0.2),ceil(Fs), Fs);
%         LFPalpha(dLoop, chan)             = median(Pxx(6:12));
%         LFPbeta(dLoop, chan)              = median(Pxx(13:25));
%         LFPgamma1(dLoop, chan)            = median(Pxx(26:45));
%         LFPgamma2(dLoop, chan)            = median(Pxx(55:100));
%         LFPgamma3(dLoop, chan)            = median(Pxx(101:300));
                    
                    
    end
end
ft_progress('close')

MER_summary.SpcMed     = MER_SpcMed;
MER_summary.MedVal     = MER_MedVal;
MER_summary.RMSVal     = MER_RMSVal;
MER_summary.MER_Spc    = MER_Spc;
MER_summary.MER_Freq    = MER_Freq;
return