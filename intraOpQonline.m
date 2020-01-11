clear all;   clc;    close all;
% tic
ft_defaults

%% Import all data  (~8 sec for 24 depths)
startTime = cputime;
directory = 'J:\Surgery\382876498';
disp(['Importing data from directory: ', directory])
disp([' ' ])

[MER_import MER_times] = MER_data_extract(directory);

disp([' ' ])
disp(['Import complete!'])
disp([' ' ])
stopTime = cputime;
disp(['Took ', num2str(stopTime-startTime), ' seconds!'])

%% Extract correponding depths (~20 ms)
disp([' ' ])
disp(['Extracting depths and channel configuration... ' ])
A       = fileread('protokoll.txt');
B       = strfind(A, 'New site no');
stim    = strfind(A, 'Stimulation');
B       = B(find(B<stim(1)));

depthEx = [];

for loop = 1:length(B)
    char_tok = tokenize(A(B(loop)+28:B(loop)+40), ': ')
    if ~isempty(char_tok{3})
        depthEx(loop) = str2num(char_tok{3});
    else
        depthEx(loop) = -20;
    end
end
% depthTrue = [min(depthEx(find(depthEx~=-20))):1:-6, -5:0.5:max(depthEx)];


% Extract channel configuration

A              = fileread('protokoll.txt');
B              = strfind(A, 'channel configuration');
what           = {A(B(1)+23:B(1)+82)};
what_tok       = tokenize(char(what),':*,');
channel_config = {what_tok{2},what_tok{4},what_tok{6},what_tok{8},what_tok{10}};

clearvars A B stim char_tok what_tok what
disp([' ' ])
disp(['Channel Configuration: ', channel_config{:}])
disp([' ' ])
% takes about 8 secs for 24 depths

%% Choose depths- compare depths vs depthEx

h1 = figure(1);
subplot 211
% stem(depthTrue), hold on,
plot(depthEx, '-o')
title('Depths')

subplot 212
h1 = figure(1);
plot(MER_times, '-ob'), hold on
plot([1 length(MER_times)],[mean(MER_times) mean(MER_times)],'k')
plot([1 length(MER_times)],[mean(MER_times)/2-5 mean(MER_times)/2-5],'m')
plot([1 length(MER_times)],[mean(MER_times)/2+5 mean(MER_times)/2+5],'m')
str1 = num2str(int16(mean(MER_times)/2+5));
text(length(MER_times),mean(MER_times)/2+5,str1,'HorizontalAlignment', 'Right',...
                               'FontSize', 12, 'FontWeight', 'bold', ...
                                        'VerticalAlignment', 'bottom')

str2 = num2str(int16(mean(MER_times)/2-5));
text(length(MER_times),mean(MER_times)/2-5,str2,'HorizontalAlignment', 'Right',...
                                  'FontSize', 12, 'FontWeight', 'bold', ...
                                        'VerticalAlignment', 'bottom')


str3 = num2str(int16(mean(MER_times)));
text(length(MER_times),mean(MER_times),str3,'HorizontalAlignment', 'Right',...
                                'FontSize', 12, 'FontWeight', 'bold', ...
                                        'VerticalAlignment', 'bottom')

[val pos] = min(MER_times);
str4 = num2str(floor(val));
text(pos,val,str4,'HorizontalAlignment', 'Right',...
                                'FontSize', 12, 'FontWeight', 'bold', ...
                                        'VerticalAlignment', 'bottom')                                    
[val pos] = max(MER_times);
str5 = num2str(floor(val));
text(pos,val,str5,'HorizontalAlignment', 'Right',...
                                'FontSize', 12, 'FontWeight', 'bold', ...
                                        'VerticalAlignment', 'bottom')                                    
                                    
xlabel('Depths'), ylabel('Data time lengths (s)')
title('Length of recordings')

proceed =1;
while(proceed)
    disp([' ' ])
    reject = input('Choose depths to reject[1 2 3...]:');
    toiChoice  = input('Choose a toi [10 20]: ');
    disp([' ' ])
    
    MER_data = MER_import(setdiff(1:length(MER_import),reject));
    
    disp([' ' ])
    disp(['MER data selected. No. of depths: ', num2str(length(MER_data))])
    disp(['MER time of interest: ',...
        num2str(toiChoice(1)), ' to ', num2str(toiChoice(2)), ' seconds'])
    disp([' ' ])
    proceed = input('Do you wish to proceed? [y/n]: ', 's');
    if proceed == 'y'
        break
    else 
        continue
    end
end


dChoice   = depthEx(setdiff(1:length(MER_import),reject));
depths = [min(dChoice):1:-6, -5:0.5:max(dChoice)];

if ~(length(MER_data) == length(depths))
   disp('***   Error: Discrepancy - MER data and depths')    
   disp('***   Using values extracted from "Protocol.txt"')    
   depths = dChoice;
end

for loop = 1:length(MER_data)
   
    MER_data(loop).label = channel_config;
    MER_data(loop).depth = depths(loop);
    
end

close(h1)

%% Filter data, pwelch power spectrum, 

startTime = cputime;

[MER_summary] = MER_proc(MER_data, toiChoice);

disp([' ' ])
stopTime = cputime;
disp(['Done!' ])
disp(['Took ', num2str(stopTime-startTime), ' seconds!'])

%% Norm channels

valMat = MER_summary.SpcMed;
normVal.SpcMed = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);
valMat = MER_summary.RMSVal;
normVal.RMSVal = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);
valMat = MER_summary.MedVal;
normVal.MedVal = (valMat-repmat(mean(valMat),length(valMat),1))...
                    ./repmat(valMat(1,:),length(valMat),1);

valMat = MER_summary.SpcMed;
normVal.SpcMed = valMat./repmat(valMat(1,:),length(valMat),1);
valMat = MER_summary.RMSVal;
normVal.RMSVal = valMat./repmat(valMat(1,:),length(valMat),1);
valMat = MER_summary.MedVal;
normVal.MedVal = valMat./repmat(valMat(1,:),length(valMat),1);
                
% Plot MER through depths

disp([' ' ])
disp('###  Plots   ###')
disp([' ' ])

h1 = figure(1);
subplot 321
plot(depths,normVal.SpcMed, '-o', 'LineWidth',2,...
    'MarkerSize',6)
title('Median of spectral power')
legend(channel_config, 'Location', 'northwest')

subplot 323
plot(depths,normVal.MedVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('Median of MER')

subplot 325
plot(depths,normVal.RMSVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('RMS of MER')

% suptitle('Normalized plots')

% h1 = figure(2);
subplot 322
plot(depths,MER_summary.SpcMed, '-o', 'LineWidth',2,...
    'MarkerSize',6)
title('Median of spectral power')
legend(channel_config, 'Location', 'northwest')

subplot 324
plot(depths,MER_summary.MedVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('Median of MER')

subplot 326
plot(depths,MER_summary.RMSVal, '-o', 'LineWidth',2,...
                        'MarkerSize',6)
title('RMS of MER')

suptitle('Normalized plots   < - - >     No-norm')

%%
dUse = [1:16];
% dUse = [1:14 15:20];

figure(2)
foi = 500:3000;
for chan = 1:5
    subplot(1,5,chan)
    imagesc(squeeze(MER_summary.MER_Freq(1,chan,foi)),...
           depths(dUse),   squeeze(MER_summary.MER_Spc(dUse,chan,foi)));
    title(channel_config{chan})
    xlabel('Frequency in Hz')
    ylabel('Depth in mm')
end
suptitle('MER')

MER_Spc_norm  =  MER_summary.MER_Spc - ...
            repmat(MER_summary.MER_Spc(1,:,:),size(MER_summary.MER_Freq,1),1);
MER_Spc_norm  = (MER_summary.MER_Spc - min(min(min(MER_summary.MER_Spc))))...
                   ./max(max(max(MER_summary.MER_Spc)));

figure(3)
for chan = 1:5
    subplot(1,5,chan)
    imagesc(squeeze(MER_summary.MER_Freq(1,chan,foi)),...
           depths(dUse),   squeeze(MER_Spc_norm(dUse,chan,foi)));
    title(channel_config{chan})
    xlabel('Frequency in Hz')
    ylabel('Depth in mm')
end
suptitle('MER norm')






