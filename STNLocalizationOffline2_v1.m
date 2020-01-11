%%%%Online STN localization based on background activity (RMS and MA) of micro and macro-tip
%%%%2013-01-10 by Simeon Knieling

clear all;      clc;        close all

cd('D:\MATLAB-2012b\samples\BAM')
options=optimset('Display','off', 'MaxIter', 300 , 'MaxFunEvals', 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 , 'TolFun', 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001, 'TolX', 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001);

load depth;
load microElectrode;
load macroElectrode;
load MicroMacro;

% % del_depths = [5, 11, 20] ;
% % retain_depths = setdiff(1:length(depth), del_depths);
% % depth = depth([retain_depths]);
% microElectrode(1,1).depth = microElectrode(1,1).depth([retain_depths]);
% % microElectrode(1,2).depth = microElectrode(1,2).depth([retain_depths]);
% microElectrode(1,1).MA = microElectrode(1,1).MA([retain_depths]);
% % microElectrode(1,2).MA = microElectrode(1,2).MA([retain_depths]);
% microElectrode(1,1).RMS = microElectrode(1,1).RMS([retain_depths]);
% % microElectrode(1,2).RMS = microElectrode(1,2).RMS([retain_depths]);
% 
% macroElectrode(1,1).depth = macroElectrode(1,1).depth([retain_depths]);
% % macroElectrode(1,2).depth = macroElectrode(1,2).depth([retain_depths]);
% macroElectrode(1,1).MA = macroElectrode(1,1).MA([retain_depths]);
% % macroElectrode(1,2).MA = macroElectrode(1,2).MA([retain_depths]);
% macroElectrode(1,1).RMS = macroElectrode(1,1).RMS([retain_depths]);
% % macroElectrode(1,2).RMS = macroElectrode(1,2).RMS([retain_depths]);

%%

channelCount=length(MicroMacro);
sr=44642.856597900390;


fmin_spikes = 300;
fmax_spikes = 3000;


start=0;
stop=9;
yamp1=-70;
yamp2=50;
handle.max_channel=length(MicroMacro);

handle.fmin_sort = 300;                					%high pass filter for sorting (default 300)
handle.fmax_sort = 3000;               					%low pass filter for sorting (default 3000)
handle.sr=sr;

xrange=[10:1000 2250:3000];

%x=round(10.^[log10(xstart):log10(xend)/xend:log10(xend)]);



depthCount=0;
for d=1:length(depth)
depthCount=depthCount+1;
for channel=1:4
    [fb,fa]=ellip(2,0.1,40,[fmin_spikes fmax_spikes]*2/sr);    % [300 3000]
    hdata=filtfilt(fb,fa,double(depth(depthCount).RAW(channel).data));
    if(MicroMacro(channel).type==1)
%         disp(MicroMacro(channel).type)
        depth(depthCount).SpikesMicro(MicroMacro(channel).pos).data= hdata;
        microElectrode(MicroMacro(channel).pos).MA(depthCount)=median(abs(hdata));
        microElectrode(MicroMacro(channel).pos).STD(depthCount)=std(hdata);
        microElectrode(MicroMacro(channel).pos).RMS(depthCount)=sqrt(mean(hdata.^2));

%         [Pxx, www] = pburg(double(depth(depthCount).RAW(channel).data),2,ceil(sr),ceil(sr));

        [Pxx,www] = pwelch(double(depth(depthCount).LFPMicro(channel).data),...
            hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
        microElectrode(MicroMacro(channel).pos).SPC(depthCount)=median(Pxx(300:3000));
        microElectrode(MicroMacro(channel).pos).spec(depthCount,:)=Pxx(1:22000);



        values=Pxx(xrange);
        pos=xrange;
        params(1)=0.1;
        Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        for j=1:3
            params(1)=Estimates(1);
            Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        end
        microElectrode(MicroMacro(channel).pos).FIT(depthCount)=Estimates(1);

        depth(depthCount).SpikesMicro(MicroMacro(channel).pos).data= hdata;
        microElectrode(MicroMacro(channel).pos).MA(depthCount)=median(abs(hdata));
        microElectrode(MicroMacro(channel).pos).STD(depthCount)=std(hdata);
        microElectrode(MicroMacro(channel).pos).RMS(depthCount)=sqrt(mean(hdata.^2));

        [Pxx, www] = pburg(double(depth(depthCount).RAW(channel).data),2,ceil(sr),ceil(sr));
%         [Pxx,www] = pwelch(double(depth(depthCount).RAW(channel).data),hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
        microElectrode(MicroMacro(channel).pos).SPC_pburg(depthCount)=median(Pxx(300:3000));
        microElectrode(MicroMacro(channel).pos).spec_pburg(depthCount,:)=Pxx(1:22000);



        values=Pxx(xrange);
        pos=xrange;
        params(1)=0.1;
        Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        for j=1:3
            params(1)=Estimates(1);
            Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        end
        microElectrode(MicroMacro(channel).pos).FIT_pburg(depthCount)=Estimates(1);        
        
        

    else
%     disp(MicroMacro(channel).type)
        depth(depthCount).SpikesMacro(MicroMacro(channel).pos).data= hdata;
        macroElectrode(MicroMacro(channel).pos).MA(depthCount)=median(abs(hdata));
        macroElectrode(MicroMacro(channel).pos).STD(depthCount)=std(hdata);
        macroElectrode(MicroMacro(channel).pos).RMS(depthCount)=sqrt(mean(hdata.^2));





        [Pxx,www] = pwelch(double(depth(depthCount).RAW(channel).data),hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
%         [Pxx,www] = pburg(double(depth(depthCount).RAW(channel).data), ceil(sr),'ConfidenceLevel',0.99);
%         [Pxx, www] = pburg(double(depth(depthCount).RAW(channel).data),2,ceil(sr),ceil(sr));

        macroElectrode(MicroMacro(channel).pos).SPC(depthCount)=median(Pxx(500:2500));
        macroElectrode(MicroMacro(channel).pos).spec(depthCount,:)=Pxx(1:22000);
        length(Pxx);




        values=Pxx(xrange);
        pos=xrange;
        params(1)=0.1;
        Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        for j=1:3
            params(1)=Estimates(1);
            Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        end
        macroElectrode(MicroMacro(channel).pos).FIT(depthCount)=Estimates(1);

%           [Pxx,www] = pwelch(double(depth(depthCount).RAW(channel).data),hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
%         [Pxx,www] = pburg(double(depth(depthCount).RAW(channel).data), ceil(sr),'ConfidenceLevel',0.99);
        [Pxx, www] = pburg(double(depth(depthCount).RAW(channel).data),2,ceil(sr),ceil(sr));

        macroElectrode(MicroMacro(channel).pos).SPC_pburg(depthCount)=median(Pxx(500:2500));
        macroElectrode(MicroMacro(channel).pos).spec_pburg(depthCount,:)=Pxx(1:22000);
        length(Pxx);




        values=Pxx(xrange);
        pos=xrange;
        params(1)=0.1;
        Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        for j=1:3
            params(1)=Estimates(1);
            Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
        end
        macroElectrode(MicroMacro(channel).pos).FIT_pburg(depthCount)=Estimates(1);



    end
end
end
%%%%%%%Select from here to regenerate plots
microElectrodeNorm=microElectrode;
macroElectrodeNorm=macroElectrode;
microElectrodeNoNorm=microElectrode;
macroElectrodeNoNorm=macroElectrode;
count=1;
for i=1:length(microElectrode)
    for j=1:length(microElectrode(i).depth)
        microMA(count)=microElectrode(i).MA(j);
        microSTD(count)=microElectrode(i).STD(j);
        microRMS(count)=microElectrode(i).RMS(j);
        microFIT_pburg(count)=microElectrode(i).FIT_pburg(j);
        microSPC_pburg(count)=microElectrode(i).SPC_pburg(j);
        microFIT(count)=microElectrode(i).FIT(j);
        microSPC(count)=microElectrode(i).SPC(j);
        count=count+1;
    end
end
count=1;
for i=1:length(macroElectrode)
    for j=1:length(macroElectrode(i).depth)
        macroMA(count)=macroElectrode(i).MA(j);
        macroSTD(count)=macroElectrode(i).STD(j);
        macroRMS(count)=macroElectrode(i).RMS(j);
        macroFIT(count)=macroElectrode(i).FIT(j);
        macroSPC(count)=macroElectrode(i).SPC(j);
        macroFIT_pburg(count)=macroElectrode(i).FIT_pburg(j);
        macroSPC_pburg(count)=macroElectrode(i).SPC_pburg(j);
        count=count+1;
    end
end
for i=1:length(microElectrodeNorm)
    microElectrodeNorm(i).MA=(microElectrodeNorm(i).MA-min(microMA))/(max(microMA)-min(microMA));
    microElectrodeNorm(i).STD=(microElectrodeNorm(i).STD-min(microSTD))/(max(microSTD)-min(microSTD));
    microElectrodeNorm(i).RMS=(microElectrodeNorm(i).RMS-min(microRMS))/(max(microRMS)-min(microRMS));
    microElectrodeNorm(i).FIT=(microElectrodeNorm(i).FIT-min(microFIT))/(max(microFIT)-min(microFIT));
    microElectrodeNorm(i).SPC=(microElectrodeNorm(i).SPC-min(microSPC))/(max(microSPC)-min(microSPC));
end
for i=1:length(macroElectrodeNorm)
    macroElectrodeNorm(i).MA=(macroElectrodeNorm(i).MA-min(macroMA))/(max(macroMA)-min(macroMA));
    macroElectrodeNorm(i).STD=(macroElectrodeNorm(i).STD-min(macroSTD))/(max(macroSTD)-min(macroSTD));
    macroElectrodeNorm(i).RMS=(macroElectrodeNorm(i).RMS-min(macroRMS))/(max(macroRMS)-min(macroRMS));
    macroElectrodeNorm(i).FIT=(macroElectrodeNorm(i).FIT-min(macroFIT))/(max(macroFIT)-min(macroFIT));
    macroElectrodeNorm(i).SPC=(macroElectrodeNorm(i).SPC-min(macroSPC))/(max(macroSPC)-min(macroSPC));
end
%%%%%Even out RMS-MA for not normalized windows
for i=1:length(microElectrodeNoNorm)
    maxMA=max(microElectrodeNoNorm(i).MA);
    maxSTD=max(microElectrodeNoNorm(i).STD);
    maxRMS=max(microElectrodeNoNorm(i).RMS);
    maxFIT=max(microElectrodeNoNorm(i).FIT);
    maxSPC=max(microElectrodeNoNorm(i).SPC);
    div=maxRMS/maxMA;
    microElectrodeNoNorm(i).RMS=microElectrodeNoNorm(i).RMS/div;
    div=maxSTD/maxMA;
    microElectrodeNoNorm(i).STD=microElectrodeNoNorm(i).STD/div;
    div=maxFIT/maxMA;
    microElectrodeNoNorm(i).FIT=microElectrodeNoNorm(i).FIT/div;
    div=maxSPC/maxMA;
    microElectrodeNoNorm(i).SPC=microElectrodeNoNorm(i).SPC/div;
end
for i=1:length(macroElectrodeNoNorm)
    maxMA=max(macroElectrodeNoNorm(i).MA);
    maxSTD=max(macroElectrodeNoNorm(i).STD);
    maxRMS=max(macroElectrodeNoNorm(i).RMS);
    maxFIT=max(macroElectrodeNoNorm(i).FIT);
    maxSPC=max(macroElectrodeNoNorm(i).SPC);
    div=maxRMS/maxMA;
    macroElectrodeNoNorm(i).RMS=macroElectrodeNoNorm(i).RMS/div;
    div=maxSTD/maxMA;
    macroElectrodeNoNorm(i).STD=macroElectrodeNoNorm(i).STD/div;
    div=maxFIT/maxMA;
    macroElectrodeNoNorm(i).FIT=macroElectrodeNoNorm(i).FIT/div;
    div=maxSPC/maxMA;
    macroElectrodeNoNorm(i).SPC=macroElectrodeNoNorm(i).SPC/div;
end

%%
for k=1:length(microElectrode)
    electrode(k).depth=microElectrode(k).depth;
end

h1 = figure(1); 
plot(-microElectrode(1).FIT); hold on
plot(microElectrodeNorm(1).MA, 'b');
plot(microElectrodeNorm(1).RMS, '-or');
plot(microElectrodeNorm(1).STD, 'k');
plot(microElectrodeNorm(1).SPC, 'c');
h1 = legend('Fit', 'MA', 'RMS', 'STD', 'SPC');
title('Micro');

%%

h2 = figure(2);
plot(electrode(1).depth,-macroElectrode(1).FIT); hold on
plot(electrode(1).depth,macroElectrodeNorm(1).MA, 'b');
plot(electrode(1).depth,macroElectrodeNorm(1).RMS, '-or');
plot(electrode(1).depth,macroElectrodeNorm(1).STD, 'k');
plot(electrode(1).depth,macroElectrodeNorm(1).SPC, 'c');
h2 = legend('Fit', 'MA', 'RMS', 'STD', 'SPC');
title('Macro');

%%
save microElectrodeNorm microElectrodeNorm;
save macroElectrodeNorm macroElectrodeNorm;
save microElectrodeNoNorm microElectrodeNoNorm;
save macroElectrodeNoNorm macroElectrodeNoNorm;

%% 

% load microElectrodeNorm
% load macroElectrodeNorm
% load microElectrodeNoNorm
% load macroElectrodeNoNorm

%%
for elecLoop = 1:length(macroElectrodeNorm)
    features_micro(elecLoop,:,:) = [macroElectrode(elecLoop).MA; ...
                                    macroElectrode(elecLoop).RMS;...
                                    macroElectrode(elecLoop).STD; 
                                    macroElectrode(elecLoop).SPC; 
                                    -macroElectrode(elecLoop).FIT];
end

%% PCA
% feature = squeeze(features_micro(1,:,:))';
% feature = feature - mean(mean(feature));
% 
% DataCov=cov(feature); % covariance matrix
% [PC,variance,explained] = pcacov(DataCov);

%%
labels = {'MA', 'RMS', 'STD', 'SPC', 'FIT'};
feature = squeeze(features_micro(1,:,:))';
[w, pc, ev, tsquare] = princomp(feature, 'econ');

mu = mean(feature);
F_hat = bsxfun(@minus,feature,mu); %# subtract the mean
norm(pc * w' - F_hat)

data_PCA = pc(:,1) * w(:,1)';

data_PCA_norm = (data_PCA-min(min(data_PCA)))./ (max(max(data_PCA))-min(min(data_PCA)));
feature_norm  = (feature-min(min(feature)))   ./ (max(max(feature))- min(min(feature)));

figure
plot(electrode(1).depth, data_PCA_norm), title('After PCA')
hold on
plot(electrode(1).depth, [0.5])

figure
plot(electrode(1).depth, feature_norm), title('Before PCA')


cumsum(ev)./sum(ev)
figure
biplot(w(:,1:3),'Scores',pc(:,1:3), 'varlabels',  labels), title('PCA')

%%

for j = 1:5
   data_PCA_norm(:,j) = (data_PCA(:,j) - min(min(data_PCA(:,j))))...
        ./(max(max(data_PCA(:,j)))- min(min(data_PCA(:,j))));
end

%%
% plot(cumsum(mean(data_PCA_norm,2)))

%%
% figure, plot(electrode(1).depth,feature(:,1:2), '-ob')
% hold on, plot(electrode(1).depth,data_PCA_norm(:,1:2), '-or')
% hold on, plot(electrode(1).depth,mean(data_PCA_norm(:,1:2),2), '-ok')

% figure
% plot(cumsum(mean(feature(:,3:5),2)))
% figure
% plot(cumsum(std(data_PCA')), 'r')
% figure
% plot(cumsum(var(data_PCA')), 'k')
% hold on 

%%

% for depth_loop = 1:length(depth)
%     data_Micro(:,depth_loop) = (depth(depth_loop).LFPMicro.data);
%     data_Macro(:,depth_loop) = (depth(depth_loop).LFPMacro.data);
% end
% 
% figure
% boxplot(data_Micro), title('Micro');
% figure
% boxplot(data_Macro),  title('Macro');

%%

% figure
% plot((macroFIT - min(macroFIT))./ (max(macroFIT)-min(macroFIT))), hold on
% plot((macroFIT_pburg - min(macroFIT_pburg))./ ...
%             (max(macroFIT_pburg)-min(macroFIT_pburg)), 'r')


[S,F,T,P] = spectrogram(depth(1).LFPMicro(1).data,ceil(sr),ceil(sr/8),ceil(sr),sr);
surf(T,F(1:5000),10*log10(P(1:5000,:)),'edgecolor','none'); axis tight; 
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');

% Concatenate all data at all depths

for loop = 1:length(depth)
   
    len_mat_micro(loop) = length(depth(loop).LFPMicro(1).data);
    len_mat_macro(loop) = length(depth(loop).LFPMacro(1).data);
    
end

max_len_mat_macro = max(len_mat_macro);
min_len_mat_macro = min(len_mat_macro);

len = min_len_mat_macro;
start = 1;
for loop = 1:length(depth)
%     len = min_len_mat_macro;
    conc_LFP_micro(start:(start+len)-1) = ...
        depth(loop).LFPMicro(1).data(1:len);
    conc_LFP_macro(start:(start+len)-1) = ...
        depth(loop).LFPMacro(1).data(1:len);
    start = start + len;

end

%%
[S,F,T,P] = spectrogram(conc_LFP_macro,ceil(sr),ceil(sr/8),ceil(sr),sr);
surf(T,F(1:5000),10*log10(P(1:5000,:)),'edgecolor','none'); axis tight; 
view(0,90); title('Macro')
xlabel('Time (Seconds)'); ylabel('Hz');

[S,F,T,P] = spectrogram(conc_LFP_micro,ceil(sr),ceil(sr/8),ceil(sr),sr);
surf(T,F(1:5000),10*log10(P(1:5000,:)),'edgecolor','none'); axis tight; 
view(0,90); title('Micro')
xlabel('Time (Seconds)'); ylabel('Hz');

%%
figure
bar(electrode(1).depth,data_PCA);
xlabel('Feature'), ylabel('Depth in mm');
title('Feature map after PCA')

figure
bar(electrode(1).depth,feature);
xlabel('Feature'), ylabel('Depth in mm');
title('Feature map before PCA')
