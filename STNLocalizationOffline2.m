%%%%Online STN localization based on background activity (RMS and MA) of micro and macro-tip
%%%%2013-01-10 by Simeon Knieling

    clear all;      clc;        close all
    
    cd('D:\MATLAB-2012b\samples\BAM')
    options=optimset('Display','off', 'MaxIter', 300 , 'MaxFunEvals', 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 , 'TolFun', 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001, 'TolX', 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001);

    load depth;
    load microElectrode;
    load macroElectrode;
    load MicroMacro;
    
    depth = depth([1:12, 14:end]);
    microElectrode(1,1).depth = microElectrode(1,1).depth([1:12, 14:end]);
    microElectrode(1,2).depth = microElectrode(1,2).depth([1:12, 14:end]);
    microElectrode(1,1).MA = microElectrode(1,1).MA([1:12, 14:end]);
    microElectrode(1,2).MA = microElectrode(1,2).MA([1:12, 14:end]);
    microElectrode(1,1).RMS = microElectrode(1,1).RMS([1:12, 14:end]);
    microElectrode(1,2).depth = microElectrode(1,2).depth([1:12, 14:end]);
    
    macroElectrode(1,1).depth = macroElectrode(1,1).depth([1:12, 14:end]);
    macroElectrode(1,2).depth = macroElectrode(1,2).depth([1:12, 14:end]);
    macroElectrode(1,1).MA = macroElectrode(1,1).MA([1:12, 14:end]);
    macroElectrode(1,2).MA = macroElectrode(1,2).MA([1:12, 14:end]);
    macroElectrode(1,1).RMS = macroElectrode(1,1).RMS([1:12, 14:end]);
    macroElectrode(1,2).depth = macroElectrode(1,2).depth([1:12, 14:end]);
    
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
        for channel=1:channelCount
            [fb,fa]=ellip(2,0.1,40,[fmin_spikes fmax_spikes]*2/sr);
            hdata=filtfilt(fb,fa,double(depth(depthCount).RAW(channel).data));
            if(MicroMacro(channel).type==1)
               
                depth(depthCount).SpikesMicro(MicroMacro(channel).pos).data= hdata;
                microElectrode(MicroMacro(channel).pos).MA(depthCount)=median(abs(hdata));
                microElectrode(MicroMacro(channel).pos).STD(depthCount)=std(hdata);
                microElectrode(MicroMacro(channel).pos).RMS(depthCount)=sqrt(mean(hdata.^2));
                

                [Pxx,www] = pwelch(double(depth(depthCount).RAW(channel).data),hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
                microElectrode(MicroMacro(channel).pos).SPC(depthCount)=median(Pxx(500:1000));
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
    

            else
                
                depth(depthCount).SpikesMacro(MicroMacro(channel).pos).data= hdata;
                macroElectrode(MicroMacro(channel).pos).MA(depthCount)=median(abs(hdata));
                macroElectrode(MicroMacro(channel).pos).STD(depthCount)=std(hdata);
                macroElectrode(MicroMacro(channel).pos).RMS(depthCount)=sqrt(mean(hdata.^2));

                
                

                
                [Pxx,www] = pwelch(double(depth(depthCount).RAW(channel).data),hanning(ceil(sr)),ceil(sr*0.2),ceil(sr), sr);
                macroElectrode(MicroMacro(channel).pos).SPC(depthCount)=median(Pxx(500:1000));
                macroElectrode(MicroMacro(channel).pos).spec(depthCount,:)=Pxx(1:22000);
                length(Pxx)
                
                
                

                values=Pxx(xrange);
                pos=xrange;
                params(1)=0.1;
                Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
                for j=1:3
                    params(1)=Estimates(1);
                    Estimates=fminsearch(@errorfunc1dx,params,options,pos,values);
                end
                macroElectrode(MicroMacro(channel).pos).FIT(depthCount)=Estimates(1);
                
                
                
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

        
        h2 = figure(2)
        plot(-macroElectrode(1).FIT); hold on
        plot(macroElectrodeNorm(1).MA, 'b');
        plot(macroElectrodeNorm(1).RMS, '-or');
        plot(macroElectrodeNorm(1).STD, 'k');
        plot(macroElectrodeNorm(1).SPC, 'c');
        h2 = legend('Fit', 'MA', 'RMS', 'STD', 'SPC');
        title('Macro');
        
        save microElectrodeNorm microElectrodeNorm;
        save macroElectrodeNorm macroElectrodeNorm;
        save microElectrodeNoNorm microElectrodeNoNorm;
        save macroElectrodeNoNorm macroElectrodeNoNorm;
% 
%             microElectrodeNorm(i).MA;
%             microElectrodeNorm(i).RMS;
%             microElectrodeNorm(i).STD;
%             microElectrodeNorm(i).FIT;
%             microElectrodeNorm(i).SPC;
%         h=figure(71);
%         hold on;
%         clf();
% %         set(gcf,'Units','normalized');
% %         set(gcf,'OuterPosition',[0 0.5 0.5 0.5]);
%         set(gcf, 'PaperSize', [100 50/4*length(microElectrodeNorm)], 'PaperPosition', [0 0 100 50/4*length(microElectrodeNorm)]);
%         for i=1:length(microElectrodeNorm)
%             subplot(length(microElectrodeNorm),5,i*5-4);
%             hold on;
%             %plot_matrix(microElectrode(i).spec,microElectrodeNorm(i).depth(1:size(microElectrode(i).spec,1)),1:10000); 
%             %axis([min(microElectrodeNorm(i).depth) max(microElectrodeNorm(i).depth) 1 10000] );
%             xlabel('DTT');
%             title(['Micro ' num2str(i) ' - PSD']);
%             set(gca,'XTick',min(microElectrodeNorm(i).depth):1:max(microElectrodeNorm(i).depth))
% 
%             subplot(length(microElectrodeNorm),5,i*5-3);
%             hold on;
%             plot_matrix(macroElectrode(i).spec,macroElectrodeNorm(i).depth(1:size(macroElectrode(i).spec,1)),1:10000); 
%             axis([min(macroElectrodeNorm(i).depth) max(macroElectrodeNorm(i).depth) 1 10000] );
%             xlabel('DTT');
%             title(['Macro ' num2str(i) ' - PSD']);
%             set(gca,'XTick',min(macroElectrodeNorm(i).depth):1:max(macroElectrodeNorm(i).depth))
% 
%             subplot(length(microElectrodeNorm),5,i*5-2);
%             hold on;
%             mymatrix(:,1)=microElectrodeNorm(i).MA;
%             mymatrix(:,2)=microElectrodeNorm(i).RMS;
%             mymatrix(:,3)=microElectrodeNorm(i).STD;
%             mymatrix(:,4)=microElectrodeNorm(i).FIT;
%             mymatrix(:,5)=microElectrodeNorm(i).SPC;
%             plot_matrix(mymatrix, microElectrodeNorm(i).depth(1:size(microElectrode(i).spec,1)),1:6); 
%             axis([min(microElectrodeNorm(i).depth) max(microElectrodeNorm(i).depth) 1 6] );
%             xlabel('DTT');
%             title(['Macro ' num2str(i) ' - All Features']);
%             set(gca,'XTick',min(microElectrodeNorm(i).depth):1:max(microElectrodeNorm(i).depth))
%             
%             
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).RMS, 'g');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).STD, 'b');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).MA, 'r');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).SPK, 'k');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).FIT, 'm');
% % %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).SPC, 'c', 'LineWidth', 2);
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).RMS, '.g');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).STD, '.b');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).MA, '.r');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).SPK, '.k');
% %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).FIT, '.m');
% % %             plot(microElectrodeNorm(i).depth, microElectrodeNorm(i).SPC, '.c', 'MarkerSize', 15);
% %             axis([min(microElectrodeNorm(i).depth) max(microElectrodeNorm(i).depth) -0.01 1.01]);
% %             title(['Micro ' num2str(i) ' - Normalized over all Micro-Tips']);
% %             set(gca,'XTick',min(microElectrodeNorm(i).depth):1:max(microElectrodeNorm(i).depth))
%             
%             subplot(length(macroElectrodeNorm),5,i*5-1);
%             hold on;
%             
%             mymatrix(:,1)=macroElectrodeNorm(i).MA;
%             mymatrix(:,2)=macroElectrodeNorm(i).RMS;
%             mymatrix(:,3)=macroElectrodeNorm(i).STD;
%             mymatrix(:,4)=macroElectrodeNorm(i).FIT;
%             mymatrix(:,6)=macroElectrodeNorm(i).SPC;
%             plot_matrix(mymatrix, macroElectrodeNorm(i).depth(1:size(macroElectrode(i).spec,1)),1:6); 
%             axis([min(macroElectrodeNorm(i).depth) max(macroElectrodeNorm(i).depth) 1 6] );
%             xlabel('DTT');
%             title(['Macro ' num2str(i) ' - All Features']);
%             set(gca,'XTick',min(macroElectrodeNorm(i).depth):1:max(macroElectrodeNorm(i).depth))
%             
%             
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).RMS, 'g');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).STD, 'b');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).MA, 'r');
% % %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).SPK, 'k');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).FIT, 'm');
% % %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).SPC, 'c', 'LineWidth', 2);
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).RMS, '.g');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).STD, '.b');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).MA, '.r');
% % %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).SPK, '.k');
% %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).FIT, '.m');
% % %             plot(macroElectrodeNorm(i).depth, macroElectrodeNorm(i).SPC, '.c', 'MarkerSize', 15);
% %             axis([min(macroElectrodeNorm(i).depth) max(macroElectrodeNorm(i).depth) -0.01 1.01]);
% %             title(['Macro ' num2str(i) ' - Normalized over all Macro-Tips']);
% %             set(gca,'XTick',min(macroElectrodeNorm(i).depth):1:max(macroElectrodeNorm(i).depth))
%             
%             subplot(length(microElectrodeNorm),5,i*5);
%             hold on;
%             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).MA,macroElectrodeNorm(i).MA,12,2);
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'r');
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.r');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'r');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'r');
%             
%             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).RMS,macroElectrodeNorm(i).RMS,12,2);
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'g');
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.g');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'g');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'g');
%             
%             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).STD,macroElectrodeNorm(i).STD,12,2);
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'b');
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.b');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'b');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'b');
%             
%             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).FIT,macroElectrodeNorm(i).FIT,12,2);
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'm');
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.m');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'm');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'm');
%             
%             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).SPC,macroElectrodeNorm(i).SPC,12,2);
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'c');
%             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.c');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'c');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'c');
%             
% %             [XCF,lags,bounds] = crosscorr(microElectrodeNorm(i).SPK,macroElectrodeNorm(i).SPK,12,2);
% %             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), 'k');
% %             plot(lags(floor(length(lags)/2):end), XCF(floor(length(lags)/2):end), '.k');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(1) bounds(1)], 'k');
%             plot([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end))], [bounds(2) bounds(2)], 'k');
% 
%             set(gca,'XTick',min(lags(floor(length(lags)/2):end)):1:max(lags(floor(length(lags)/2):end)))
%             axis([min(lags(floor(length(lags)/2):end)) max(lags(floor(length(lags)/2):end)) -0.01 1.01]);
%         end
%         %mtit(h, 'STN Localization - Left Hemisphere - 28.01.2013', 'xoff',0.0001,'yoff',0.025);
%         %print(gcf,'-dpng', '2013-01-28-left-STN.png');%2013-01-25-right-STN
% 
% %         figure(72)
% %         clf();
% %         set(gcf,'Units','normalized');
% %         set(gcf,'OuterPosition',[0 0 0.5 0.5]);
% %         for i=1:length(macroElectrodeNorm)
% %             subplot(length(macroElectrodeNorm),1,i);
% % 
% %         end
% 
% % 
% %         figure(81)
% %         clf();
% %         set(gcf,'Units','normalized');
% %         set(gcf,'OuterPosition',[0.5 0.5 0.5 0.5]);
% %         for i=1:length(microElectrodeNoNorm)
% %             subplot(length(microElectrodeNoNorm),1,i);
% %             hold on;
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).RMS, 'g', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).STD, 'b', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).MA, 'r', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).SPK, 'k', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).FIT, 'm', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).SPC, 'c', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).RMS, '.g', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).STD, '.b', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).MA, '.r', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).SPK, '.k', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).FIT, '.m', 'LineWidth', 1);
% %             plot(microElectrodeNoNorm(i).depth, microElectrodeNoNorm(i).SPC, '.c', 'LineWidth', 1);
% %             title(['Micro ' num2str(i)]);
% %         end
% % 
% % 
% %         figure(82)
% %         clf();
% %         set(gcf,'Units','normalized');
% %         set(gcf,'OuterPosition',[0.5 0 0.5 0.5]);
% %         for i=1:length(macroElectrodeNoNorm)
% %             subplot(length(macroElectrodeNoNorm),1,i);
% %             hold on;
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).RMS, 'g', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).STD, 'b', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).MA, 'r', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).SPK, 'k', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).FIT, 'm', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).SPC, 'c', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).RMS, '.g', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).STD, '.b', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).MA, '.r', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).SPK, '.k', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).FIT, '.m', 'LineWidth', 1);
% %             plot(macroElectrodeNoNorm(i).depth, macroElectrodeNoNorm(i).SPC, '.c', 'LineWidth', 1);
% %             title(['Macro ' num2str(i)]);
% %         end
% %             
% %           
% 
%         