close all;
clc

% DATA = {NREMSoundallTrialsControl*100, NREMSoundallTrials*100};
DATA = {NREMSoundallTrials*100};

s2n = [];
incLatency = [];

Fs = round(1/(eventT(2)-eventT(1)));


for i = 1:length(DATA)
    data = DATA{i};
   
    w = size(data,1);
    for j = 2:2:w
        %initialize
        s2n(i, j/2) = nan;s2n(2+i, j/2) = nan;
        incLatency(i, j/2) = nan;incLatency(2+i, j/2) = nan;
        
        startP = 1600;
        eventT = eventT-eventT(1600);
        
        
        rtemp = smoothdata(data(j-1,:), 'movmean',20)';
        gtemp = smoothdata(data(j,:), 'movmean',20)';
        
        rbaseline = rtemp(startP-5*Fs+1:startP);
        gbaseline = gtemp(startP-5*Fs+1:startP);
        
        rtemp = rtemp - mean(rbaseline);
        gtemp = gtemp - mean(gbaseline);
        rbaseline = rtemp(1:1000);
        gbaseline = gtemp(1:1000);
        
        diffR = [0; diff(rtemp)./diff(eventT)];
        diffG = [0; diff(gtemp)./diff(eventT)];
        diffdiffR = [0; [0; diff(rtemp,2)]./diff(eventT)];
        diffdiffG = [0; [0; diff(gtemp,2)]./diff(eventT)];
        waketime = size(eventT,1)/2;
        
        f = figure; f.Position = [400, 500, 900, 350];
        subplot(1,2,1);hold on;
%         figure(i); hold on; 
        plot(eventT,rtemp, '-r');
%         plot(eventT,diffR, '-');
%         plot(zscore(smoothdata(diffdiffR, 'movmean',20)),'-k');
%         plot(diffG, '-');
%         plot(gtemp, '-g');
        ylim([-1.2 1.5]);
        xlim([-5 10]);


        tmean = mean(diffR(startP-5*Fs+1:startP)); tdev = std(diffR(startP-5*Fs+1:startP+5*Fs));
%         figure;
        [iupper, ilower] = cusum(diffR(1601:end),5,1,tmean, tdev,'all');
%         [~,~] = cusum(diffdiffR);
        [~, mini] = min(rtemp(1601:2600));
%         lowi = ilower(1);
%         mini = max(mini, lowi);
%         plot(eventT(mini+1600), rtemp(mini+1600), 'rs', 'MarkerFaceColor', 'k');
        itargetR = iupper(iupper>mini);
        
        tmean = mean(diffdiffR(startP+1:startP+5*Fs)); tdev = std(diffdiffR(startP+1:startP+5*Fs));
%         figure;cusum(diffdiffR(1601:end),10,1,tmean, tdev,'all')
        [iupper2, ilower2] = cusum(diffdiffR(1601:end),7,1,tmean, tdev,'all');
        Lia = ismember(itargetR, iupper2);
        if ~isempty(itargetR(Lia))
            itargetR = itargetR(Lia);
        end
        figure(f);
        subplot(1,2,1);hold on;
        plot(eventT(itargetR(1)+1600), rtemp(itargetR(1)+1600), 'rs', 'MarkerFaceColor', 'r');
%         legend('upper sum', 'lower sum');
%         figure; cusum(diffR(1601:end),5,1,tmean, tdev)
%         
        
        figure(f); 
        subplot(1,2,2);hold on;
%         figure(i); hold on; 
        plot(eventT,gtemp, '-g');
%         plot(eventT,diffG, '-');
%         plot(zscore(smoothdata(diffdiffG, 'movmean',20)),'-k');
        ylim([-1.2 1.5]);
        xlim([-5 10]);
        
        tmean = mean(diffG(startP-5*Fs+1:startP)); tdev = std(diffG(startP-5*Fs+1:startP+5*Fs));
%         tmean = mean(diffG(1:1600)); tdev = std(diffG(1:1600));
%         tmean = mean(diffG(1601:2601)); tdev = std(diffG(1601:2601));
%         figure;
        [iupper, ilower] = cusum(diffG(1601:end),5,1,tmean, tdev,'all');
%         [~,~] = cusum(diffdiffG);
        [~, mini] = min(gtemp(1601:2600));
%         lowi = ilower(1);
%         mini = max(mini, lowi);
        
%         plot(eventT(mini+1600), gtemp(mini+1600), 'rs', 'MarkerFaceColor', 'k');
        itargetG = iupper(iupper>mini);
        
        tmean = mean(diffdiffG(startP+1:startP+5*Fs)); tdev = std(diffdiffG(startP+1:startP+5*Fs));
%         figure;
        [iupper2, ilower2] = cusum(diffdiffG(1601:end),7,1,tmean, tdev,'all');
        Lia = ismember(itargetG, iupper2);
        if ~isempty(itargetG(Lia))
            itargetG = itargetG(Lia);
        end
        figure(f);
        subplot(1,2,2);hold on;
        plot(eventT(itargetG(1)+1600), gtemp(itargetG(1)+1600), 'rs', 'MarkerFaceColor', 'r');
%         figure; cusum(diffG(1601:end),5,1,tmean, tdev)


        incLatency(i+1,j/2) = eventT(itargetR(1)+startP);
        incLatency(i+3,j/2) = eventT(itargetG(1)+startP);
        
        rsd = std(rbaseline); % standard deviation of baseline period
        gsd = std(gbaseline);
        
        signalR = rtemp(startP+1:3000);
        signalG = gtemp(startP+1:3000);
        signalT = eventT(startP+1:3000);
        
        incR = [diff(signalR)>0; 0];
        incG = [diff(signalG)>0; 0];        
        incSigR = signalR.*incR;
        incSigG = signalG.*incG;
        
        
%         rthreshold = 3*rsd;
%         gthreshold = 3*gsd;
        rthreshold = 2*rsd;
        gthreshold = 2*gsd;
        
        ris = islocalmin(signalR); ris = find(ris==1); % find latency
        gis = islocalmin(signalG); gis = find(gis==1);
        
        highR = find(abs(signalR) >rthreshold);
        highG = find(abs(signalG) >gthreshold);
        
        l = 19;
%         l = 1;
        
        if length(highR)>l
%         if length(ris)>l
            jj = 1;
%             for jj = 1:length(ris)
%                 tmpris = ris(jj)+startP;
%                 if sum(diffR(tmpris:tmpris+l)>0) >= (l+1)*95/100
%                     break;
%                 end
%             end
            is1 = ris(jj) + startP;
            
            % calculate signal to noise ratio. only accept snr>2.5
            tempSNR = snr(rtemp(startP-5*Fs+1:startP+10*Fs),Fs);
            s2n(i,j/2) = tempSNR;
            
            if tempSNR>2.5
                incLatency(i,j/2) = eventT(itargetR(1)+startP);
                
            end
        end
        
        if length(highG)>l
%         if length(gis)>l
            jj = 1;
%             for jj = 1:length(gis)
%                 tmpgis = gis(jj)+startP;
%                 if sum(diffR(tmpgis:tmpgis+l)>0) >= (l+1)*95/100
%                     break;
%                 end
%             end
            is1 = gis(jj) + startP;
            
            tempSNR = snr(gtemp(startP-5*Fs+1:startP+10*Fs),Fs);
            s2n(2+i,j/2) = tempSNR;
            
            incLatency(2+i, j/2) = eventT(itargetG(1)+startP);
            if tempSNR>2.5                
                incLatency(2+i, j/2) = eventT(itargetG(1)+startP);
                
            end
        end

        
    end
    
    
end
clc