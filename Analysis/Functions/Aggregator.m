function[] = Aggregator(tracks,filename)

close all

stimframe = 7;
cutofftime = 30;
stimulation = 7;

figure()
hold on
avgpeakdist = [];
stdlist =[];
highlist = [];
cdistance = [];
firstpeakheight = [];
firstpeakdist =[];

valuedata = NaN(length(tracks),100);

for i = 1:length(tracks)
    if length(tracks{i,1}) > cutofftime
        tracks{i} = tracks{i,1}(1:cutofftime,:);
    end
    
    tracks{i}(:,7) = tsmovavg(tracks{i,1}(:,6)','s',4)';
    tracks{i}(:,8) = (tracks{i,1}(:,7)-min(tracks{i,1}(:,7)))./(max(tracks{i,1}(:,7))-min(tracks{i,1}(:,7)));

    if max(tracks{i,1}(:,7)) > 7000
        highlist = [highlist; i];
    end
    [pks, locs] = findpeaks(tracks{i,1}(:,8),5*tracks{i,1}(:,2));
    tempfirstpeak = locs(find(locs > stimulation));
    cdistance = [cdistance; trapz(tracks{i,1}(tempfirstpeak(1)/5-4:locs(1)/5+4,7)),tracks{i,2}];
    firstpeakheight = [firstpeakheight; tracks{i,1}(tempfirstpeak(1)/5,7),tracks{i,2}];
    firstpeakdist = [firstpeakdist; tempfirstpeak(1)-stimulation,tracks{i,2}];
    peakInterval = diff(tempfirstpeak);
    stdlist = [stdlist;std(peakInterval),tracks{i,2}];
    avgpeakdist = [avgpeakdist; mean(diff(tempfirstpeak)),tracks{i,2}];
    
    for timeofx = 1:length(tracks{i}(:,2))
        valuedata(i,tracks{i}(timeofx,2)) = tracks{i}(timeofx,7);
    end

    
    if isequal(tracks{i,2},1)        
        g1= plot(5*tracks{i}(:,2), tracks{i}(:,7),'r');
        g1.Color(4) = 0.1;
        text(5*tracks{i}(end,2),tracks{i}(end,7),num2str(i))
    end
    
    if isequal(tracks{i,2},2)        
%         g2 = plot3(5*tracks{i}(:,2),i*ones(length(tracks{i}(:,7))), tracks{i}(:,7));
        %g2.Color(4) = 0.1;
        %text(5*tracks{i}(end,2),tracks{i}(end,7),num2str(i))  
    end

end

g1reporteravg = [];
g2reporteravg = [];
for j = 1:cutofftime
    g1tempavg = [];
    g2tempavg = [];
    for i = 1:length(tracks)
        if isequal(tracks{i,2},1) 
            g1tempavg = [g1tempavg; tracks{i,1}(tracks{i,1}(:,2)==j,7)];
        end
        if isequal(tracks{i,2},2) 
            g2tempavg = [g2tempavg; tracks{i,1}(tracks{i,1}(:,2)==j,7)];
        end
    end
    g1reporteravg = [g1reporteravg; [mean(g1tempavg(:))]];
    g2reporteravg = [g2reporteravg; [mean(g2tempavg(:))]];
end
    alpha(.2)
plot(1:5:cutofftime*5,g1reporteravg,'r')
plot(1:5:cutofftime*5,g2reporteravg,'b')
title('Raw Cell Nuclear Fluorescence Traces')
xlabel('time(mins)')
ylabel('Raw Fluorescence')

saveas(gcf,strcat('singlecelltraces_',filename,'.png'))

countersinglecell = 1;

 for i = 1:ceil(length(tracks)/25)
    figure()
    for ii = 1:25
        if countersinglecell == length(tracks)+1
            break
        end
        subplot(5,5,ii);
        if isequal(tracks(countersinglecell,2),1)        
            plot(5*tracks{countersinglecell,1}(:,2), tracks{countersinglecel,1}(:,7),'r')
        end
        if isequal(tracks(countersinglecell,2),2)        
            plot(5*tracks{countersinglecell,1}(:,2), tracks{countersinglecell,1}(:,7),'b')
        end

        plot(5*tracks{countersinglecell,1}(:,2), tracks{countersinglecell,1}(:,7))
        text(5*tracks{countersinglecell,1}(end,2),tracks{countersinglecell,1}(end,7),num2str(countersinglecell))
        %bla
        countersinglecell =  countersinglecell + 1;

    end
    saveas(gcf,strcat('singlecells_',filename,'_', num2str(i), '.png'))
 end
title('Raw Cell Nuclear Fluorescence Traces')
xlabel('time(mins)')
ylabel('Raw Fluorescence')




g1maxprob = [];
g2maxprob = [];


g1reporteravg = [];
g2reporteravg = [];
for j = 1:cutofftime
    g1tempavg = [];
    g2tempavg = [];
    for i = 1:length(tracks)
        if isequal(tracks{i,2},1) 
            g1tempavg = [g1tempavg; tracks{i}(tracks{i,1}(:,2)==j,8)];
        end
        if isequal(tracks{i,2},2) 
            g2tempavg = [g2tempavg; tracks{i}(tracks{i,1}(:,2)==j,8)];
        end
    end
    g1reporteravg = [g1reporteravg; [mean(g1tempavg(:))]];
    g2reporteravg = [g2reporteravg; [mean(g2tempavg(:))]];
    g1maxprob = [g1maxprob; [mean(g1tempavg(:))/max(g1tempavg(:))]];
    g2maxprob = [g2maxprob; [mean(g2tempavg(:))/max(g2tempavg(:))]];

end
figure()
hold on
plot(1:5:cutofftime*5,g1reporteravg,'r')
plot(1:5:cutofftime*5,g2reporteravg,'b')

title('Average Cell Nuclear Fluorescence Traces')
xlabel('time(mins)')
ylabel('Relative Fluorescence')
saveas(gcf,strcat('Average_singlecelltraces_',filename,'.png'))

figure()
hold on
plot(1:5:cutofftime*5,g1maxprob,'r')
plot(1:5:cutofftime*5,g2maxprob,'b')

title('Probability Cell Nuclear Fluorescence Traces')
xlabel('time(mins)')
ylabel('Relative Fluorescence')
saveas(gcf,strcat('ProbAvg_singlecelltraces_',filename,'.png'))

figure()
hold on
n = find(avgpeakdist(:,2)==1);
nn = find(avgpeakdist(:,2)==2);
histogram(avgpeakdist(n,1),15,'facecolor','r')
histogram(avgpeakdist(nn,1),15,'facecolor','b')
alpha(.3)
title('Average Peak Distance')
ylabel('Counts')
saveas(gcf,strcat('Peakdist_',filename,'.png'))


figure()
hold on
n = find(stdlist(:,2)==1);
nn = find(stdlist(:,2)==2);
histogram(stdlist(n,1),15,'facecolor','r')
histogram(stdlist(nn,1),15,'facecolor','b')
alpha(.3)
title('Standard Deviation of Peaks')
ylabel('Counts')
saveas(gcf,strcat('STDPeakdist_',filename,'.png'))

figure()
hold on
n = find(cdistance(:,2)==1);
nn = find(cdistance(:,2)==2);
histogram(cdistance(n,1),15,'facecolor','r')
histogram(cdistance(nn,1),15,'facecolor','b')
alpha(.3)
title('Integrated Area of First Peak')
ylabel('Counts')
saveas(gcf,strcat('FirstPeakArea_',filename,'.png'))

figure()
hold on
n = find(firstpeakheight(:,2)==1);
nn = find(firstpeakheight(:,2)==2);
histogram(firstpeakheight(n,1),15,'facecolor','r')
histogram(firstpeakheight(nn,1),15,'facecolor','b')
alpha(.3)
title('Height of First Peak')
ylabel('Counts')
saveas(gcf,strcat('FirstPeakHeight_',filename,'.png'))

figure()
hold on
n = find(firstpeakdist(:,2)==1);
nn = find(firstpeakdist(:,2)==2);
histogram(firstpeakdist(n,1),15,'facecolor','r')
histogram(firstpeakdist(nn,1),15,'facecolor','b')
alpha(.3)
title('Time to First Peak')
ylabel('Counts')
saveas(gcf,strcat('FirstPeakDist_',filename,'.png'))

end
