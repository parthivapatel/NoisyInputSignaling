function[maxlist] = graphs(tracks)

maxcell = max(tracks(:,111));
cellvals = [];
fn = [];

hold on
sumvals = zeros(18);
countvals = zeros(18);
integrated = [];
maxlist = [];
peaktime = [];

for i = 1:maxcell
    tempindex = tracks(find(tracks(:,111)==i),:);
    %/max(tempindex(:,3))
    tempmovmean = movmean(tempindex(:,3)-tempindex(4,3),3);
    plot(tempindex(:,110)*5-20,tempmovmean,'Color', [.3 .3 .3 .3])
    for j = 1:length(tempindex(:,110))
        sumvals(tempindex(j,110)) = sumvals(tempindex(j,110)) + tempmovmean(j);
        countvals(tempindex(j,110)) = countvals(tempindex(j,110)) +1;
    end
    
    integrated = [integrated; trapz(5*linspace(1,length(tempmovmean),length(tempmovmean)),tempmovmean/6500)];
    maxlist = [maxlist; max(tempmovmean)];
    peaktime = [peaktime; 5*find(tempmovmean == max(tempmovmean))];
    
    %,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    %'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05
%     active = input('active?');
%     cellvals = [cellvals; active];
%     
%     frame = input('fn?');
%     fn = [fn; frame];
end
plot(-15:5:70, sumvals./countvals, 'Color', [1 0 0])
out = horzcat(cellvals,fn);

