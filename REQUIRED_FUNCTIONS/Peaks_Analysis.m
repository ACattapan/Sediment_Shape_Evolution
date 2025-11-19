clear Data N_Data edges_Data bins_Data pks_Data locs_Data pks_Data_Val x_peaks y_peaks
BinWidth=0.02;
Fig_Peaks=figure
x_peaks=[];
y_peaks=[];
for i=1:size(LocPath{PathNum,1},2)

    Data{PathNum,i}=T.NormCirc(T.Location==LocPath{PathNum,1}(i)&T.Lithology=="ARENITES"&T.Origin==OrigPath{PathNum,1}(i));

    if length(Data{PathNum,i})>0
        [N_Data{PathNum,i}(1,:),edges_Data{PathNum,i}(1,:)] = histcounts(Data{PathNum,i},'BinWidth',BinWidth);%'BinLimits',[0,1],'BinMethod','auto',
        bins_Data{PathNum,i}=conv(edges_Data{PathNum,i},[0.5,0.5],'valid');
        if length(N_Data{PathNum,i}(1,:))>=3
            [pks_Data{PathNum,i},locs_Data{PathNum,i}] = findpeaks(N_Data{PathNum,i}(1,:));
            pks_Data_Val{PathNum,i}=bins_Data{PathNum,i}(locs_Data{PathNum,i});
            for j=1:length(pks_Data_Val{PathNum,i})
                x_peaks=[x_peaks X{PathNum,1}.SourceDist(X{PathNum,1}.Location==LocPath{PathNum,1}(i))];%LocDist(LocPath{PathNum,1}(i),1)
                y_peaks=[y_peaks pks_Data_Val{PathNum,i}(j)];
            end
        end
    end

end

scatter(x_peaks,y_peaks, 'filled')
ylim([0.7 1])
xlabel('Distance from sediment source [km]');
ylabel('Peaks of Normalized Isoperimetric Ratio distribution [-]');
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Loc18_AR=T.NormCirc(T.Location==18&T.Lithology=="ARENITES");
Loc18_AR_Source=T.NormCirc(T.Location==18&T.Lithology=="ARENITES"&T.Origin=="SOURCE");
Loc18_AR_River=T.NormCirc(T.Location==18&T.Lithology=="ARENITES"&T.Origin=="RIVER");

Loc18_AR_BetaParam=betafit(Loc18_AR);
Loc18_Source_BetaParam=betafit(Loc18_AR_Source);
Loc18_River_BetaParam=betafit(Loc18_AR_River);

% BinWidth=0.03;
[N_All,edges_All] = histcounts(Loc18_AR,'BinWidth',BinWidth);%'BinLimits',[0,1],'BinMethod','auto',
bins_All=conv(edges_All,[0.5,0.5],'valid');
[pks_All,locs_All] = findpeaks(N_All);
pks_All_Val=bins_All(locs_All);

[N_Source,edges_Source] = histcounts(Loc18_AR_Source,'BinWidth',BinWidth);%'BinLimits',[0,1],'BinMethod','auto',
bins_Source=conv(edges_Source,[0.5,0.5],'valid');
[pks_Source,locs_Source] = findpeaks(N_Source);
pks_Source_Val=bins_Source(locs_Source);

[N_River,edges_River] = histcounts(Loc18_AR_River,'BinWidth',BinWidth);%'BinLimits',[0,1],'BinMethod','auto',
bins_River=conv(edges_River,[0.5,0.5],'valid');
[pks_River,locs_River] = findpeaks(N_River);
pks_River_Val=bins_River(locs_River);

x_data=[0:0.01:1];

Fig_Test=figure
RGB = orderedcolors("gem");
hold on

histogram(Loc18_AR_Source,'BinEdges',edges_Source,'normalization','pdf','FaceColor',RGB(2,:),'DisplayName',"Source Data")
plot(x_data,betapdf(x_data,Loc18_Source_BetaParam(1),Loc18_Source_BetaParam(2)),'Color',RGB(2,:),'DisplayName',"Source Fit")

histogram(Loc18_AR_River,'BinEdges',edges_River,'normalization','pdf','FaceColor',RGB(1,:),'DisplayName',"River Data")
plot(x_data,betapdf(x_data,Loc18_River_BetaParam(1),Loc18_River_BetaParam(2)),'Color',RGB(1,:),'DisplayName',"River Fit")

% histogram(Loc18_AR,'BinEdges',edges_All,'normalization','pdf','FaceColor',RGB(3,:),'DisplayName',"All Data") %Loc18_AR,,'BinCounts',N_All
% plot(x_data,betapdf(x_data,Loc18_AR_BetaParam(1),Loc18_AR_BetaParam(2)),'Color',RGB(3,:),'DisplayName',"All Fit")

plot_diff=betapdf(x_data,Loc18_River_BetaParam(1),Loc18_River_BetaParam(2))-betapdf(x_data,Loc18_Source_BetaParam(1),Loc18_Source_BetaParam(2));
plot(x_data,plot_diff,'--','Color',RGB(4,:),'DisplayName',"pdf difference")

xlim([0.5 1])
xlabel('Normalized Isoperimetric Ratio [-]');
ylabel('Normalized Isoperimetric Ratio pdf [-]');
legend('Location','northwest')
legend('boxoff')
hold off;