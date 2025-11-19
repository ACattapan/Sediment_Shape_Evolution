FigNum=length(findobj('type','figure'));
Fig_1=figure(FigNum+1);
Tileplot=tiledlayout(size(Locations_Analysis,2),1);

%--------------------------------------------------------------------------
% DEFINE COLOUR MAP TO USE COLOURS THAT SCALE WITH DISTANCE FROM BASIN OUTLET
map=zeros(size(Locations_Analysis,2),3);
for m=1:size(Locations_Analysis,2)
    map(m,1)=LocDist(Locations_Analysis(m))/LocDist(min(Locations_Analysis));
    map(m,2)=0;
    map(m,3)=0;
end
cmap=colormap(map);
%--------------------------------------------------------------------------
StatFlag="Local_a_Ka"; %"Average", "Distribution", "Local_a_Ka"
Opt_Val_Flag="Local_a_Ka";%"Avg", "M", "Dist", "Local_a_Ka"
Fitting_Param="Y";  %"Y", "N"

FitTest="AD";       %"KS", "AD"
alpha=0.05;
Plotting_flag="Y";  %"Y", "N"
PlotType="pdf";     %"pdf", "cdf"

for m=1:size(Locations_Analysis,2)
    % subplot(size(Locations_Analysis,2),1,m)
    nexttile
    ax=gca;
    ax.FontSize = 11;
    RasterData_TDtest=[];

    RasterData_TDtest=TD_Rasters_Analysis{1,m}(:);
    RasterSize=size(RasterData_TDtest,1);
    histogram(RasterData_TDtest,'Normalization','pdf','FaceColor',[0 0.4470 0.7410],'DisplayName','Lidar raster');
    hold on;
    Sed_Sample_Size=size(T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)),1);
    if Sed_Sample_Size>0
        C_plot=[];
        mu_plot=[];
        Ls_plot=[];
        TD_plot_Size=max(Sed_Sample_Size,RasterSize);
        if Opt_Val_Flag=="Avg"
            Opt_Val_Set=Opt_Param_Avg;
        elseif Opt_Val_Flag=="M"
            Opt_Val_Set=Opt_Param_M;
        elseif Opt_Val_Flag=="Local_a_Ka"
            Opt_Val_Set=Opt_Param_Dist_Loc_aKa;
        elseif Opt_Val_Flag=="Dist"
            Opt_Val_Set=Opt_Param_Dist;
        end
        if StatFlag=="Average"
            if or(Fitting_Param=="Y",Fitting_Param=="y")
                [beta1,beta2,h,p]=BetaFit(T,Parameter,FitTest,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                C_plot{m}=betarnd(beta1,beta2,TD_plot_Size,1);
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2});
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
            else
                C_plot{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m));
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2});
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
            end
        elseif StatFlag=="Local_a_Ka"
            if or(Fitting_Param=="Y",Fitting_Param=="y")
                [beta1,beta2,h,p]=BetaFit(T,Parameter,FitTest,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                C_plot{m}=betarnd(beta1,beta2,TD_plot_Size,1);
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2}(1,m));
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
                % ax1=gca;
                % ax1.FontSize = 11;
            else
                C_plot{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m));
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2}(1,m));
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
            end
        elseif StatFlag=="Distribution"
            if or(Fitting_Param=="Y",Fitting_Param=="y")
                [beta1,beta2,h,p]=BetaFit(T,Parameter,FitTest,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                C_plot{m}=betarnd(beta1,beta2,TD_plot_Size,1);
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2});
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
                % ax1=gca;
                % ax1.FontSize = 11;
            else
                C_plot{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m));
                mu_plot{m}=mu_func(C_plot{m}(:),mu_max,Opt_Val_Set{1},C_max,Opt_Val_Set{2});
                Ls_plot{m}=Length_func(mu_plot{m}(:),10^-Opt_Val_Set{3}(1,m));
                histogram(Ls_plot{m}(:),'Normalization','pdf','FaceColor',[map(m,1) map(m,2) map(m,3)],'DisplayName','Model estimate');
            end
        end

    end
    if m<=2
        xlim([0 100]);
    elseif m<=5
        xlim([0 800]);
    elseif m<=7
        xlim([0 1600]);
    else
        xlim([1000 6000]);
    end
    if m==1 %size(LocPath{PathNum},2)
        Lgnd = legend('FontSize',14);
        %         Lgnd.Position=[0.755 0.825 0.2 0.1];
        Lgnd.Location='northeast';
        %         Lgnd.Position(1) = 0;
        %         Lgnd.Position(2) = 0;
    end
end
% labelx=xlabel('Sediments'' Travel distance [km]');
% labely=ylabel('Travel distance pdf [km^{-1}]');
% labely.Position(1)=-0;
% labely.Position(2)=10;
% set(Fig_1,'Colormap',cmap);
% caxis([0 7000])
% colorbar;
%Set the current Colormap
% cb=colorbar('Ticks',linspace(0,1,11),'TickLabels',string(linspace(0,1,11).*7000)); %linspace(0,1,10).*7000
% cb.Layout.Tile='east';
% colormap(cmap)

% han=axes(Fig_1,'visible','off');
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Travel distance pdf [km^{-1}]','FontSize',14);
% xlabel(han,'Sediments'' Travel distance [km]','FontSize',14);
% title(han,'yourTitle');

ylabel(Tileplot,'Travel distance pdf [km^{-1}]','FontSize',14);
xlabel(Tileplot,'Sediments'' Travel distance [km]','FontSize',14);