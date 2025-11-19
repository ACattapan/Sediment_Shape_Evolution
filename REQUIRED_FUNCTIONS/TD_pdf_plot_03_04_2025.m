% function TD_pdf_plot_03_04_2025(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
% h =  findobj('type','figure');
% n = length(h);

FigNum=length(findobj('type','figure'));
Fig1=figure(FigNum+1)

FigNum=length(findobj('type','figure'));
Fig2=figure(FigNum+1)
ax2=axes;
Xlim2=[0.5 1];

FigNum=length(findobj('type','figure'));
Fig3=figure(FigNum+1)
ax3=axes;
Xlim3=[0.85 1];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for m=1:size(Locations_Analysis,2)%size(LocPath{PathNum},2)
    RasterData_TDtest=[];
    set(0,'CurrentFigure',Fig1)
    subplot(size(Locations_Analysis,2),1,m)% subplot(size(LocPath{PathNum},2),1,m)
    ax1=gca;
    RasterData_TDtest=RastersData{1,RastersLocations==Locations_Analysis(m)}(:);%RastersData{1,RastersLocations==LocPath{PathNum}(m)}(:);
    RasterData_TDtest(isnan(RasterData_TDtest))=[];
    SampSz=size(RasterData_TDtest,1);
    StatFlag="Distribution";
    if size(T.Pebble_Id(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==OrigPath{PathNum}(m)),1)>0%size(T.Pebble_Id(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m)),1)>0
        C_plot=[];
        mu_plot=[];
        Ls_plot=[];
        histogram(ax1,RasterData_TDtest,'normalization','pdf');
        hold on;
        if StatFlag=="Average"
            if or(Fitting_Param=="Y",Fitting_Param=="y")
                [beta1,beta2,h,p]=BetaFit(T,Parameter,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                C=betarnd(beta1,beta2,SampSz,1);
                mu_plot=mu_func(C,mu_max,Opt_Param_avg(1,m),C_max,Opt_Param_avg(2,m));
                Ls_plot=Length_func(mu_plot,Opt_Param_avg(3,m));
                histogram(ax1,Ls_plot,'normalization','pdf');
            else
                C_plot=T.(Parameter)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m));
                mu_plot=mu_func(C_plot,mu_max,Opt_Param_avg(1,m),C_max,Opt_Param_avg(2,m));
                Ls_plot=Length_func(mu_plot(:),Opt_Param_avg(3,m));
                histogram(ax1,Ls_plot,'normalization','pdf');
            end
        elseif StatFlag=="Distribution"
            if or(Fitting_Param=="Y",Fitting_Param=="y")
                [beta1,beta2,h,p]=BetaFit(T,Parameter,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                C=betarnd(beta1,beta2,SampSz,1);
                mu_plot=mu_func(C,mu_max,Output_file_Travel_Dist(1,m),C_max,Output_file_Travel_Dist(2,m));
                Ls_plot=Length_func(mu_plot,Output_file_Travel_Dist(3,m));
                histogram(ax1,Ls_plot,'normalization','pdf');
            else
                C_plot=T.(Parameter)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m));
                %----------------------------------------------------------
                set(0,'CurrentFigure',Fig2)
                ax2=subplot(size(Locations_Analysis,2),1,m)% subplot(size(LocPath{PathNum},2),1,m)
                histogram(ax2,C_plot,'normalization','pdf');
                xlim(ax2,Xlim2);
                hold on;
                %----------------------------------------------------------
                % mu_plot=mu_func(C_plot,mu_max,Output_file_Travel_Dist(1,m),C_max,Output_file_Travel_Dist(2,m));
                mu_plot=mu_func(C_plot,mu_max,Opt_Param_Dist{1},C_max,Opt_Param_Dist{1});
                %----------------------------------------------------------
                set(0,'CurrentFigure',Fig3)
                ax3=subplot(size(Locations_Analysis,2),1,m)% subplot(size(LocPath{PathNum},2),1,m)
                histogram(ax3,mu_plot,'normalization','pdf');
                xlim(ax3,Xlim3);
                hold on;
                %----------------------------------------------------------
                Ls_plot=Length_func(mu_plot(:),Opt_Param_Dist{3}(m));
                set(0,'CurrentFigure',Fig1)
                histogram(ax1,Ls_plot,'normalization','pdf');
            end
        end
    end
    % xlim([0 max([Ls_plot])]); %max([RasterData_TDtest;Ls_plot])])
end
