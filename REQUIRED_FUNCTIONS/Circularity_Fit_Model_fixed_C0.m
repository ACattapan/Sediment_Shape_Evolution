function Circularity_Fit_Model_fixed_C0(X,Parameter,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file,mu_max,C_max,a,Ka)
% 
%--------------------------------------------------------------------------
arguments
    X (:,:) cell
    Parameter (1,1) string              %
    Distance_Units (1,:) string         % e.g. "km", or "m"
    Distance_Metric (1,:) string        % e.g. "SourceDist" or "TravDist"
    Fitting_Param (1,1) string
    FitTest (1,1) string
    CLowerC0 (1,1) string
    PathNum double
    Circularity_Fit_Model_file (:,:) string
    mu_max double
    C_max double

    a (:,:) double %=[0.2:0.1:18.0];
    Ka (:,:) double %=[1.0:0.1:12];
end
%--------------------------------------------------------------------------
tic
disp('Running Circularity Fit model')
%--------------------------------------------------------------------------
matfile=Circularity_Fit_Model_file;
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
F1=figure(FigNum+1);
hold on;
FSize=12;
ColOrder=parula(size(Parameter,2));
% errorbar(X{PathNum,1}.TravDist,X{PathNum,1}.Average,X{PathNum,1}.SEM,"diamond","MarkerFaceColor",[0.85 0.33 0.10],"MarkerSize",5,"LineWidth", 1)
% L_plot=linspace(min(X{PathNum,1}.TravDist)*1000,max(X{PathNum,1}.TravDist)*1000,100);
%--------------------------------------------------------------------------
if or(Distance_Units=="km",Distance_Units=="Km")
    L=X{ PathNum,1}.(Distance_Metric);                      % Length between source and measurement location (L is in kilometers [km])
elseif or(Distance_Units=="m",Distance_Units=="m")
    L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
else
    L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
end
%--------------------------------------------------------------------------
%                   OPTIMIZATION PARAMETERS MODEL BEGINS HERE
%--------------------------------------------------------------------------
C_E=X{PathNum,1}.Average;                                   % Isoperimetric Ratio measured (Empirical) in the field in each location along a certain path
RMSE_C_L=zeros(size(a,2),size(Ka,2));
C0=X{PathNum,1}.Average(1);
clear mu_L C_L
    for j=1:size(a,2)
        for k=1:size(Ka,2)
            % mu_L=1-exp(-Ka(k).*L);
            mu_L=mu_of_L_func(L,mu_max,Ka(k));
            % C_L=C0(i)+(C_max-C0(i)).*(1-(1-mu_L./mu_max).^a(j)).^(1/a(j));
            C_L=C_of_mu_func(mu_L,mu_max,C0,C_max,a(j));

            RMSE_C_L(j,k)=rmse(C_L,C_E);
            clear mu_L C_L
        end
    end

%--------------------------------------------------------------------------
[min_RMSE_C_L,I_C_L]=min(RMSE_C_L,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_C_L,I2_C_L]=ind2sub([size(RMSE_C_L,1) size(RMSE_C_L,2)],I_C_L);

Opt_Param_C_L{1}=a(I1_C_L);
Opt_Param_C_L{2}=Ka(I2_C_L);
Opt_Param_C_L{3}=RMSE_C_L(I1_C_L,I2_C_L);
%--------------------------------------------------------------------------
Table_OptPar=cell2table(Opt_Param_C_L);
Table_OptPar.Properties.VariableNames={'a','Ka','RMSE'};
writetable(Table_OptPar,'Table_OptPar.xlsx');
%--------------------------------------------------------------------------
% PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA 
L_plot=min(L):ceil(max(L)-min(L))/100:max(L);
% mu_L_plot=1-exp(-Ka(I3_C_L).*L);
mu_L_plot=mu_of_L_func(L_plot,mu_max,Ka(I2_C_L));

% C_L=C0(I1_C_L)+(C_max-C0(I1_C_L)).*(1-(1-mu_L_plot./mu_max).^a(I2_C_L)).^(1/a(I2_C_L));
C_L_plot=C_of_mu_func(mu_L_plot,mu_max,C0,C_max,a(I1_C_L));

errorbar(L,C_E,X{PathNum,1}.SEM,"diamond","MarkerSize",5,"LineWidth", 1,'DisplayName','Field data')%"MarkerFaceColor",[0.85 0.33 0.10],

plot(L_plot,C_L_plot,'HandleVisibility','off');%'Color',ColOrder(n,:),...,'MarkerSize',10
%--------------------------------------------------------------------------
xlabel(sprintf('Distance from sediment source [%s]',Distance_Units),'fontsize',FSize);
ylabel('Isoperimetric Ratio [-]','fontsize',FSize);
legend('Location','southeast')
legend('boxoff')
%--------------------------------------------------------------------------
% yyaxis right
% plot(L_plot,mu_L_plot)
% ylabel('Relative Mass Loss [-]','fontsize',15);
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
F2=figure(FigNum+1);
hold on;
ColOrder=parula(size(Parameter,2));
%--------------------------------------------------------------------------
% PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
 C_E_scatter=C_E;

mu_L_scatter=mu_of_L_func(L,mu_max,Ka(I2_C_L));

C_L_scatter=C_of_mu_func(mu_L_scatter,mu_max,C0,C_max,a(I1_C_L));

scatter(C_E_scatter,C_L_scatter,[],'filled','o');
%--------------------------------------------------------------------------

plot([min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], [min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], '--k','HandleVisibility','off')
xlabel(sprintf('Measured %s [-]',Parameter),'fontsize',FSize);
ylabel(sprintf('Modelled %s [-]',Parameter),'fontsize',FSize);
% legend('Location','Southeast')
% legend('boxoff')
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% FigName='Modelled_%s_of_%s_along_path_number_%d.fig';
% saveas(F2,[Output_folder_path+'\'+sprintf(FigName,Parameter,Lithology,PathNum)])

%--------------------------------------------------------------------------
% if size(C0,2)>1
%     FigNum=length(findobj('type','figure'));
%     F2=figure(FigNum+1)
% 
%     RMSE_plot=RMSE_C_L(:,:,I3_C_L);
%     surf(RMSE_plot)
%     hold on
%     scatter3(I2_C_L,I1_C_L,RMSE_plot(I1_C_L,I2_C_L),100,'r','filled')
% end
%--------------------------------------------------------------------------
save(matfile);
toc