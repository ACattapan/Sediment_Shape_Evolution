function Circularity_Fit_Model_Quantiles_Ka(X,Parameters,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_Quantiles_file,mu_max,C_max,C0,a,Ka)
%
%--------------------------------------------------------------------------
arguments
    X (:,:) cell
    Parameters (1,:) string             % e.g. ["Average" "q_5"]
    Distance_Units (1,:) string         % e.g. "km", or "m"
    Distance_Metric (1,:) string        % e.g. "SourceDist" or "TravDist"
    Fitting_Param (1,1) string
    FitTest (1,1) string
    CLowerC0 (1,1) string
    PathNum double
    Circularity_Fit_Model_Quantiles_file (:,:) string
    mu_max double
    C_max double
    C0 (:,:) double
    a (:,:) double 
    Ka (:,:) double
end
%--------------------------------------------------------------------------
tic
disp('Running Circularity Fit Quantiles Model')
%--------------------------------------------------------------------------
matfile=Circularity_Fit_Model_Quantiles_file;
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
Fig_Quant1=figure(FigNum+1);
hold on;
FSize=12;
ColOrder=parula(size(Parameters,2));
% errorbar(X{ PathNum,1}.(Distance_Metric),X{PathNum,1}.Average,X{PathNum,1}.SEM,"diamond","MarkerFaceColor",[0.85 0.33 0.10],"MarkerSize",5,"LineWidth", 1)
% X{2,1}{:,Quantiles_Names(1)}
% L_plot=linspace(min(X{ PathNum,1}.(Distance_Metric)*1000,max(X{ PathNum,1}.(Distance_Metric)*1000,100);
%--------------------------------------------------------------------------
Variables_Names=Parameters;
Variables_col_indx=zeros(1,size(Variables_Names,2));
for i=1:size(Variables_Names,2)
    Variables_col_indx(i)=find(string(X{PathNum,1}.Properties.VariableNames) == Variables_Names(1,i));
end
%--------------------------------------------------------------------------
if or(Distance_Units=="km",Distance_Units=="Km")
    L=X{ PathNum,1}.(Distance_Metric);                      % Length between source and measurement location (L is in kilometers [km])
elseif or(Distance_Units=="m",Distance_Units=="m")
    L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
else
    L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
end

for n=1:size(Variables_Names,2)
    C_E_quant{n}=X{PathNum,1}.(Variables_Names(n));         % Isoperimetric Ratio measured (Empirical) in the field in each location along a certain path
    RMSE_C_L_quant=zeros(size(C0,2),size(a,2),size(Ka,2));  % Root Mean Square Error as a function of the model parameters: C0, a, Ka
    for i=1:size(C0,2)
        for j=1:size(a,2)
            for k=1:size(Ka,2)
                % mu_L_quant{n}=1-exp(-Ka(k).*L);
                mu_L_quant{n}=mu_of_L_func(L,mu_max,Ka(k));

                % C_L_quant{n}=C0(i)+(C_max-C0(i)).*(1-(1-mu_L_quant{n}(:)./mu_max).^a(j)).^(1/a(j));
                C_L_quant{n}=C_of_mu_func(mu_L_quant{n},mu_max,C0(i),C_max,a(j));

                RMSE_C_L_quant(i,j,k)=rmse(C_E_quant{n},C_L_quant{n});
            end
        end
    end
    %----------------------------------------------------------------------
    [min_RMSE_C_L_quant,I_C_L_quant]=min(RMSE_C_L_quant,[],'all','linear');
    % Convert linear indices I(m) into row,col indices basedon the
    % dimensions of the matrix R2
    [I1_C_L_quant,I2_C_L_quant,I3_C_L_quant]=ind2sub([size(RMSE_C_L_quant,1) size(RMSE_C_L_quant,2) size(RMSE_C_L_quant,3)],I_C_L_quant);
    Opt_Param_C_L_quant{n,1}=C0(I1_C_L_quant);
    Opt_Param_C_L_quant{n,2}=a(I2_C_L_quant);
    Opt_Param_C_L_quant{n,3}=Ka(I3_C_L_quant);
    Opt_Param_C_L_quant{n,4}=RMSE_C_L_quant(I1_C_L_quant,I2_C_L_quant,I3_C_L_quant);

    %----------------------------------------------------------------------
    % PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
    L_plot=min(L):ceil(max(L)-min(L))/100:max(L);
    % mu_L_quant_plot=1-exp(-Ka(I3_C_L_quant).*L_plot);
    mu_L_quant_plot=mu_of_L_func(L_plot,mu_max,Ka(I3_C_L_quant));

    % C_L_quant_plot=C0(I1_C_L_quant)+(C_max-C0(I1_C_L_quant)).*(1-(1-mu_L_quant_plot./mu_max).^a(I2_C_L_quant)).^(1/a(I2_C_L_quant));
    C_L_quant_plot=C_of_mu_func(mu_L_quant_plot,mu_max,C0(I1_C_L_quant),C_max,a(I2_C_L_quant));

    scatter(L,C_E_quant{n},[],ColOrder(n,:),'filled','o','DisplayName',string(erase(Variables_Names(n),"q_"))+'^{th} percentile');%+'^{th} percentile'

    plot(L_plot,C_L_quant_plot,'Color',ColOrder(n,:),'HandleVisibility','off');%'MarkerSize',10
    %----------------------------------------------------------------------
end
Table_OptPar=cell2table(Opt_Param_C_L_quant);
Table_OptPar.Properties.VariableNames={'C0','a','Ka','RMSE'};
writetable(Table_OptPar,'Table_OptPar.xlsx');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

xlabel(sprintf('Distance from sediment source [%s]',Distance_Units),'fontsize',FSize);
ylabel('Isoperimetric Ratio [-]','fontsize',FSize);
legend('Location','southeast')
legend('boxoff')
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
Fig_Quant2=figure(FigNum+1);
hold on;
ColOrder=parula(size(Parameters,2));
for n=1:size(Variables_Names,2)
    %----------------------------------------------------------------------
    % PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
    C_E_quant_scatter{n}=X{PathNum,1}.(Variables_Names(n));

    % mu_L_quant_scatter{n}=1-exp(-Opt_Param_C_L_quant{n,3}.*L);

    mu_L_quant_scatter{n}=mu_of_L_func(L,mu_max,Opt_Param_C_L_quant{n,3});

    % C_L_quant_scatter{n}=Opt_Param_C_L_quant{n,1}+(C_max-Opt_Param_C_L_quant{n,1}).*(1-(1-mu_L_quant_scatter{n}./mu_max).^Opt_Param_C_L_quant{n,2}).^(1/Opt_Param_C_L_quant{n,2});
    C_L_quant_scatter{n}=C_of_mu_func(mu_L_quant_scatter{n},mu_max,Opt_Param_C_L_quant{n,1},C_max,Opt_Param_C_L_quant{n,2});
    
    scatter(C_E_quant_scatter{n},C_L_quant_scatter{n},[],ColOrder(n,:),'filled','o','DisplayName',string(erase(Variables_Names(n),"q_"))+'^{th} percentile');%+'^{th} percentile' string(Quantiles_Used(i))+'^{th} percentile' ,string(Variables_Names(n))

    %----------------------------------------------------------------------
end
plot([min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], [min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], '--k','HandleVisibility','off')
xlabel(sprintf('Measured %s [-]',"Parameter"),'fontsize',FSize);
ylabel(sprintf('Modelled %s [-]',"Parameter"),'fontsize',FSize);
legend('Location','Southeast')
legend('boxoff')
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% FigName='Modelled_quantiles_of_%s_of_%s_along_path_number_%d.fig';
% saveas(Fig_Quant1,[Output_folder_path+'\'+sprintf(FigName,Parameter,Lithology,PathNum)])


% yyaxis right
% plot(X{ PathNum,1}.(Distance_Metric),mu_L_quant_plot)
% ylabel('Relative Mass Loss [-]','fontsize',15);
%--------------------------------------------------------------------------
% if size(C0,2)>1    
%     FigNum=length(findobj('type','figure'));
%     Fig_Quant2=figure(FigNum+1)
% 
%     RMSE_quant_plot=RMSE_C_L_quant(:,:,I3_C_L_quant);
%     surf(RMSE_quant_plot)
%     hold on
%     scatter3(I2_C_L_quant,I1_C_L_quant,RMSE_quant_plot(I1_C_L_quant,I2_C_L_quant),100,'r','filled')
% end
%--------------------------------------------------------------------------
save(matfile);
toc