function Circularity_Fit_Model_Quantiles_Common_Ka_fixed_C0(Lithology,LithologyColor,X,Parameters,Distance_Units,Distance_Metric,PathNum,Circularity_Fit_Model_Quantiles_file,mu_max,C_max,a,Ka)
% Fitting_Param,FitTest,CLowerC0,                               
%--------------------------------------------------------------------------
arguments
    Lithology (1,:) string
    LithologyColor (1,3) double
    X (:,:) cell
    Parameters (1,:) string             % e.g. ["Average" "q_5"]
    Distance_Units (1,:) string         % e.g. "km", or "m"
    Distance_Metric (1,:) string        % e.g. "SourceDist" or "TravDist"
    % Fitting_Param (1,1) string
    % FitTest (1,1) string
    % CLowerC0 (1,1) string
    PathNum double
    Circularity_Fit_Model_Quantiles_file (:,:) string
    mu_max double
    C_max double
    a (:,:) double 
    Ka (:,:) double
end
%--------------------------------------------------------------------------
tic
disp('Running Circularity Fit Quantiles COmmon Ka Model')
%--------------------------------------------------------------------------
matfile=Circularity_Fit_Model_Quantiles_file;
Eq_num=10;
%--------------------------------------------------------------------------
% Lithology with propoer font to use in figures and table name
Str=Lithology;
Lith_case=strcat(upper(Str{1}(1)),lower(Str{1}(2:end)));
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
%--------------------------------------------------------------------------
%                   OPTIMIZATION PARAMETERS MODEL BEGINS HERE
%--------------------------------------------------------------------------
% C_E_quant=cell(1,size(Variables_Names,2));
% mu_L_quant=cell(1,size(Variables_Names,2));
% C_L_quant=cell(1,size(Variables_Names,2));
RMSE_Ka=zeros(size(Ka));
RMSE_Ka_q=zeros(size(Ka,2),size(Variables_Names,2));

for n=1:size(Variables_Names,2)
    C_E_quant{n}=X{PathNum,1}.(Variables_Names(n));         % Isoperimetric Ratio measured (Empirical) in the field in each location along a certain path
    C0(n)=X{PathNum,1}.(Variables_Names(n))(1);
end

for k=1:size(Ka,2)
    % mu_L_quant{n}=1-exp(-Ka(k).*L);
    mu_L_quant{k}=mu_of_L_func(L,mu_max,Ka(k));
    
    for n=1:size(Variables_Names,2)
        RMSE_A=zeros(size(a,2),1);              % Root Mean Square Error for each combination of model parameters: C0, a for a fixed combination of Ka and for quantile n
        MAPE_A=zeros(size(a,2),1);
        C_L_a_quant=cell(size(a,2));
            for j=1:size(a,2)

                C_L_a_quant{j}(:,1)=C_of_mu_func(mu_L_quant{k},mu_max,C0(n),C_max,a(j));

                RMSE_A(j)=rmse(C_E_quant{n},C_L_a_quant{j});
                MAPE_A(j)=mape(C_E_quant{n},C_L_a_quant{j});

            end
        %------------------------------------------------------------------
        [RMSE_Ka_q(k,n),I_a(k,n)]=min(RMSE_A,[],'all','linear');
        [MAPE_Ka_q(k,n),I_MAPE_a(k,n)]=min(MAPE_A,[],'all','linear');
    end
    RMSE_Ka(k)=sum(RMSE_Ka_q(k,n),2);
end
[RMSE_Global,I_Ka]=min(RMSE_Ka,[],'all','linear');
Ka_opt=Ka(I_Ka);

for n=1:size(Variables_Names,2)
    % Convert linear indices I_C0_a into row,col indices basedon the
    % dimensions of the parameters space C0 and a, for the Ka that gives
    % the minimum RMSE over all quantiles considered
    % [I_a(n)]=ind2sub([size(a,2)],I_a(I_Ka,n));
    Opt_Param_C_L_quant_common_Ka{n,1}=C0(n);
    Opt_Param_C_L_quant_common_Ka{n,2}=a(I_a(I_Ka,n));%a(I_a(n));
    Opt_Param_C_L_quant_common_Ka{n,3}=Ka(I_Ka);
    Opt_Param_C_L_quant_common_Ka{n,4}=RMSE_Ka_q(I_Ka,n);%RMSE_Global;
    Opt_Param_C_L_quant_common_Ka{n,5}=MAPE_Ka_q(I_Ka,n);
    Opt_Param_C_L_quant_common_Ka{n,6}=a(I_a(I_Ka,n));
end
%--------------------------------------------------------------------------
Table_OptPar=cell2table(Opt_Param_C_L_quant_common_Ka);
Table_OptPar.Properties.VariableNames={'C0','a','Ka','RMSE','MAPE','a_MAPE'};
TableFormat='.xlsx';
TableName=[sprintf('%s_Optimum_Parameters_Quantiles_',Lith_case),datestr(now,'mmmm dd yyyy hhMM'),TableFormat];
writetable(Table_OptPar,TableName);
%--------------------------------------------------------------------------
%                   OPTIMIZATION PARAMETERS MODEL ENDS HERE
%--------------------------------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PLOTTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
%                   PLOT Fig_Quant1
addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
addpath 'C:\Program Files\gs\gs10.05.1'
color_palette="imola";
color_palette_path=sprintf("D:\\02_DATA\\CODES\\MATLAB Codes\\ScientificColourMaps8\\%s\\%s.mat",color_palette,color_palette);
load(color_palette_path)
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
Fig_Quant1=figure(FigNum+1);
hold on;
ax_Quant1=gca;
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(0,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=8.89;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=4;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(Fig_Quant1, 'PaperPosition', [positionVect]) %Set paper total dimension
%--------------------------------------------------------------------------
%                   FONT NAMEcolormap
% listfonts
fontname = 'Times New Roman';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(0,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation

set(findall(Fig_Quant1,'-property','Box'),'Box','off') % optional
set(findall(Fig_Quant1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Fig_Quant1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%--------------------------------------------------------------------------
%                   FONT SIZE
FontSize = 17;
set(findall(Fig_Quant1,'-property','FontSize'),'FontSize',FontSize) % adjust fontsize to your document
% set(0,'DefaultAxesFontSize',FontSize); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
% colororder(imola)
colororder_name=parula;
Brightest_remove=50;
% colororder=imola(round(linspace(1,size(imola,1),size(Parameters,2))),:);
colororder=colororder_name(round(linspace(1,size(colororder_name,1)-Brightest_remove,size(Parameters,2))),:);
colormap(colororder_name(1:size(colororder_name,1)-Brightest_remove,:));
%--------------------------------------------------------------------------
L_plot=min(L):ceil(max(L)-min(L))/1000:max(L);
% mu_L_quant_plot=1-exp(-Ka(I3_C_L_quant).*L_plot);
mu_L_quant_plot=mu_of_L_func(L_plot,mu_max,Ka(I_Ka));

for n=1:size(Variables_Names,2)
    scatter(L,C_E_quant{n},[],colororder(n,:),'filled','o','DisplayName',string(erase(Variables_Names(n),"q_"))+'^{th} percentile','HandleVisibility','off');%+'^{th} percentile'
    %--------------------------------------------------------------------------
    % PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA    

    % C_L_quant_plot=C0(I1_C_L_quant)+(C_max-C0(I1_C_L_quant)).*(1-(1-mu_L_quant_plot./mu_max).^a(I2_C_L_quant)).^(1/a(I2_C_L_quant));
    C_L_quant_plot{n}=C_of_mu_func(mu_L_quant_plot,mu_max,C0(n),C_max,a(I_a(I_Ka,n)));
    %--------------------------------------------------------------------------

    plot(L_plot,C_L_quant_plot{n},'Color',colororder(n,:),'HandleVisibility','off');%'MarkerSize',10
end
%--------------------------------------------------------------------------
% UNCOMMENT TO ADD COLORBAR
% for z=1:size(Parameters,2)
%     YTicks(z)=str2num(Parameters{z}(3:end))/100;
%     TicksLabels{z}=strcat(Parameters{z}(3:end),' %');%replace(Parameters{z},'_',' ')
% end
% c=colorbar('Ticks',YTicks,...
%     'TickLabels',TicksLabels);
%--------------------------------------------------------------------------
ax_Quant1.YAxis.TickLabelFormat = '%.2f';
ax_Quant1.XAxis.TickLabelFormat = '%.2f';


xlabel(sprintf('Distance from sediment source [%s]',Distance_Units),'fontsize',FontSize);
ylabel('$IR_{n}$ [-]','fontsize',FontSize);
legend('Location','southeast')
legend('boxoff')
%--------------------------------------------------------------------------
hold off;
%--------------------------------------------------------------------------
%                   SAVING Fig_Quant1 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Long_Profile_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Long_Profile_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(Fig_Quant1,sprintf('%s_Calibration_Long_Profile_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), '.fig'));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   PLOT Fig_Quant2
FigNum=length(findobj('type','figure'));
Fig_Quant2=figure(FigNum+1);
hold on;
ax_Quant2=gca;
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(Fig_Quant2,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=8.89;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=4;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(Fig_Quant2, 'PaperPosition', [positionVect]) %Set paper total dimension
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   FONT NAME
% listfonts
fontname = 'Times New Roman';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(0,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation

set(findall(Fig_Quant2,'-property','Box'),'Box','off') % optional
set(findall(Fig_Quant2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Fig_Quant2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%--------------------------------------------------------------------------
%                   FONT SIZE
FontSize = 17;
set(findall(Fig_Quant2,'-property','FontSize'),'FontSize',FontSize) % adjust fontsize to your document
% set(0,'DefaultAxesFontSize',FontSize); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
% colororder(imola)
colororder_name=parula;
Brightest_remove=50;
% colororder=imola(round(linspace(1,size(imola,1),size(Parameters,2))),:);
colororder=colororder_name(round(linspace(1,size(colororder_name,1)-Brightest_remove,size(Parameters,2))),:);
colormap(colororder_name(1:size(colororder_name,1)-Brightest_remove,:));
%--------------------------------------------------------------------------
for n=1:size(Variables_Names,2)
    %----------------------------------------------------------------------
    % PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
    C_E_quant_scatter{n}=X{PathNum,1}.(Variables_Names(n));

    % mu_L_quant_scatter{n}=1-exp(-Opt_Param_C_L_quant{n,3}.*L);

    mu_L_quant_scatter{n}=mu_of_L_func(L,mu_max,Opt_Param_C_L_quant_common_Ka{n,3});

    % C_L_quant_scatter{n}=Opt_Param_C_L_quant{n,1}+(C_max-Opt_Param_C_L_quant{n,1}).*(1-(1-mu_L_quant_scatter{n}./mu_max).^Opt_Param_C_L_quant{n,2}).^(1/Opt_Param_C_L_quant{n,2});
    C_L_quant_scatter{n}=C_of_mu_func(mu_L_quant_scatter{n},mu_max,Opt_Param_C_L_quant_common_Ka{n,1},C_max,Opt_Param_C_L_quant_common_Ka{n,2});
    
    scatter(C_E_quant_scatter{n},C_L_quant_scatter{n},[],colororder(n,:),'filled','o','DisplayName',string(erase(Variables_Names(n),"q_"))+'^{th} percentile','HandleVisibility','off');%+'^{th} percentile' string(Quantiles_Used(i))+'^{th} percentile' ,string(Variables_Names(n))

    %----------------------------------------------------------------------
end
plot([min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], [min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], '--k','HandleVisibility','off')

%--------------------------------------------------------------------------
% UNCOMMENT TO ADD COLORBAR
for z=1:size(Parameters,2)
    YTicks(z)=str2num(Parameters{z}(3:end))/100;
    TicksLabels{z}=strcat(Parameters{z}(3:end),'%');%replace(Parameters{z},'_',' ')
end
c=colorbar('Ticks',YTicks,...
    'TickLabels',TicksLabels);
%--------------------------------------------------------------------------
% axis equal
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
xtl = linspace(min(yt), max(yt), numel(yt));
set(gca,'XTick',xtl,'XTickLabel',compose('%.2f',xtl))

ax_Quant2.YAxis.TickLabelFormat = '%.2f';
ax_Quant2.XAxis.TickLabelFormat = '%.2f';
grid on

xlabel(sprintf('Measured $IR_{n}$ [-]'),'fontsize',FontSize);%sprintf('Measured %s [-]',"Parameter")
ylabel(sprintf('Modelled $IR_{n}$ [-]'),'fontsize',FontSize);%sprintf('Modelled %s [-]',"Parameter")
legend('Location','Southeast')
legend('boxoff')
hold off;
%--------------------------------------------------------------------------
%                   SAVING Fig_Quant1 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Scatterplot_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Scatterplot_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(Fig_Quant1,sprintf('%s_Calibration_Scatterplot_Quantiles_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), '.fig'));
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