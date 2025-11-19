function Circularity_Fit_Model_fixed_C0_Distances_Samples(Lithology,LithologyColor,X,Parameter,Distance_Units,Distance_Metric,DistanceWeigthing,RastersData, RastersLocations,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file,mu_max,C_max,a,Ka,lambda)
% 
%--------------------------------------------------------------------------
arguments
    Lithology (1,:) string
    LithologyColor (1,3) double
    X (:,:) cell
    Parameter (1,1) string              %
    Distance_Units (1,:) string         % e.g. "km", or "m"
    Distance_Metric (1,:) string        % e.g. "SourceDist" or "TravDist"
    DistanceWeigthing (1,1) string
    RastersData cell
    RastersLocations (:,:) double
    Fitting_Param (1,1) string
    FitTest (1,1) string
    CLowerC0 (1,1) string
    PathNum double
    Circularity_Fit_Model_file (:,:) string
    mu_max double
    C_max double

    a (:,:) double %=[0.2:0.1:18.0];
    Ka (:,:) double %=[1.0:0.1:12];
    lambda (:,:) double 
end
%--------------------------------------------------------------------------
tic
disp('Running Circularity Fit model')
%--------------------------------------------------------------------------
matfile=Circularity_Fit_Model_file;
Eq_num=12;
%--------------------------------------------------------------------------
% Lithology with propoer font to use in figures and table name
Str=Lithology;
Lith_case=strcat(upper(Str{1}(1)),lower(Str{1}(2:end)));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if Distance_Metric=="SourceDist"
    L=cell(size(X{PathNum,1}(:,1)));
    for i=1:length(L)
        L{i,1}=X{ PathNum,1}.(Distance_Metric)(i,1);
    end
    % if or(Distance_Units=="km",Distance_Units=="Km")
    %     L=X{ PathNum,1}.(Distance_Metric);                      % Length between source and measurement location (L is in kilometers [km])
    % elseif or(Distance_Units=="m",Distance_Units=="m")
    %     L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
    % else
    %     L=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
    % end    
elseif Distance_Metric=="TravDist"
    idx=0;
    for m=1:length(RastersLocations)
        if ismember(RastersLocations(m),X{PathNum,1}.Location)
            idx=idx+1;
            L{idx,1}=RastersData{m,1};
        end
    end
end
%--------------------------------------------------------------------------
%                   OPTIMIZATION PARAMETERS MODEL BEGINS HERE
%--------------------------------------------------------------------------
C_E=X{PathNum,1}.Average;                                   % Isoperimetric Ratio measured (Empirical) in the field in each location along a certain path
C0=X{PathNum,1}.Average(1);
clear mu_L C_L_AVG 
if DistanceWeigthing=="Y"
    lambda_N=length(lambda);
else
    lambda_N=1;
end
RMSE_C_L=zeros(lambda_N,size(a,2),size(Ka,2));
MAPE_C_L=zeros(lambda_N,size(a,2),size(Ka,2));
for i=1:lambda_N
    for j=1:size(a,2)
        for k=1:size(Ka,2)
            C_L_AVG=[];
            idx=0;
            for m=1:length(X{PathNum,1}.Location) %length(RastersLocations)
                % if ismember(RastersLocations(m),X{PathNum,1}.Location)
                    idx=idx+1;
                    mu_L{m,1}(1,:)=1-exp(-Ka(k).*L{m,1}(1,:));

                    C_L_Distrib_Sources{m,1}(1,:)=C0+(C_max-C0).*(1-(1-mu_L{m,1}(1,:)./mu_max).^a(j)).^(1/a(j));

                    if DistanceWeigthing=="Y"
                        
                        P_L{m,1}(1,:)=1./lambda(i).*exp(-L{m,1}(1,:)./lambda(i));
                        
                        C_L_AVG(idx,1)=sum(C_L_Distrib_Sources{m,1}(1,:).*P_L{m,1}(1,:))/sum(P_L{m,1}(1,:));

                    else
                        C_L_AVG(idx,1)=mean(C_L_Distrib_Sources{m,1}(1,:));
                    end
                % end

            end

            RMSE_C_L(i,j,k)=rmse(C_L_AVG,C_E);
            MAPE_C_L(i,j,k)=mape(C_L_AVG,C_E);
            clear mu_L C_L_Distrib_Sources P_L C_L_AVG
        end
    end
end
%--------------------------------------------------------------------------
[min_RMSE_C_L,I_C_L]=min(RMSE_C_L,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_C_L,I2_C_L,I3_C_L]=ind2sub([size(RMSE_C_L,1) size(RMSE_C_L,2) size(RMSE_C_L,3)],I_C_L);

Opt_Param_C_L{1}=C0;

if DistanceWeigthing=="Y"
    Opt_Param_C_L{2}=lambda(I1_C_L);
    lambda_Opt=lambda(I1_C_L);
else
    Opt_Param_C_L{2}=-999.9;
    lambda_Opt=-999.9;
end
Opt_Param_C_L{3}=a(I2_C_L);
a_Opt=a(I2_C_L);
Opt_Param_C_L{4}=Ka(I3_C_L);
Ka_Opt=Ka(I3_C_L);
Opt_Param_C_L{5}=RMSE_C_L(I1_C_L,I2_C_L,I3_C_L);
Opt_Param_C_L{6}=MAPE_C_L(I1_C_L,I2_C_L,I3_C_L);
%--------------------------------------------------------------------------
Table_OptPar=cell2table(Opt_Param_C_L);
Table_OptPar.Properties.VariableNames={'C0','lambda','a','Ka','RMSE','MAPE'};
TableFormat='.xlsx';
if DistanceWeigthing=="Y"
    Dist_wei="Lambda";
else
    Dist_wei="No_Lambda";
end
TableName=[sprintf('%s_Optimum_Parameters_AVG_%s_',Lith_case,Dist_wei),datestr(now,'mmmm dd yyyy hhMM'),TableFormat];
writetable(Table_OptPar,TableName);
%--------------------------------------------------------------------------
%                   OPTIMIZATION PARAMETERS MODEL ENDS HERE
%--------------------------------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PLOTTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
%                   PLOT Fig_1
addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
addpath 'C:\Program Files\gs\gs10.05.1'
color_palette="imola";
color_palette_path=sprintf("D:\\02_DATA\\CODES\\MATLAB Codes\\ScientificColourMaps8\\%s\\%s.mat",color_palette,color_palette);
load(color_palette_path)
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
fig_1=figure(FigNum+1);
hold on;
ax_1=gca;
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(0,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=8.89;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=4;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(fig_1, 'PaperPosition', [positionVect]) %Set paper total dimension
%--------------------------------------------------------------------------
%                   FONT NAME
% listfonts
fontname = 'Times New Roman';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(0,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation

set(findall(fig_1,'-property','Box'),'Box','off') % optional
set(findall(fig_1,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig_1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%--------------------------------------------------------------------------
%                   FONT SIZE
FontSize = 17;
set(findall(fig_1,'-property','FontSize'),'FontSize',FontSize) % adjust fontsize to your document
% set(0,'DefaultAxesFontSize',FontSize); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
colororder(imola)
ColOrder=imola(size(Parameter,2));
%--------------------------------------------------------------------------
% PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
if or(Distance_Units=="km",Distance_Units=="Km")
    X_plot=X{ PathNum,1}.(Distance_Metric);                      % Length between source and measurement location (L is in kilometers [km])
elseif or(Distance_Units=="m",Distance_Units=="m")
    X_plot=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
else
    X_plot=X{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
end

field_data=errorbar(X_plot,C_E,X{PathNum,1}.SEM,"o","MarkerSize",8,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Field data');%"s", "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Field data';    "o","MarkerSize",8,"MarkerFaceColor",imola(1,:),'DisplayName','Field data','HandleVisibility','on'
%--------------------------------------------------------------------------
C_L_AVG=[];
idx=0;
for m=1:length(X{PathNum,1}.Location)
    % if ismember(RastersLocations(m),X{PathNum,1}.Location)
        idx=idx+1;
        L_plot(idx,1)=mean(L{m,1}(1,:));
        mu_L{idx,1}(1,:)=1-exp(-Ka_Opt.*L{m,1}(1,:));

        C_L_Distrib_Sources{idx,1}(1,:)=C0+(C_max-C0).*(1-(1-mu_L{idx,1}(1,:)./mu_max).^a_Opt).^(1/a_Opt);
        
        if DistanceWeigthing=="Y"

            P_L{idx,1}(1,:)=1./lambda_Opt.*exp(-L{m,1}(1,:)./lambda_Opt);

            C_L_AVG(idx,1)=sum(C_L_Distrib_Sources{idx,1}(1,:).*P_L{idx,1}(1,:))/sum(P_L{idx,1}(1,:));

        else
            C_L_AVG(idx,1)=mean(C_L_Distrib_Sources{idx,1}(1,:));%[C_L_AVG; mean(C_L_Distrib_Sources{idx,1}(1,:))]
        end

    % end
end

if Distance_Metric=="SourceDist"
    X_plot=[0:0.001:X{ PathNum,1}.(Distance_Metric)(end)];
    mu_plot=1-exp(-Ka_Opt.*X_plot);
    C_L_plot=C0+(C_max-C0).*(1-(1-mu_plot./mu_max).^a_Opt).^(1/a_Opt);
    Y_plot=C_L_plot;
elseif Distance_Metric=="TravDist"
    X_plot=X_plot;
    Y_plot=C_L_AVG;
end

% plot(L_plot,C_L_AVG)
plot(X_plot,Y_plot,'LineWidth',linewidth,'Color',LithologyColor,'DisplayName',sprintf('Equation %0.f',Eq_num),'HandleVisibility','on')
%--------------------------------------------------------------------------
ax_1.YAxis.TickLabelFormat = '%.2f';
ax_1.XAxis.TickLabelFormat = '%.2f';

xlabel(sprintf('Distance from sediment source [%s]',Distance_Units),'fontsize',FontSize);
ylabel('$IR_{n}$ [-]','fontsize',FontSize);
legend('Location','southeast')
legend('boxoff')
%--------------------------------------------------------------------------
% yyaxis right
% plot(L_plot,mu_L_plot)
% ylabel('Relative Mass Loss [-]','fontsize',15);
hold off;
%--------------------------------------------------------------------------
%                   SAVING fig_1 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Long_Profile_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Long_Profile_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,sprintf('%s_Calibration_Long_Profile_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), '.fig'));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   PLOT Fig_2
FigNum=length(findobj('type','figure'));
fig_2=figure(FigNum+1);
hold on;
ax_2=gca;
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(0,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=8.89;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=4;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(fig_2, 'PaperPosition', [positionVect]) %Set paper total dimension
%--------------------------------------------------------------------------
%                   FONT NAME
% listfonts
fontname = 'Times New Roman';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(0,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation

set(findall(fig_2,'-property','Box'),'Box','off') % optional
set(findall(fig_2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%--------------------------------------------------------------------------
%                   FONT SIZE
FontSize = 17;
set(findall(fig_2,'-property','FontSize'),'FontSize',FontSize) % adjust fontsize to your document
% set(0,'DefaultAxesFontSize',FontSize); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
colororder(imola)
ColOrder=imola(size(Parameter,2));
%--------------------------------------------------------------------------
% PLOT OF OPTIMIZATION MODEL RESULTS TO VISUALLY COMPARE WITH FIELD DATA
C_E_scatter=C_E;

% mu_L_scatter=mu_of_L_func(L,mu_max,Ka(I2_C_L));
% 
% C_L_scatter=C_of_mu_func(mu_L_scatter,mu_max,C0,C_max,a(I1_C_L));

"o","MarkerSize",8,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Model vs. Measured','HandleVisibility','off'
Calib_Scatter=scatter(C_E_scatter,C_L_AVG,50,"o","MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Model vs. Measured','HandleVisibility','off'); %'filled','o'
%--------------------------------------------------------------------------

Calib_RefPlot=plot([min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], [min([xlim; ylim],[],'all') max([xlim; ylim],[],'all')], '--k','HandleVisibility','off')

% axis equal
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
xtl = linspace(min(yt), max(yt), numel(yt));
set(gca,'XTick',xtl,'XTickLabel',compose('%.2f',xtl))

ax_2.YAxis.TickLabelFormat = '%.2f';
ax_2.XAxis.TickLabelFormat = '%.2f';

axis tight
grid on

xlabel(sprintf('Measured $IR_{n}$ [-]'),'fontsize',FontSize);
ylabel(sprintf('Modelled $IR_{n}$ [-]'),'fontsize',FontSize);
% legend('Location','Southeast')
% legend('boxoff')
hold off;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   SAVING fig_1 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Scatterplot_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_Scatterplot_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,sprintf('%s_Calibration_Scatterplot_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), '.fig'));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   PLOT Fig_2
FigNum=length(findobj('type','figure'));
fig_3=figure(FigNum+1);
hold on;
ax_3=gca;
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(fig_3,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=16.00;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=4;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(fig_3, 'PaperPosition', [positionVect]) %Set paper total dimension

fig_3.Position=[12.8852, 6.3765, 16, 8];
%--------------------------------------------------------------------------
%                   FONT NAME
% listfonts
fontname = 'Times New Roman';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(0,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation

set(findall(fig_3,'-property','Box'),'Box','off') % optional
set(findall(fig_3,'-property','Interpreter'),'Interpreter','latex') 
set(findall(fig_3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')

%--------------------------------------------------------------------------
%                   FONT SIZE
FontSize = 17;
set(findall(fig_3,'-property','FontSize'),'FontSize',FontSize) % adjust fontsize to your document
% set(0,'DefaultAxesFontSize',FontSize); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
colororder(imola)
ColOrder=imola(size(Parameter,2));
%--------------------------------------------------------------------------
RMSE_plot(:,:)=RMSE_C_L(1,:,:);

RMSE_Surf=surf(Ka',a',RMSE_plot);
% RMSE_Surf=surf(a,ka,RMSE_plot);

RMSE_Min_Point=scatter3(Ka_Opt,a_Opt,min_RMSE_C_L,100,'r','filled');
% RMSE_Min_Point=scatter3(I3_C_L,I2_C_L,min_RMSE_C_L,100,'r','filled');

view(-140,40);
colorbar;

set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1.0]);
set(gca,'YTick',[unique(round(a))]);

ax_3.XAxis.TickLabelFormat = '%.1f';
ax_3.YAxis.TickLabelFormat = '%.1f';

xlabel(sprintf('$k_{a}$ [-]'),'fontsize',FontSize);
ylabel(sprintf('$a$ [-]'),'fontsize',FontSize);
zlabel(sprintf('$RMSE$ [-]'),'fontsize',FontSize);

hold off
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   SAVING fig_3 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_RMSE_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Calibration_RMSE_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), fileformat); %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,sprintf('%s_Calibration_RMSE_%s%s',Lithology,datestr(now,'mmmm dd yyyy hhMM'), '.fig'));

%--------------------------------------------------------------------------
save(matfile);
toc