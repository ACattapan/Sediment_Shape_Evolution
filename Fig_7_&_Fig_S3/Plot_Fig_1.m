% load('CIRC_FIT_MODEL_AVG_NormCirc_METABASALTS_PATH_2_SourceDist_PAPER.mat','C0','C_E','a_Opt','Ka_Opt','lambda_Opt','mu_L','L','L_plot')
load('CIRC_FIT_MODEL_AVG_NormCirc_METABASALTS_PATH_2.mat')
load('STATS_OVERVIEW_METABASALTS_NormCirc.mat')
load('DATA_IMPORT_METABASALTS_NormCirc.mat')
load('TRAVEL_DISTANCES_METABASALTS_NormCirc.mat')
%                   PLOT Fig_1
addpath(genpath('..'))
color_palette="imola";
color_palette_path=sprintf("%s\\ScientificColourMaps8\\%s\\%s.mat",fileparts(fullfile(pwd)),color_palette,color_palette);
load(color_palette_path)
%--------------------------------------------------------------------------
Eq_num=12;
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
    X_plot=X_calibration{ PathNum,1}.(Distance_Metric);                      % Length between source and measurement location (L is in kilometers [km])
elseif or(Distance_Units=="m",Distance_Units=="m")
    X_plot=X_calibration{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
else
    X_plot=X_calibration{ PathNum,1}.(Distance_Metric).*1000;                % Length between source and measurement location (L is in meters [m])
end

field_data=errorbar(X_plot,C_E,X_calibration{PathNum,1}.SEM,"o","MarkerSize",8,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Field data');%"s", "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','Field data';    "o","MarkerSize",8,"MarkerFaceColor",imola(1,:),'DisplayName','Field data','HandleVisibility','on'
%--------------------------------------------------------------------------
C_L_AVG=[];
idx=0;
for m=1:length(X_calibration{PathNum,1}.Location)
    % if ismember(RastersLocations(m),X_calibration{PathNum,1}.Location)
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
    X_plot=[0:0.001:X_calibration{ PathNum,1}.(Distance_Metric)(end)];
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