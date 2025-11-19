Import_Summary_Table_Modified
addpath(genpath('..'))

Lithology="ARENITES"; %"METABASALTS" "ARENITES"

T=T_Summary(T_Summary.Lithology==Lithology,:); 

if Lithology=="METABASALTS"
    LithologyColor=[19, 141, 117]./255;
elseif Lithology=="ARENITES"
    LithologyColor=[192, 57, 43]./255;
end
Colors=[LithologyColor;[165, 105, 189]./255;[243, 156, 18]./255];

FigNum=length(findobj('type','figure'));
fig_1=figure(FigNum+1);
hold on;
ax_1=gca;
set(ax_1, 'XDir','reverse')

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
%--------------------------------------------------------------------------
%                   LINE WIDTH
linewidth=2;
%--------------------------------------------------------------------------
%                   COLORS
set(gcf,'Color','white')
% colororder(imola)
% ColOrder=imola(size(Parameter,2));
%--------------------------------------------------------------------------

if Lithology=="METABASALTS"
    plot_max_Dist=6.0; % Distance in km of the most upstream location plotted
    plot_min_Dist=3.5; % Distance in km of the most downstream location plotted
elseif Lithology=="ARENITES"
    plot_max_Dist=6.5; % Distance in km of the most upstream location plotted
    plot_min_Dist=4.5; % Distance in km of the most downstream location plotted
end

set(ax_1,"XLim",[plot_min_Dist plot_max_Dist])

Tab_Dit_idx=find(T.Distance>=plot_min_Dist,1,'last');

e_IRn=errorbar(ax_1,T.Distance(1:Tab_Dit_idx),T.NormCirc(1:Tab_Dit_idx),T.NormCirc_StDev(1:Tab_Dit_idx)./sqrt(T.SampleSize(1:Tab_Dit_idx)),"-o","Color",Colors(1,:),"MarkerSize",6,"MarkerEdgeColor","blue","MarkerFaceColor",Colors(1,:),'DisplayName',sprintf('IR_{n}'));

e_IR=errorbar(ax_1,T.Distance(1:Tab_Dit_idx),T.Circularity(1:Tab_Dit_idx),T.Circularity_StDev(1:Tab_Dit_idx)./sqrt(T.SampleSize(1:Tab_Dit_idx)),"--o","Color",Colors(2,:),"MarkerSize",6,"MarkerEdgeColor","blue","MarkerFaceColor",Colors(2,:),'DisplayName','IR');

e_Rnd=errorbar(ax_1,T.Distance(1:Tab_Dit_idx),T.Roundness(1:Tab_Dit_idx),T.Roundness_StDev(1:Tab_Dit_idx)./sqrt(T.SampleSize(1:Tab_Dit_idx)),"-.o","Color",Colors(3,:),"MarkerSize",6,"MarkerEdgeColor","blue","MarkerFaceColor",Colors(3,:),'DisplayName','Roundness');


ax_1.YAxis.TickLabelFormat = '%.2f';
ax_1.XAxis.TickLabelFormat = '%.2f';

xlabel({'$\leftarrow$ Distance from outlet [km]'},'fontsize',FontSize, 'FontWeight', 'bold');%\Leftarrow
ylabel('Shape parameter [-]','fontsize',FontSize, 'FontWeight', 'bold');
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
fileName = ['Shape_Parameter_Selection_Paper_', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = ['Shape_Parameter_Selection_Paper_', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,['Shape_Parameter_Selection_Paper_',datestr(now,'mmmm dd yyyy hhMM'),'.fig']);
%--------------------------------------------------------------------------
MatFileName=['Shape_Parameter_Selection_Paper_',datestr(now,'mmmm dd yyyy hhMM'),'.mat'];
save(MatFileName);



% %--------------------------------------------------------------------------
% fig_1=figure()
% % =gcf;
% ax_1=gca;
% set(ax_1, 'XDir','reverse')
% left=0;                         %in cm
% bottom=0;                       %in cm
% width=8.89;                     %in cm 8.89 for IEEE journals
% hw_ratio = 0.65;                % feel free to play with this ratio
% height=4;                       %in cm
% positionVect = [left bottom width height];  %Size and position of figure
% set(fig,'Units','centimeters')
% set(fig, 'PaperPosition', [positionVect]) %Set paper total dimension
% 
% set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
% leg=legend({'Roundness', '$IR$', '$IR_{n}$'})
% set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% 
% ax_1.FontSize=17
% addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
% 
% 
% xlabel({'$\leftarrow$ Distance from outlet [km]'},  'FontWeight', 'bold')
% ylabel('Shape parameter [-]',  'FontWeight', 'bold')
% 
% leg=legend({'Roundness', '$IR$', '$IR_{n}$'})
% 
% Xvalues=[0:1000:6000];
% set(ax_1,'XTick',Xvalues)
% 
% Xvalues_new=string;
% for i=1:size(Xvalues,2)
%     Xvalues_new(i)=sprintf('%d.0', Xvalues(i)/1000);
% end
% set(ax_1,'XTickLAbel',Xvalues_new)
% set(gcf,'Color','white')
% fileformat='.png';%'.eps' '.PDF'
% filename = ['Shape_Parameter_Selection_Paper ', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format
% export_fig (filename, '-nocrop')%'-m3',
