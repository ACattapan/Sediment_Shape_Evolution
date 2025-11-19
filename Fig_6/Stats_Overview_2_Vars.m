function Stats_Overview_2_Vars(A,A_Label,B,B_Label,Lithology,alpha)
% 
%--------------------------------------------------------------------------
arguments
    A (:,1) double
    A_Label string
    B (:,1) double
    B_Label string
    Lithology string
    alpha(1,1) double
end
%--------------------------------------------------------------------------
tic
disp('Running Stats_Overview_2_Vars')
%--------------------------------------------------------------------------
% addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
% addpath 'C:\Program Files\gs\gs10.05.1'
addpath(genpath('..'))
%--------------------------------------------------------------------------
if Lithology=="ARENITES"
    LithologyColor=[192, 57, 43]./255;
elseif Lithology=="METABASALTS"
    LithologyColor=[19, 141, 117]./255;
elseif Lithology=="W100"
    LithologyColor=[162, 162, 162]./255;
end
%--------------------------------------------------------------------------
% Remove empty values and NaNs
A=rmmissing(A);
% A(isnan(A))=[];
B=rmmissing(B);
%--------------------------------------------------------------------------
% Descriptive statistics
A_SplSz=length(A);
B_SplSz=length(B);
A_mean=mean(A);
B_mean=mean(B);
A_std=std(A);
B_std=std(B);
A_stErr=A_std/sqrt(A_SplSz);
B_stErr=B_std/sqrt(B_SplSz);
A_skw=skewness(A);
B_skw=skewness(B);
%--------------------------------------------------------------------------
% Test for Normality

[h_A,p_A,adstat_A,cv_A] = adtest(A,'Alpha',alpha);
[h_B,p_B,adstat_B,cv_B] = adtest(B,'Alpha',alpha);

%--------------------------------------------------------------------------

[p_Wilcoxon,h_Wilcoxon] = ranksum(A,B);

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
%--------------------------------------------------------------------------
h1=normplot(A);
h2=normplot(B);
h1(1).Marker="o";
% h1(1).Color="blue";
h1(1).Color=[192, 57, 43]./255;
h2(1).Marker="o";
% h2(1).Color=LithologyColor;
h2(1).Color=[19, 141, 117]./255;

h1(2).LineStyle = '-.';
h1(2).LineWidth = 2;
h2(2).LineStyle = '-.';
h2(2).LineWidth = 2;

%--------------------------------------------------------------------------
ax_1.XAxis.TickLabelFormat = '%.2f';

xlabel('$IR_{n}$ [-]','fontsize',FontSize);

legend([h1(1) h2(1)],{A_Label,B_Label},'Location','southeast','Interpreter','latex')
legend('boxon')
ax_1.Title.String='';
set(gcf,'Color','white')
%--------------------------------------------------------------------------
hold off;
%--------------------------------------------------------------------------
%                   SAVING fig_1 AS PDF, PNG AND AS .fig
fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Normality_Plot_Source_NonSource_%s',Lithology,fileformat); 
export_fig (fileName, '-nocrop')%'-m3',

fileformat='.png';%'.eps' '.PDF', '.png'
fileName = sprintf('%s_Normality_Plot_Source_NonSource_%s',Lithology,fileformat); 
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,sprintf('%s_Normality_Plot_Source_NonSource_%s',Lithology,'.fig'));
%--------------------------------------------------------------------------
% Save results
FileName=[sprintf('%s_Statistical_Test_%s',Lithology,'.mat')];
save(FileName)