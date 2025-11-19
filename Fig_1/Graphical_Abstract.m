%--------------------------------------------------------------------------
%                           IMPORT DATA
%--------------------------------------------------------------------------
OptParamAVG = importfile("Opt_Param_AVG.xlsx", "Sheet1", [2, Inf]);

load('T.mat');
%--------------------------------------------------------------------------
Lithologies=OptParamAVG.Lithology;
C0=OptParamAVG.C0;
a=OptParamAVG.a;
ka=OptParamAVG.ka;
%--------------------------------------------------------------------------
%-------------------FIGURE SET UP------------------------------------------
addpath(genpath('..'))
fig_1=figure;hold on
ax_1=subplot(3,1,1);hold on
ax_2=subplot(3,1,2);hold on
ax_3=subplot(3,1,3);hold on
hold on
set(fig_1,'DefaultFigureColor',[1 1 1])         %Figure background color
%--------------------------------------------------------------------------
%                   DIMENSIONS
set(0,'Units','centimeters')
left=0;                         %in cm
bottom=0;                       %in cm
width=18.00;                     %in cm 8.89 for IEEE journals
hw_ratio = 0.65;                % feel free to play with this ratio
height=6;                       %in cm

positionVect = [left bottom width height];  %Size and position of figure

set(fig_1, 'PaperPosition', [positionVect]) %Set paper total dimension

%--------------------------------------------------------------------------
%                   FONT NAME
% listfonts
fontname = 'times';                         %'Helvetica'
set(0, 'DefaultAxesFontName', fontname);     
set(0, 'DefaultTextFontName', fontname);    
set(fig_1,'defaulttextinterpreter','latex')     %Allows us to use LaTeX maths notation


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
linewidth=1.5;
%--------------------------------------------------------------------------
%                   COLOR PALETTE
PS.DRed1 = [205, 97, 85]./255;
PS.DRed2 = [192, 57, 43]./255;
PS.DRed3 = [169, 50, 38]./255;
PS.DRed4 = [146, 43, 33]./255;
PS.DRed5 = [123, 36, 28]./255;

PS.DGreen1 = [69, 179, 157]./255;
PS.DGreen2 = [19, 141, 117]./255;
PS.DGreen3 = [17, 122, 101]./255;
PS.DGreen4 = [14, 102, 85]./255;
PS.DGreen5 = [11, 83, 69]./255;

PS.Green1 = [82, 190, 128]./255;
PS.Green2 = [39, 174, 96]./255;
PS.Green3 = [34, 153, 84]./255;
PS.Green4 = [30, 132, 73]./255;
PS.Green5 = [20, 90, 50]./255;

PS.Orange1 = [245, 176, 65]./255;
PS.Orange2 = [243, 156, 18]./255;
PS.Orange3 = [214, 137, 16]./255;
PS.Orange4 = [185, 119, 14]./255;
PS.Orange5 = [156, 100, 12]./255;

PS.Blue1 = [133, 193, 233]./255;
PS.Blue2 = [93, 173, 226]./255;
PS.Blue3 = [52, 152, 219]./255;
PS.Blue4 = [40, 116, 166]./255;
PS.Blue5 = [27, 79, 114]./255;
ColorPalette=PS;
%--------------------------------------------------------------------------
%                   ACTUAL PLOTTING STARTS HERE
MetabasaltsColor = PS.DGreen2;
ArenitesColor =PS.DRed2;
MarkerSize = 30;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
C_max=1;
mu_max=1;
mu_plot=linspace(0,1,1000000);
Locations={[3,4,5,10,13,14,15];[6,7,11,13,14,15,17,18,21,22,23,24]}; %[3,4,5,10,13,14,15,17,18,21,22,23,24]
for i=1:length(Lithologies)
    Str=Lithologies(i);
    Lith_case=strcat(upper(Str{1}(1)),lower(Str{1}(2:end)));
    HandleVis='on';
    if Lithologies(i)=='ARENITES'
        Color=ArenitesColor;
    elseif Lithologies(i)=='METABASALTS'
        Color=MetabasaltsColor;
    end
    for j=1:length(Locations{i,1})
        
        LocDist{i,1}(j)=unique(T.Distance(T.Location==Locations{i,1}(j)));
        SourceDist{i,1}(j)=(LocDist{i,1}(1)-unique(T.Distance(T.Location==Locations{i,1}(j))))/1000;

        sz=length(T.NormCirc(T.Lithology==Lithologies(i)&T.Location==Locations{i,1}(j)));

            Y_plot_SSize{i,j}=sz;
            Y_plot_Data{i,j}=T.NormCirc(T.Lithology==Lithologies(i)&T.Location==Locations{i,1}(j));
            Y_plot_AVG{i,j}=mean(Y_plot_Data{i,j});
            Y_plot_STD{i,j}=std(Y_plot_Data{i,j});
            if sz>0
                Y_plot_STD_Err{i,j}=Y_plot_STD{i,j}./sqrt(Y_plot_SSize{i,j});
            else
                Y_plot_STD_Err{i,j}=0;
            end
    end
    %----------------------------------------------------------------------
    L_plot{i}=1/ka(i).*log(1./(1-mu_plot));

    C_plot{i}=C_of_mu_func(mu_plot,mu_max,C0(i),C_max,a(i));
    %----------------------------------------------------------------------
    plot_1{i}=plot(ax_1,L_plot{i},mu_plot);
    set(plot_1{i}, 'LineStyle', '-','LineWidth',linewidth,'color',Color,'Marker', 'none','DisplayName',sprintf('%s ',Lith_case),'HandleVisibility',HandleVis); %$\mu$ of $\mathcal{L}$ for 
    plot_2{i}=plot(ax_2,mu_plot,C_plot{i});
    set(plot_2{i}, 'LineStyle', '-','LineWidth',linewidth,'color',Color,'Marker', 'none','DisplayName',[sprintf('IR_{n} of '),'\mu',sprintf(' for %s ',Lith_case)],'HandleVisibility',HandleVis);
    plot_3{i}=plot(ax_3,L_plot{i},C_plot{i});
    set(plot_3{i}, 'LineStyle', '-','LineWidth',linewidth,'color',Color,'Marker', 'none','DisplayName',[sprintf('IR_{n} of L for %s ',Lith_case)],'HandleVisibility',HandleVis); %[sprintf('IR_{n} of '),'\mathcal{L}',sprintf(' for %s',Lith_case)]
    X_coord=cell2mat(SourceDist(i,1));
    Y_coord=cell2mat(Y_plot_AVG(i,:));
    Y_err=cell2mat(Y_plot_STD_Err(i,:));
    error3{i}=errorbar(ax_3,X_coord,Y_coord,Y_err,'vertical',"o","LineStyle","none", "MarkerEdgeColor",Color,"MarkerFaceColor",Color,'HandleVisibility','off');
end

ax_1.YAxis.TickLabelFormat = '%.1f';
ax_1.XAxis.TickLabelFormat = '%.1f';
ax_1.XLim=[0 500];
ax_1.XLabel.String =({'Travel distance L [km]'}); %,  'FontWeight', 'bold')
ax_1.YLabel.String =({'$\mu$ [-]'}); %,  'FontWeight', 'bold')
lg1=legend(ax_1,'location','northeast');
set(lg1,'Box','off');

ax_2.YAxis.TickLabelFormat = '%.1f';
ax_2.XAxis.TickLabelFormat = '%.1f';
ax_2.XLabel.String =({'Relative mass-loss $\mu$ [-]'}); %,  'FontWeight', 'bold')
ax_2.YLabel.String =({'$IR_{n}$ [-]'}); %,  'FontWeight', 'bold')
lg2=legend(ax_2,'location','southeast');
set(lg2,'Box','off');

ax_3.YAxis.TickLabelFormat = '%.1f';
ax_3.XAxis.TickLabelFormat = '%.1f';
ax_3.XLim=[0 4.5];
ax_3.XLabel.String =({'Travel distance L [km]'}); %,  'FontWeight', 'bold')
ax_3.YLabel.String =({'$IR_{n}$ [-]'}); %,  'FontWeight', 'bold')
lg3=legend(ax_3,'location','southeast');
set(lg3,'Box','off');
set(gcf,'Color','white');
%--------------------------------------------------------------------------
%                   SAVE GRAPH TO FILE
%Now ready for export
% print(hfig,fname,'-dpng','-painters')
fileformat='.png'; %'.eps' '.PDF'
fileName = ['Graphical_Abstract_', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format

export_fig (fileName, '-nocrop'); %'-m3',

savefig(fig_1,['Graphical_Abstract_', datestr(now,'mmmm dd yyyy hhMM'), '.fig']);

dataName=['Graphical_Abstract_',datestr(now,'mmmm dd yyyy hhMM'),'.mat'];
save(dataName)






%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function OptParamAVG = importfile(workbookFile, sheetName, dataLines)
%% Input handling

% If no sheet is specified, read from Sheet1
if nargin == 1 || isempty(sheetName)
    sheetName = "Sheet1";
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["Lithology", "C0", "a", "ka"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Lithology", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Lithology", "EmptyFieldRule", "auto");

% Import the data
OptParamAVG = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    OptParamAVG = [OptParamAVG; tb]; %#ok<AGROW>
end

end