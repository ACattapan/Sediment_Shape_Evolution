%--------------------------------------------------------------------------
Import_Data_Long_Profile; %Import all the data from the spreadsheet
%--------------------------------------------------------------------------
addpath(genpath('..'))
%--------------------------------------------------------------------------
% Remove locations 1, 2 and 12 and renumber them from 1 to 21
Summary_Data_Plot=Summary_Data;
Summary_Data_Plot(Summary_Data_Plot.Location==1,:)=[];
Summary_Data_Plot(Summary_Data_Plot.Location==2,:)=[];
Summary_Data_Plot(Summary_Data_Plot.Location==12,:)=[];
Locations_List=Summary_Data_Plot.Location;
[Locations_unique,ia,ic]=unique(Locations_List);
New_Locations=[1:size(Locations_unique,1)]';
Locations_List_new=New_Locations(ic);
Locations_List_new=reshape(Locations_List_new,size(Locations_List));
Summary_Data_Plot.Location=Locations_List_new;
%--------------------------------------------------------------------------
StreamSedDist=Summary_Data_Plot.Distance(Summary_Data_Plot.Lithology=='MIXED SAMPLE'&Summary_Data_Plot.Origin=="RIVER");

StreamSedSize=Summary_Data_Plot.bAverage(Summary_Data_Plot.Lithology=='MIXED SAMPLE'&Summary_Data_Plot.Origin=="RIVER");
StreamSedSizeErr=Summary_Data_Plot.bStDev(Summary_Data_Plot.Lithology=='MIXED SAMPLE'&Summary_Data_Plot.Origin=="RIVER")./sqrt(Summary_Data_Plot.SampleSize(Summary_Data_Plot.Lithology=='MIXED SAMPLE'&Summary_Data_Plot.Origin=="RIVER"));

%--------------------------------------------------------------------------
%-------------------FIGURE SET UP------------------------------------------
fig_1=figure;
ax_1=gca;
hold on
set(fig_1,'DefaultFigureColor',[1 1 1])         %Figure background color
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
fontname = 'times';                         %'Helvetica'
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
StreamColor = PS.Blue3;
LocationColor = PS.Orange2;
MetabasaltsColor = PS.DGreen2;
ArenitesColor =PS.DRed2;
MarkerSize = 30;

% Plotting Longitudinal Profiles of different paths
% Long_Profile=plot(Profile_Path_1.Distance,Profile_Path_1.ZMeters,'LineWidth',linewidth,'color',StreamColor,'DisplayName','$Longitudinal profile$');
Long_Profile=plot(Profile_Path_1.Distance,Profile_Path_1.ZMeters);
Path_3=plot(Profile_Path_3.Distance,Profile_Path_3.ZMeters);
Path_6=plot(Profile_Path_6.Distance,Profile_Path_6.ZMeters);
Path_8=plot(Profile_Path_8.Distance,Profile_Path_8.ZMeters);
Path_19=plot(Profile_Path_19.Distance,Profile_Path_19.ZMeters);
% Plotting Sampling Locations
SamplLocDots=scatter(Summary_Data_Plot.Distance(ia),Summary_Data_Plot.ZMeters(ia));%'MarkerFaceColor'

% Changing the stacking order of elements
% (https://blogs.mathworks.com/graphics/2014/11/04/sortmethod/)
ax_1.Children=ax_1.Children([1 6 2 3 4 5]);

%--------------------------------------------------------------------------
%             ADJUST LINE PROPERTIES (FUNCTIONAL + AESTHETICS)
%
set(Long_Profile, 'LineStyle', '-','LineWidth',linewidth, 'color',StreamColor,'Marker', 'none','DisplayName','Longitudinal profile','HandleVisibility','off');
set(Path_3, 'LineStyle', '-','LineWidth',linewidth,'color',ArenitesColor,'Marker', 'none','DisplayName','Arenites path');
set(Path_6, 'LineStyle', '-','LineWidth',linewidth,'color',MetabasaltsColor,'Marker', 'none','DisplayName','Metabasalts path');
set(Path_8, 'LineStyle', '-','LineWidth',linewidth,'color',MetabasaltsColor,'Marker', 'none','HandleVisibility','off');
set(Path_19, 'LineStyle', '-','LineWidth',linewidth,'color',ArenitesColor,'Marker', 'none','HandleVisibility','off');
set(SamplLocDots,'SizeData',MarkerSize,'MarkerEdgeColor', LocationColor, 'MarkerFaceColor',LocationColor,'DisplayName','Sampling Locations');%,'color',LocationColor,'filled'
% set(textHandle,'Vert','bottom', 'Horiz','left', 'FontSize',12)

set(ax_1, 'XDir','reverse')

%--------------------------------------------------------------------------
%                   AXES AND LABELS
 
axis([0 9000 500 1500]) 
 
set(ax_1,'YTick',[500:250:1500])

set(ax_1,'YTickLAbel',string([500:250:1500]))

Xvalues=[0:1000:9000];
set(ax_1,'XTick',Xvalues)

Xvalues_new=string;
for i=1:size(Xvalues,2)
    Xvalues_new(i)=sprintf('%d.0', Xvalues(i)/1000);
end
set(ax_1,'XTickLAbel',Xvalues_new)

xlabel({'$\leftarrow$ Distance from outlet [km]'},  'FontWeight', 'bold') % Note cell matrices for line breaks
ylabel('Elevation [m.a.sl.]',  'FontWeight', 'bold')

grid on
%--------------------------------------------------------------------------
%                   SECONDARY AXIS PLOT
colororder({'k'});
yyaxis(ax_1,'right')

errHandle=errorbar(StreamSedDist,StreamSedSize,StreamSedSizeErr,'vertical',"o","LineStyle","none");

set(errHandle,"MarkerSize",5,'Color',"blue",...
    "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'DisplayName','D_{50}')

ylabel('$D_{50} [mm]$',  'FontWeight', 'bold') %Intermediate size b [mm]

%--------------------------------------------------------------------------

lgd=legend('FontSize',11,'Location','southwest');
legend('boxoff')
set(gcf,'Color','white')
%--------------------------------------------------------------------------
%                   SAVE GRAPH TO FILE
fileformat='.png';%'.eps' '.PDF'
fileName = ['Longitudinal_Profile_', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format

export_fig (fileName, '-nocrop')%'-m3',

fileformat='.PDF';%'.eps' '.PDF', '.png'
fileName = ['Longitudinal_Profile_', datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format
export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,['Longitudinal_Profile_', datestr(now,'mmmm dd yyyy hhMM'), '.fig']);

dataName=['Plot_Long_Profile_Image_',datestr(now,'mmmm dd yyyy hhMM'), '.mat'];
save(dataName)