%--------------------------------------------------------------------------
% For AGU (American Geophysical Union) journals, 
% figures should generally be sized between 50 mm and 170 mm in width. 
% One-column figures typically range from 50 to 85 mm wide. 
% Figures should be submitted at a high resolution, with 300 dpi for images and photos, 
% and 600 dpi for line art. Acceptable file formats include EPS, TIFF, JPEG, or PNG. 
% Key Details:
% 
%     Width: Figures can be sized between 50 mm and 170 mm wide. 
% 
% One-Column Figures: A single-column figure should be between 50 and 85 mm wide. 
% Resolution:
% 
%     Images and photos: 300 dpi 
% 
% Line art: 600 dpi
%--------------------------------------------------------------------------
%                   LOADING THE DATA
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
load Plot_cdf_Data.mat
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
addpath(genpath('..'))
fig_1=figure;
ax_1=gca;
hold on
set(fig_1,'DefaultFigureColor',[1 1 1])         %Figure background color
set(gcf,'Color','white')
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
LenLith=length(Lithologies);
ExcludeLocations=[1;2;12];
for i=1%1:length(Lithologies)
    Str=Lithologies(i);
    Lith_case=strcat(upper(Str{1}(1)),lower(Str{1}(2:end)));
    if Lithologies(i)=='ARENITES'
        Color=ArenitesColor;
    elseif Lithologies(i)=='METABASALTS'
        Color=MetabasaltsColor;
    end

    Source_Locations{i,1}=unique(T.Location(T.Lithology==Lithologies(i)&T.Origin==Origins(3)&~ismember(T.Location,ExcludeLocations)));

    NumSources=length(Source_Locations{i,1});

    for j=1:NumSources
        if j==1
            HandleVis='on';
        else
            HandleVis='off';
        end
        Sourceplot{j,1}=cdfplot(T.NormCirc(T.Lithology==Lithologies(i)&T.Origin==Origins(3)&T.Location==Source_Locations{i,1}(j)));

        set(Sourceplot{j,1}, 'LineStyle', '-','LineWidth',linewidth,'color',Color,'Marker', 'none','DisplayName',sprintf('Source %s',Lith_case),'HandleVisibility',HandleVis);
    end
    NonSource_Locations{i,1}=unique(T.Location(T.Lithology==Lithologies(i)&T.Origin~=Origins(3)));
    NumNonSources=length(NonSource_Locations{i,1});
    for j=1:NumNonSources
        if j==1
            HandleVis='on';
        else
            HandleVis='off';
        end
        Nonsourceplot{j,1}=cdfplot(T.NormCirc(T.Lithology==Lithologies(i)&T.Origin~=Origins(3)&T.Location==NonSource_Locations{i,1}(j)));

        set(Nonsourceplot{j,1}, 'LineStyle', '-','LineWidth',linewidth,'color',[184 184 184]./255,'Marker', 'none','DisplayName',sprintf('Non-Source %s',Lith_case),'HandleVisibility',HandleVis);
    end
end
legend('location','northwest')
legend boxoff
title('')

ax_1.YAxis.TickLabelFormat = '%.1f';
ax_1.XAxis.TickLabelFormat = '%.1f';

xlabel({'$IR_{n}$ [-]'},  'FontWeight', 'bold') % Note cell matrices for line breaks
ylabel('cdf($IR_{n}$) [-]',  'FontWeight', 'bold')

%--------------------------------------------------------------------------
%                   INSERT LEGEND
% legend({ '$y$','$z$'},'Interpreter','latex','Location','best');
%--------------------------------------------------------------------------
%                   SAVE THE FIGURE IN A PLOT
% Save a color *.eps file
% print('figue_of_my_data.eps','-depsc');
% print('figue_of_my_data.eps','-djpgc');
fileformat='.png';%'.eps' '.PDF'
fileName = [sprintf('IRn_cdf_ %s',Lith_case), datestr(now,'mmmm dd yyyy hhMM'), fileformat]; %Descriptive name timestamp and .png file format

export_fig (fileName, '-nocrop')%'-m3',

savefig(fig_1,[sprintf('IRn_cdf_ %s',Lith_case), datestr(now,'mmmm dd yyyy hhMM'), '.fig'])



































% set(0,'DefaultAxesFontSize',10); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
% set(0,'DefaultFigureColor','w')
% set(0,'defaulttextinterpreter','latex') %Allows us to use LaTeX maths notation
% set(0, 'DefaultAxesFontName', 'times');
% 
% fig=figure  %Let's make a simple time series plot of notional data
% set(gcf, 'Units','centimeters')
% 
% %Set figure total dimension
% set(gcf, 'Position',[0 0 8.89 4]) %Absolute print dimensions of figure. 8.89cm is essential here as it is the linewidth of a column in IEEE format
% %Height can be adusted as suits, but try and be consistent amongst figures for neatness
% %[pos_from_left, pos_from_bottom, fig_width, fig_height]
% 
% 
% 
% stairs(Time, Voltage, 'LineWidth', 1.5, 'Color', [190,174,212]/255); %Plot as paired data, so we're explicity stipulating the time index
% %Use a stairstep plot so as not to linearly interpolate 
% 
% axis([0 48 18 22]) %Note symetric vertical axis +- 2kV around the nominal 20kV level
% 
% set(gca,'YTick',[18, 19.2, 20, 21.25, 22]) %Now impose sensible tickmark locations
% %Let's pretend that the undervoltage limit is 19.2 kV, and the overvoltage is 21.25 kV
% 
% set(gca,'YTickLAbel',{'18 kV', 'V_{min}', '20 kV', 'V_{max}', '22 kV'}) %Now put in informative labels at these tickmarks
% 
% %Now sort out the horizontal axes: it needs to be shown in wallclock units
% 
% set(gca,'XTick',[0, 8, 16, 24, 32, 40, 48]) %Now impose sensible tickmark locations
% 
% set(gca,'XTickLAbel',{'00:00', '08:00', '16:00','00:00', '08:00', '16:00', '00:00'}) %consistent tick interval
% 
% 
% %Set size and position of axes plotting area within figure dimensions
% %It is nice to keep the vertical axes aligned for multiple figures, so be consistent with the horizontal positioning of axes 
% set(gca, 'Units','centimeters')
% set(gca, 'Position',[2 0.9 6.5 2.9]) %This is the relative positioning of the axes within the frame. 
% %[inset_from_left, inset_from_bottom, axes_width, axes_height]
% 
% box off %Removes the borders of the plot area
% 
% ylabel({'↑'; 'Bus 2'; 'Voltage';'(kV)'},  'FontWeight', 'bold') %Note cell matrices for line breaks
% 
% set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
% 
% 
% xlabel('Time →',  'FontWeight', 'bold') %Note use of unicode arrow for clarity
% 
% %Now ready for export
% 
% filename = ['VoltageTimeSeries ', datestr(now,'mmmm dd yyyy hhMM'), '.png'] %Descriptive name timestamp and .png file format
% 
% export_fig (filename, '-m3', '-nocrop')
