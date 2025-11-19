function [RastersData,RastersLocations]=Travel_Distance_Extract_func(Form,Current_folder_path,Input_folder_path,Output_file_path,LocDist,Plot_Flag,Plot_Type)
%--------------------------------------------------------------------------
% THIS CODE IMPORTS .tif IMAGES OF TRAVEL DISTANCES FOR EACH SAMPLING LOCATION
% AND SAVES RESULTS TO OUTPUT FILE PATH
%--------------------------------------------------------------------------
arguments
    Form string
    Current_folder_path string
    Input_folder_path string

    Output_file_path string
    LocDist
    Plot_Flag 
    Plot_Type
end
%---------CHECK WHETHER THE FORMAT FOR INPUT IMAGES WAS INTRODUCED---------
if isempty(Form)
    Form='tif';
end
Format=strcat('*.',Form);
%size(char(Form),2)
match="."+digitsPattern(size(char(Form),2));
tic
disp('Running Travel Distance Extraction')
%--------------------------------------------------------------------------
%                           DATA IMPORT
%--------------------------------------------------------------------------
% List all files with the specified format from the Input folder path
cd (Input_folder_path)
tifFiles = dir(Format);
numfiles=length(dir(Format));
if numfiles==0
    error('There are no files with extension '+Form+' in the folder '+Input_folder_path)
    return
end
% Files are sorted in alphabetic order inside tiffiles. In order to sorted,
% file names must be converted to string first, the file format is erased
% and their name is converted to double to sort them in ascending order.
% Finally their name is converted to string back and the file format is
% attached back.
for i=1:numfiles
    Filename(i,1)=string(tifFiles(i).name);
    FilenameDouble(i,1)=str2double(extractBefore(Filename(i,1),"."));
end
Filename=[];
FilenameDouble=sort(FilenameDouble,'ascend');
Filename=strcat(string(FilenameDouble),".",Form);
%--------------------------------------------------------------------------
% Create an anonymous function that compares the 
% 
for i=1:numfiles
    fun = @(x) tifFiles(x).name == Filename(i);
    indx(i,1)=find(arrayfun(fun, 1:numfiles));
end
RastersLocations = Filename;
RastersLocations = str2double(extractBefore(RastersLocations,"."));
%--------------------------------------------------------------------------
% Store all distance rasters within the RastersData cell array, change all
% negative values to NaN
% 
% Print their histograms in pdf or cdf format according to Plot_Type
%--------------------------------------------------------------------------
% Pre allocate the cell with the proper size
RastersData = cell(numfiles,1);
for k = 1:numfiles
  RastersData{k} = imread(tifFiles(indx(k)).name);
  RastersData{k}(RastersData{k}(:,:)<0)=NaN;
  
  if LocDist(RastersLocations(k))<min(RastersData{k},[],'all')
      % COMMENT THE LINES WHERE THE DIFFERENCE IS TAKEN IF TRAVEL DISTANCE
      % RASTERS HAVE BEEN CREATED WITH  RESPECT TO THE MEASURING LOCATION
      % AND NOT TO THE CATCHMENT OUTLET
      % RastersData{k}=RastersData{k}-LocDist(RastersLocations(k));
      RastersData{k}(isnan(RastersData{k}))=[];
  else
      % RastersData{k}=RastersData{k}-min(RastersData{k},[],'all');
      RastersData{k}(isnan(RastersData{k}))=[];
  end
  
  % if or(Plot_Flag=="Y",Plot_Flag=="y")
      % h =  findobj('type','figure');
      % n = length(h);
      % figure(n+1)
      % nbins = 50;
      % Nmin = 30;
      % [N,edges] = histcounts(RastersData{k},nbins);
      % histogram(RastersData{k}/1000,nbins,'normalization',Plot_Type);
  %     title(sprintf('Travel Distances Distribution at Location '+string(RastersLocations(k))));
  %     xlabel("Distribution of distances [km]");
  % end
end

if or(Plot_Flag=="Y",Plot_Flag=="y")
    figure
    ColOrder=parula;
    hold on;
    Legend_Texts=[];
    for k = 1:numfiles
        ColIndex=ceil(k*(size(ColOrder,1)/numfiles));
        histogram(RastersData{k}/1000,'normalization',Plot_Type,'FaceColor',ColOrder(ColIndex,:));%,'FaceColor',ColOrder(k+25)
        Legend_Texts=[Legend_Texts,'Location '+string(RastersLocations(k))];
    end
    title(sprintf('Travel Distances Distribution at Multiple Locations'));
    xlabel("Distribution of distances [km]");
    legend(Legend_Texts);
end
save(Output_file_path);
cd(Current_folder_path);
toc