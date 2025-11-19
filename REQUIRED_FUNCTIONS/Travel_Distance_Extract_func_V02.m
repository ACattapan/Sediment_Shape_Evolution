function [RastersData,RastersLocations]=Travel_Distance_Extract_func_V02(Form,Current_folder_path,Input_folder_path,Output_file_path,LocDist,Plot_Flag,Plot_Type)
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
tic
disp('Running Travel Distance Extraction')
%--------------------------------------------------------------------------
%                           DATA IMPORT
%--------------------------------------------------------------------------
% Import all tif files from the folder Distance_Rasters
cd (Input_folder_path)
tifFiles = dir(Format);
numfiles=length(dir(Format));
for i=1:numfiles
    Filename(i,1)=string(tifFiles(i).name);
    FilenameDouble(i,1)=str2double(erase(Filename(i,1),strcat(".",Form)));
end
Filename=[];
FilenameDouble=sort(FilenameDouble,'ascend');
Filename=strcat(string(FilenameDouble),".",Form);
% FUNCTION FOUND ON:
% https://it.mathworks.com/matlabcentral/answers/410795-how-to-read-series-of-images-from-folder-in-specific-order
% fun = @(x) tifFiles(x).name == Filename(i);
% tf2 = arrayfun(fun, 1:numel(tifFiles));
% index2 = find(tf2);
for i=1:numfiles
    %    Filename(i,1)=strcat(string(i),'.tif');
%     Filename(i,1)=string(tifFiles(i).name);
    fun = @(x) tifFiles(x).name == Filename(i);
    indx(i,1)=find(arrayfun(fun, 1:numel(tifFiles)));
end
RastersLocations = Filename;
RastersLocations = str2double(erase(RastersLocations,strcat(".",Form)));
%--------------------------------------------------------------------------
RastersData = cell(1, numfiles);
%--------------------------------------------------------------------------
% Store all distance rasters within the RastersData cell array, change all
% negative values to NaN
% and print their cdf as histograms, if needed
for k = 1:numfiles
  RastersData{k} = imread(tifFiles(indx(k)).name);
  RastersData{k}(RastersData{k}(:,:)<0)=NaN;
  
  if LocDist(RastersLocations(k))<min(RastersData{k},[],'all')
      RastersData{k}=RastersData{k}-LocDist(RastersLocations(k));  
  else
      RastersData{k}=RastersData{k}-min(RastersData{k},[],'all');
  end
  
  if or(Plot_Flag=="Y",Plot_Flag=="y")
      figure(k)
      nbins = 50;
      Nmin = 30;
      [N,edges] = histcounts(RastersData{k},nbins);
      histogram(RastersData{k}/1000,nbins,'normalization',Plot_Type);
      title(sprintf('Travel Distances Distribution at Location '+string(k)));
      xlabel("Distribution of distances [km]");
  end
end
save(Output_file_path);
cd(Current_folder_path);
toc