function [T,Lithologies,Origins,Locations,LocDist,SampleSize]=Data_Import(Input_file_path,Output_file_path,Detect_Thresh,Plottting_flag,Fig_number)
%--------------------------------------------------------------------------
% THIS CODE IMPORTS THE DATA FROM THE DEFINED FOLDER AND FILE,
% IT ALSO FILTERS DATA AND COMPUTES SAMPLES SIZES FOR EACH ORIGIN TYPE, LITHOLOGY AND LOCATION
%--------------------------------------------------------------------------
clc
tic
disp('Running Data Import')
%--------------------------------------------------------------------------
% N.B. Data in the imoprt spreadsheet are supposed to be in pixels. These data 
% will be transpoformed in the desired unit for all output based on the average 
% scaling expressed in the spreadsheet.
% Available options for output units are: 
%                   "px" = pixels
%                   "m"  = meters
%                   "cm" = centimeters
%                   "mm" = millimiters
%                   "inch" = inches
%                   "ft" = feet
Output_data_units="mm";
%--------------------------------------------------------------------------
%               DEFINITION OF INPUTS AND OUTPUT FOLDERS NAMES
%--------------------------------------------------------------------------
% Add all sub-folders of the current folder to MATLAB paths
% addpath(genpath(pwd));
%--------------------------------------------------------------------------
%                           DATA IMPORT
%--------------------------------------------------------------------------
% T0=ImportDataset(Input_file_path);
T0=readtable(Input_file_path,'PreserveVariableNames',true);
T=ScaleANDClean_Vars(T0,Detect_Thresh,Plottting_flag,Fig_number);
% h=(T.a-T.b).^2./(T.a+T.b).^2;
% EllPerim=1/2*pi*(T.a+T.b).*(1+(3.*h./(10+(4-3.*h).^(1/2))));
% EllArea=1/4*pi.*T.a.*T.b;
% T.TheorCirc=4*pi.*EllArea./(EllPerim.^2);
% T.NormCirc=T.Circularity./T.TheorCirc;
% T = movevars(T,"TheorCirc",'After',"Circularity");
% T = movevars(T,"NormCirc",'After',"TheorCirc");
%--------------------------------------------------------------------------
%                           DATA ANALYSIS
%--------------------------------------------------------------------------
% SAMPLING LOCATIONS ARE NAMED FROM UPSTREAM TO DOWNSTEAM ACCORDING TO
% THEIR DISTANCE FROM THE BASIN OUTLET;
% Possible Lithologies are :    ARENITES=1, METABASALTS=2,  W100=3,     WERFER=4;
% Possible Origines are    :    CREEK=1,    RIVER=2,        SOURCE=3;
% Metab_SampleSize=SampleSize(:,2,[6 7 11 13 14 15 17 18 21 22 23 24]);
% Code to extract the number of particles collected in locations belonging
% to path 2 for Metabasalts
Lithologies = string(unique(T.Lithology(:)));
Origins     = string(unique(T.Origin(:)));

% DEFINE SAMPLING LOCATIONS AND THEIR DISTANCES FROM BASIN OUTLET
Locations=unique(T.Location);
LocDist=sort(unique(T.Distance),'descend');
%--------------------------------------------------------------------------
% DEFINITION OF SAMPLE SIZE FOR EACH LOCATION AND FOR EACH LITHOLOGY
% AND ORIGIN
%--------------------------------------------------------------------------
SampleSize=[];
for i=1:size(Origins,1)
    for j=1:size(Lithologies,1)
        for k=1:size(Locations,1)
            SampleSize(i,j,k)=size(T.Pebble_Id(T.Lithology==Lithologies(j)&T.Location==Locations(k)&T.Origin==Origins(i)),1);
        end
    end
end
%--------------------------------------------------------------------------
save(Output_file_path);
toc