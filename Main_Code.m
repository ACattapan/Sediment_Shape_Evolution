%--------------------------------------------------------------------------
% This code was developed by 
%       Alessandro Cattapan (August 2024)
%       IHE DElft, Institute for Water Education
%       TU Delft
%       a.cattapan@un-ihe.org
%       a.cattapan@tudelft.nl
%   If you use this code in your work, please cite:
%
%--------------------------------------------------------------------------
%                       INPUTS DEFINITION
%--------------------------------------------------------------------------
% This code requires two sets of input:
% - a table with morphometric properties of each particle present in each
% image and its associated properties: lithology, origin type and location ID
% - a set of images (raster layers produced by GIS) representing the distance 
% of each pixel from the basin outlet
%--------------------------------------------------------------------------
%                       MODEL DESCRIPTION
%--------------------------------------------------------------------------
% THE CODE PERFORMS A SERIES OF TASKS:
% 1 IMPORT THE TABLE CONTAINING SEDIMENT DATA, CONVERT FROM #PIXELS TO 
%   DIMENSIONAL UNITS, REMOVE PARTICLES WITH SIZE SMALLER THAN A THRESHOLD 
%   AND COMPUTE SAMPLES SIZES FOR EACH ORIGIN TYPE, LITHOLOGY AND LOCATION
%   THIS CODE PRIDUCES: T,Lithologies,Origins,Locations,LocDist,SampleSize
%
% 2 IMPORT THE IMAGES OF TRAVEL DISTANCES FOR EACH SAMPLING LOCATION AND 
%   SAVE RESULTS TO OUTPUT FILE PATH
%   THIS CODE PRODUCES: RastersData,RastersLocations
%
% 3 IMPORT THE TABLE CONTAINING THE SERIES OF LOCATIONS WHERE SAMPLES WERE
%   TAKEN FOR EACH LITHOLOGY, ORGANIZED IN SETS OF PATHS
%   THIS CODE PRODUCES: Paths
%
% 4 ORGANIZE RESULTS IN A TABLE (T) AND SUBSTITUTE EMPTY VALUES WITH NANs
%
% Additional actions that it might perform if requested include: 
% 5 IMPORT IMAGES PROPERTIES FROM THE "Spreadsheet_Input" WORKBOOK
% 6 IDENTIFY COLUMNS ATTRIBUTES AND PASTE THEM INTO THE OUTPUT VARIABLES: T
%   AND PropertINPUT
% 7 WRITE THE FULL TABLE (T) INTO THE "Spreadsheet_Output" WORKBOOK AND 
%   SAVE THE MAT FILE
% 8 IF THE FLAG Run_TD_Model IS SET TO "y", THE MODEL FITS THE PROPOSED
%   THEORETICAL EQUATION FOR SEDIMENT SHAPE AS A FUNCTION OF TRAVEL DISTANCE
%   AND PRODUCES THE OPTIMIZATION PARAMETERS ASSOCIATED. THIS STEP CAN BE
%   PERFORMED IN MULTIPLE WAYS, DEPENDING ON THE PARAMETERS AND FLAGS SELECTED.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% clc
tic
disp('Running Main Analysis File')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 1 SET CURRENT, INPUT, OUTPUT AND MAT FOLDERS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                       LITHOLOGY AND PATH TO ANALYZE
Lithology="METABASALTS";                % Alternative values are: "ARENITES", "METABASALTS", "W100"

if Lithology=="ARENITES"
    PathNum=1;                           % 1 for ARENITES; 2 or 3 for METABASALTS
    Last_Location_Calibration=15;        % This should be 15 for ARENITES and 22 for METABASALTS
    LithologyColor=[192, 57, 43]./255;
elseif Lithology=="METABASALTS"
    PathNum=2;                           % 1 for ARENITES; 2 or 3 for METABASALTS
    Last_Location_Calibration=22;        % This should be 15 for ARENITES and 22 for METABASALTS
    LithologyColor=[19, 141, 117]./255;
elseif Lithology=="W100"
    PathNum=1;                           % 1 for ARENITES; 2 or 3 for METABASALTS
    Last_Location_Calibration=23;        % This should be 15 for ARENITES and 22 for METABASALTS
end

HillsNetw_Flag="N";                      % Flag to use Hillslope sources rasters (H) or Network sources rasters (N)
CaliVali_Flag="C";                       % Flag for Calibration (C) or Validation (V) modality
%--------------------------------------------------------------------------
%                       SHAPE PARAMETER AND DISTANCE FOR MODEL INPUT
%--------------------------------------------------------------------------
Parameter="NormCirc";                   % Alternative values are: "Circularity", "NormCirc", "Roundness", "elongation",...
Distance_Metric="SourceDist";           % Alternative values are: "SourceDist", or "TravDist"
DistanceWeigthing="N";                  % Use or not exponential distance weighting "Y" or "N"
%--------------------------------------------------------------------------
% INPUT AND OUTPUT FOLDERS, FILES AND SETTINGS
%-----------------------INPUTS---------------------------------------------
% Define the paths where the functions used by the program are stored
% N.B. Distances in the Excel file and in the rasters should be in coherent
% units: [m].
restoredefaultpath;
Functions_folder="REQUIRED_FUNCTIONS";
Functions_folder_path=fullfile(pwd,Functions_folder);
addpath(genpath(Functions_folder_path));
% addpath 'D:\02_DATA\CODES\MATLAB Codes\EXPORT_FIG'
% addpath 'C:\Program Files\gs\gs10.05.1'
addpath(genpath('..'))
% color_palette="imola";
% color_palette_path=sprintf("D:\\02_DATA\\CODES\\MATLAB Codes\\ScientificColourMaps8\\%s\\%s.mat",color_palette,color_palette);
% load(color_palette_path)
% colormap(imola);
colormap(parula);
% Define current folder and the folders for the inputs: sediment properties
% data and Travel distance images data
Current_folder_path=pwd;
Input_folder="INPUTS";
Input_file="INPUT_DATASET_2025.xlsx";
Input_folder_path=fullfile(pwd,Input_folder);
Input_file_path=fullfile(Input_folder_path,Input_file);
if HillsNetw_Flag=="H"
    Input_folder_Travel_Dist_Images=sprintf("Travel_Distances_%s",Lithology);
elseif HillsNetw_Flag=="N"
    Input_folder_Travel_Dist_Images=sprintf("Travel_Distances_Network_%s",Lithology);
end
Input_folder_Travel_Dist_Images_path=fullfile(pwd,Input_folder_Travel_Dist_Images);
%--------------------------------------------------------------------------
%                       PATHS FILE AND PATH NUMBER
Paths_Input_file=sprintf("%s_PATHS.xlsx",Lithology);
Paths_Input_file_path=fullfile(Input_folder_path,Paths_Input_file);
%-----------------------OUTPUTS--------------------------------------------
Output_folder="PAPER_RESULTS";%PAPER_RESULTS RESULTS
Output_folder_path=fullfile(pwd,Output_folder);
Output_file=sprintf("DATA_IMPORT_%s_%s.mat",Lithology,Parameter);
Output_file_path=fullfile(Output_folder_path,Output_file);
Output_file_Travel_Dist=sprintf("TRAVEL_DISTANCES_%s_%s.mat",Lithology,Parameter);
Output_file_Travel_Dist_path=fullfile(Output_folder_path,Output_file_Travel_Dist);
Stats_Overview_file=sprintf("STATS_OVERVIEW_%s_%s.mat",Lithology,Parameter);
Stats_Overview_path=fullfile(Output_folder_path,Stats_Overview_file);
%-----------------------SAVE MAT FILES FOR MODEL RESULTS-------------------
Circularity_Fit_Model_file=sprintf("CIRC_FIT_MODEL_AVG_%s_%s_PATH_%s_%s_PAPER.mat",Parameter,Lithology,string(PathNum),Distance_Metric);
Circularity_Fit_Model_file_matfile=fullfile(Output_folder_path,Circularity_Fit_Model_file);

Circularity_Fit_Model_Quantiles_file=sprintf("CIRC_FIT_QUANT_MODEL_AVG_C0_AVG_%s_%s_PATH_%s_6876.mat",Lithology,Parameter,string(PathNum));
Circularity_Fit_Model_Quantiles_file_matfile=fullfile(Output_folder_path,Circularity_Fit_Model_Quantiles_file);

Circularity_Fit_Model_Quantiles_common_Ka_file=sprintf("CIRC_FIT_QUANT_COMM_Ka_MODEL_AVG_C0_AVG_%s_%s_PATH_%s_6876.mat",Lithology,Parameter,string(PathNum));
Circularity_Fit_Model_Quantiles_common_Ka_file_matfile=fullfile(Output_folder_path,Circularity_Fit_Model_Quantiles_common_Ka_file);

Circularity_Fit_Model_Quantiles_common_function_file=sprintf("CIRC_FIT_QUANT_COMM_FUNCTION_MODEL_AVG_C0_AVG_%s_%s_PATH_%s_6876.mat",Lithology,Parameter,string(PathNum));
Circularity_Fit_Model_Quantiles_common_function_file_matfile=fullfile(Output_folder_path,Circularity_Fit_Model_Quantiles_common_function_file);

TD_Model_file_Variab_Parameters=sprintf("TD_MODEL_AVERAGES_VARIAB_PARAM_%s_%s_PATH_%s_6876.mat",Lithology,Parameter,string(PathNum));
TD_Model_Variab_Parameters_matfile=fullfile(Output_folder_path,TD_Model_file_Variab_Parameters);
%--------------------------------------------------------------------------
%                       SETTINGS
%--------------------------------------------------------------------------
%                       MODEL PARAMETERS SET-UP
C_max=1.0;
mu_max=1.0;
% C0=0.5:0.01:0.95;%[0.40:0.01:0.95];   % C0=0.48; % C0=0.75;   Initial value of the Parameter considered, at the source [-]
a=0.1:0.01:10;                          % Exponent of the elliptical function relating relative mass loss and the shape Parameter chosen [-]
Ka_exp=0.0:0.01:10;                     % Exponent of 10^-( ) in the coefficient Ka controlling Sternberg's Law
Ka=10.^-Ka_exp;                         % Coefficient Ka controlling Sternberg's Law the units are 1/Distance_Unit
lambda=[50:10:2000]./1000;              % parameter of exponential probability weighting function in [km]
Distance_Units="km";                    % Controls in which units distances are reported. It can be "km", or "m".
%--------------------------------------------------------------------------
%      CHOOSE WHETHER TO FIT SHAPE PARAMETER FOR TRAVEL DISTANCE MODEL AND
%      WHETHER C VALUES LOWER THAN C0 HAVE TO BE DELTED (Y) ON NOT (N)
Fitting_Param="N";
CLowerC0="N";
%--------------------------------------------------------------------------
Detect_Thresh=7;                        % Size in mm of the smallest grain considered.
Plotting_flag="N";                      % Controls whether plots of models results will be produced or not
Plot_Type="pdf";                        % Can be "pdf" for probability density function or "cdf" for cumulative density function
Fig_number=1;                           % Number of the first figure produced
Form='tif';                             % Format of rasters containing Travel Distances (GIS produced)
FitTest="AD";                           % FitTest: "KS" or "AD" for Kolmogorov-Smirnov and Anderson-Darling respectively
alpha=0.05;                             % Significance level of statistic tests performed
Quantiles_Used=[5,10:10:90,95];         % Quantiles used in the discretization of the shape parameters distribution
Run_TD_Model="N";                       % Flag to run the series of model listed at the end (indidual ones can be omitted by commenting them)
%--------------------------------------------------------------------------
%                       MODEL BEGIN
%--------------------------------------------------------------------------
% IMPORT SEDIMENT DATA, FILTER THEM BASED ON SMALLER SIZE
% N.B.: in the table T, distances are in meters [m] but the function
% Stats_Extract_Path transforms them to km [km] and that is how they are
% stored in X. Each model converts them to meters, if needed.
[T,Lithologies,Origins,Locations,LocDist,SampleSize]=Data_Import_V1(Input_file_path,Output_file_path,Detect_Thresh,Plotting_flag,Fig_number);
%--------------------------------------------------------------------------
% THIS PART IS FOR CALIBRATION
% IMPORT RASTES DATA OF TRAVEL DISTANCE FOR EACH PATH, REMOVE NEGATIVE AND NaN VALUES
[RastersData, RastersLocations]=Travel_Distance_Extract_func(Form,Current_folder_path,Input_folder_Travel_Dist_Images_path,Output_file_Travel_Dist_path,LocDist,Plotting_flag,Plot_Type);
%--------------------------------------------------------------------------
% THE FOLLOWING LOOP CONVERTS RASTERS DATA DISTANCES TO Km. COMMENT IF DISTANCES ARE
% NEEDED IN METERS.
if or(Distance_Units=="km",Distance_Units=="Km")
    for i=1:numel(RastersData)
        RastersData{i,1}=RastersData{i,1}./1000;
    end
end
%--------------------------------------------------------------------------
%               IMPORT PATHS MATRIX FOR A GIVEN LITHOLOGY
%--------------------------------------------------------------------------
Paths = ImportPaths(Paths_Input_file_path);
%--------------------------------------------------------------------------
%               EXTRACT INDIVIDUAL PATHS FROM THE MATRIX
%--------------------------------------------------------------------------
[LocPath, OrigPath, NumPaths]=Path_Extract(Origins,Paths);
%--------------------------------------------------------------------------
% REMOVE THE SAMPLING LOCATIONS THAT ARE DRAINING MULTIPLE
% SOURCES OF THE LITHOLOGY CONSIDERED AND USE ONLY THOSE FOR CALIBRATION
%--------------------------------------------------------------------------
OrigPath_Calibration=OrigPath;
LocPath_calibration=LocPath;
OrigPath_Calibration{PathNum,1}(find(LocPath_calibration{PathNum,1}(:)>Last_Location_Calibration,1):end)=[];
LocPath_calibration{PathNum,1}(find(LocPath_calibration{PathNum,1}(:)>Last_Location_Calibration,1):end)=[];
Num_Locations_calibration=max(size(LocPath_calibration{PathNum}));
RastersLocations_calibration=zeros(Num_Locations_calibration,1);
RastersData_calibration=cell(Num_Locations_calibration,1);
for i=1:max(size(LocPath_calibration{PathNum}))
    RastersLocations_calibration(i,1)=RastersLocations(RastersLocations==LocPath_calibration{PathNum,1}(i));
    RastersData_calibration{i,1}=RastersData{RastersLocations==LocPath_calibration{PathNum,1}(i)};
    RastersData_calibration{i,1}(isnan(RastersData_calibration{i}))=[];                         % Remove NaN values
    RastersData_calibration{i,1}=sort(RastersData_calibration{i});                              % Sort Travel Distances for Location m
end
%--------------------------------------------------------------------------
%               COMPUTE PARAMETER STATISTICS FOR THE CHOSEN PATH
%--------------------------------------------------------------------------
X=cell(NumPaths,1);
FreqDistr=cell(NumPaths,1);
for k=1:NumPaths
    %----------------------------------------------------------------------
    X{k,1}=Stats_Extract_Path(T,Parameter,alpha,Quantiles_Used,Lithology,OrigPath{k},LocPath{k},RastersLocations,RastersData);
    FreqDistr{k,1}=Freq_Distr_Extract(T,Parameter,FitTest,alpha,Lithology,OrigPath{k},LocPath{k},Plotting_flag,Plot_Type);
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
% REMOVE FROM THE CALIBRATION STRUCTURE, THE LOCATIONS WITH MULTIPLE
% SOURCES DOWNSTREAM
X_calibration=X;
X_calibration{PathNum,1}(X_calibration{PathNum,1}.Location>Last_Location_Calibration,:)=[];
save(Stats_Overview_path);
%--------------------------------------------------------------------------
%               TRAVEL DISTANCE MODEL
%--------------------------------------------------------------------------
if or(Run_TD_Model=="Y",Run_TD_Model=="y")

            % Circularity_Fit_Model(X_calibration,Parameter,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file_matfile,mu_max,C_max,C0,a,Ka)
            % Circularity_Fit_Model_fixed_C0(X_calibration,Parameter,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file_matfile,mu_max,C_max,a,Ka)
    Circularity_Fit_Model_fixed_C0_Distances_Samples(Lithology,LithologyColor,X_calibration,Parameter,Distance_Units,Distance_Metric,DistanceWeigthing,RastersData, RastersLocations,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file_matfile,mu_max,C_max,a,Ka,lambda)

    Quantiles_Names="q_"+string(Quantiles_Used);
    Parameters_used=Quantiles_Names;
            % Circularity_Fit_Model_Quantiles(X,Parameters_used,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_Quantiles_file_matfile,mu_max,C_max,C0,a,Ka)%Ka_exp
            % Circularity_Fit_Model_Quantiles_Ka(X_calibration,Parameters_used,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_Quantiles_file_matfile,mu_max,C_max,C0,a,Ka)%Ka_exp
            % Circularity_Fit_Model_Quantiles_Ka_V2(X_calibration,Parameters_used,Distance_Units,Distance_Metric,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_Quantiles_file_matfile,mu_max,C_max,C0,a,Ka)%Ka_exp
                      % (Lithology,LithologyColor,X_calibration,Parameter      ,Distance_Units,Distance_Metric,DistanceWeigthing,RastersData, RastersLocations,Fitting_Param,FitTest,CLowerC0,PathNum,Circularity_Fit_Model_file_matfile,mu_max,C_max,a,Ka,lambda) 
            % Circularity_Fit_Model_Quantiles_Common_Ka(X_calibration,Parameters_used,Distance_Units,Distance_Metric                                                                               ,PathNum,Circularity_Fit_Model_Quantiles_common_Ka_file_matfile,mu_max,C_max,C0,a,Ka)%Ka_exp, Fitting_Param,FitTest,CLowerC0,
    Circularity_Fit_Model_Quantiles_Common_Ka_fixed_C0(Lithology,LithologyColor,X_calibration,Parameters_used,Distance_Units,Distance_Metric                                                                               ,PathNum,Circularity_Fit_Model_Quantiles_common_Ka_file_matfile,mu_max,C_max,a,Ka)
            % Circularity_Fit_Model_Quantiles_Common_function(X_calibration,Parameters_used,Distance_Units,Distance_Metric,PathNum,Circularity_Fit_Model_Quantiles_common_function_file_matfile,mu_max,C_max,C0,a,Ka)
    %----------------------------------------------------------------------
    % Trav_Dist_Model_V03(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_Variab_Parameters_matfile);
    %----------------------------------------------------------------------
    %           LOAD RESULTS PRODUCED BY THE OPTIMIZATION ALGORITHM
    % load(Circularity_Fit_Model_file_matfile)
    % load(Circularity_Fit_Model_Quantiles_file_matfile)
    % load(Circularity_Fit_Model_Quantiles_common_Ka_file_matfile)
    % load(Circularity_Fit_Model_Quantiles_common_function_file_matfile)
    % load(TD_Model_Variab_Parameters_matfile)

    %----------------------------------------------------------------------
%     scatter(TDe_avg(TDe_avg>0),TD_M_avg(TD_M_avg>0),"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
%     title('Comparison of DEM derived and modelled Travel Distances');
%     xlabel('Travel Distances obtained from DEM [m]');
%     ylabel('Travel Distances obtained from sediment morphometry model [m]');
toc
end