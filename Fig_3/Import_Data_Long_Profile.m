%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 20);

% Specify sheet and range
opts.Sheet = "SUMMARY_DATA";
opts.DataRange = "A2:T37";

% Specify column names and types
opts.VariableNames = ["Location", "Lithology", "Origin", "Distancekm", "SampleSize", "bAverage", "bStDev", "St_Error", "elongationAverage", "elongationStDev", "AreaAverage", "AreaStDev", "RoundnessAverage", "RoundnessStDev", "CircularityAverage", "CircularityStDev", "NormalizedCircularityAverage", "NormalizedCircularityStDev", "Distance", "ZMeters"];
opts.VariableTypes = ["double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Lithology", "Origin"], "EmptyFieldRule", "auto");

% Import the data
Summary_Data = readtable("Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Profile_Path_1";
opts.DataRange = "A2:E1425";

% Specify column names and types
opts.VariableNames = ["Vertex", "XMeters", "YMeters", "ZMeters", "Distance"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Profile_Path_1 = readtable("D:\08_PhD\10_PAPERS\Circularity_evolution_model\MATERIAL\FIGURES\Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Profile_Path_3";
opts.DataRange = "A2:E1210";

% Specify column names and types
opts.VariableNames = ["Vertex", "XMeters", "YMeters", "ZMeters", "Distance"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Profile_Path_3 = readtable("D:\08_PhD\10_PAPERS\Circularity_evolution_model\MATERIAL\FIGURES\Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Profile_Path_6";
opts.DataRange = "A2:E1251";

% Specify column names and types
opts.VariableNames = ["Vertex", "XMeters", "YMeters", "ZMeters", "Distance"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Profile_Path_6 = readtable("D:\08_PhD\10_PAPERS\Circularity_evolution_model\MATERIAL\FIGURES\Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Profile_Path_8";
opts.DataRange = "A2:E1077";

% Specify column names and types
opts.VariableNames = ["Vertex", "XMeters", "YMeters", "ZMeters", "Distance"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Profile_Path_8 = readtable("D:\08_PhD\10_PAPERS\Circularity_evolution_model\MATERIAL\FIGURES\Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Profile_Path_19";
opts.DataRange = "A2:E590";

% Specify column names and types
opts.VariableNames = ["Vertex", "XMeters", "YMeters", "ZMeters", "Distance"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
Profile_Path_19 = readtable("D:\08_PhD\10_PAPERS\Circularity_evolution_model\MATERIAL\FIGURES\Data_Long_Profile.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts