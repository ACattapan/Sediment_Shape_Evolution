load('T.mat');
load('X.mat');
Lithology_Normality="ARENITES";     % "ARENITES","METABASALTS", "W100"
A_Label="Source";
B_Label="Non-Source";
alpha=0.05;
A=T.NormCirc(T.Lithology==Lithology_Normality&T.Origin=="SOURCE");
B=T.NormCirc(T.Lithology==Lithology_Normality&T.Origin~="SOURCE");
Stats_Overview_2_Vars(A,A_Label,B,B_Label,Lithology_Normality,alpha)