function Trav_Dist_Model_Averages_Const_Parameters(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
% 
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    Fitting_Param (1,1) string
    FitTest (1,1) string
    CLowerC0 (1,1) string
    C0 (:,:) double
    a (:,:) double
    Ka_exp (:,:) double
    Lithology (:,:) string
    Origins (:,:) string
    Paths (:,:) double
    PathNum double
    RastersData cell
    RastersLocations (:,:) double
    alpha (1,1) double
    TD_Model_file (:,:) string
    % OrigPath (:,:) string
    % LocPath (1,:) %{mustBeNumeric(LocPath)}
end
%--------------------------------------------------------------------------
tic
disp('Running Travel Distance Model Constant Ka')
% This function uses a conceptual model to derive sediment travel distances
% distributions from distributions of circularity based on two assumptions:
% curvature flow attrition (modelled by an elliptical relationship between
% relative mass loss and travel distance) and
% Sternberg Law relationship between relative mass loss and travel distance
matfile=TD_Model_file;
%--------------------------------------------------------------------------
%                      EXTRACT PATHS FOR A GIVEN LITHOLOGY
[LocPath, OrigPath, NumPaths]=Path_Extract(Origins,Paths);
% Consider the list of locations that form a certain path: PathNum
disp(['Path number: ',num2str(PathNum)])
Locations=LocPath{PathNum};
disp(['Locations: ',num2str(Locations)])
%--------------------------------------------------------------------------
% CONVERT TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM GENERATED) IN A
% COLUMN VECTOR, REMOVE NaNs AND DEFINE THE SAMPLE SIZE
for m=1:size(Locations,2)
    RasterData_TDtest=RastersData{1,RastersLocations==Locations(m)}(:); % Select the RasterData of Travel Distances of Location m
    RasterData_TDtest(isnan(RasterData_TDtest))=[];                     % Remove NaN values
    RasterData_TDtest=sort(RasterData_TDtest);                          % Sort Travel Distances for Location m
    Emp_SampSz=size(RasterData_TDtest,1);                               % Number of pixels in the RasterData of Location m
    %----------------------------------------------------------------------
    % COMPUTE THE AVERAGE AND STANDARD DEVIATION OF TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM
    % GENERATED) IN [m]
    %----------------------------------------------------------------------
    TD_E_Stats(1,m)=mean(RasterData_TDtest,'all','omitnan');            % Mean of empirical Travel Distance from RasterData for location m
    TD_E_Stats(2,m)=std(RasterData_TDtest);                             % St.Dev. of empirical Travel Distance from RasterData for location m
end
%--------------------------------------------------------------------------
%                      DEFINE RANGES FOR MODEL PARAMETERS VALUES
mu_max=1.0;
C_max=1.0;
% C0=[0.45:0.01:0.99];
% % C0=0.48;
% a=1.0:0.1:8;
% Ka_exp=1.0:0.1:10;
Ka=10.^-(Ka_exp);
ParamSpace=size(C0,2)*size(a,2)*size(Ka,2);
%--------------------------------------------------------------------------
%                      INITIALIZE VARIABLES FOR EACH LOCATION
RMSE_Avg_M=zeros(size(C0,2),size(a,2),size(Ka,2));          % Root Mean Square Error for the Average of the Modelled Travel Distances
RMSE_M_Avg=zeros(size(C0,2),size(a,2),size(Ka,2));          % Root Mean Square Error for the Model of the Average fo Travel Distances
% I_Avg_M=zeros(1,size(Locations,2));         % Linear Index for the minimum RMSA for the Average of Modelled Travel Distances
% I_M_Avg=zeros(1,size(Locations,2));         % Linear Index for the minimum RMSA for the Model of the Average Travel Distance
% Opt_Param_Avg=zeros(4,size(Locations,2));   % Optimum parameters for the Average of Modelled Travel Distances (the fourth row contains the minimum RMSE)
% Opt_Param_M=zeros(4,size(Locations,2));     % Optimum parameters for the Model of the Average Travel Distance (the fourth row contains the minimum RMSE)
%--------------------------------------------------------------------------
Locations_Analysis=Locations;
Locations_Analysis(1:3)=[];
Locations_Analysis(end-1:end)=[];
Origin_Analysis=OrigPath{PathNum};

Origin_Analysis(1:3)=[];
Origin_Analysis(end-1:end)=[];
TD_E_Stats(:,1:3)=[];
TD_E_Stats(:,end-1:end)=[];
%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
%           ACTUAL MODEL TEST FOR EACH PARAMETER
%--------------------------------------------------------------------------
%           INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE CYCLE
%--------------------------------------------------------------------------
Iter_N=0;
for i=1:size(C0,2)
    for j=1:size(a,2)
        for k=1:size(Ka,2)
            Iter_N=Iter_N+1;
            % disp(['Iteration: ',num2str(Iter_N)])
            % Delete variables for each location to avoid rewriting on them
            % and indexing errors
            TD_Avg_M=[];
            TD_M_Avg=[];
            ERR_Avg_M=[];
            ERR_M_Avg=[];
            % Cycle through all locations of the path considered
            for m=1:size(Locations_Analysis,2)
                %----------------------------------------------------------
                % SHAPE PARAMETER USED TO ESTIMATE TRAVEL DISTANCE
                Param=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
                Sample_Size=size(Param,1); % Size of the sample in location m
                if Sample_Size>0
                    C=sort(Param);
                    %------------------------------------------------------
                    % CHOOSE WHETHER C VALUES LOWER THAN C0 HAVE TO BE DELTED (Y) ON NOT (N)
                    if or(CLowerC0=="Y",CLowerC0=="y")
                        C(C<C0(i))=[];
                    end
                    C_Avg=mean(C);
                    %------------------------------------------------------
                    mu=mu_func(C,mu_max,C0(i),C_max,a(j));
                    mu_Avg=mu_func(C_Avg,mu_max,C0(i),C_max,a(j));
                    Ls=Length_func(mu,Ka(k));
                    Ls_Avg=Length_func(mu_Avg,Ka(k));
                    TD_Avg_M(m)=mean(Ls,'all');
                    TD_M_Avg(m)=Ls_Avg;
                    ERR_Avg_M(m)= TD_E_Stats(1,m)-TD_Avg_M(m);
                    ERR_M_Avg(m)= TD_E_Stats(1,m)-TD_M_Avg(m);
                end
                clear C C_Avg mu mu_Avg Ls Ls_Avg
            end
            RMSE_Avg_M(i,j,k)= rmse(TD_Avg_M(1:m),TD_E_Stats(1,1:m));
            RMSE_M_Avg(i,j,k)= rmse(TD_M_Avg(1:m),TD_E_Stats(1,1:m));
        end
    end
end
%--------------------------------------------------------------------------
% clear TD_Avg_M TD_M_Avg ERR_Avg_M ERR_M_Avg
%--------------------------------------------------------------------------
%           OPTIMUM PARAMETERS FOR AVERAGES
%--------------------------------------------------------------------------
[min_RMSE_Avg_M,I_Avg_M]=min(RMSE_Avg_M,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_AM,I2_AM,I3_AM]=ind2sub([size(RMSE_Avg_M,1) size(RMSE_Avg_M,2) size(RMSE_Avg_M,3)],I_Avg_M);
Opt_Param_Avg(1)=C0(I1_AM);
Opt_Param_Avg(2)=a(I2_AM);
Opt_Param_Avg(3)=Ka(I3_AM);
Opt_Param_Avg(4)=RMSE_Avg_M(I1_AM,I2_AM,I3_AM);
%--------------------------------------------------------------------------
[min_RMSE_M_Avg,I_M_Avg]=min(RMSE_M_Avg,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_MA,I2_MA,I3_MA]=ind2sub([size(RMSE_M_Avg,1) size(RMSE_M_Avg,2) size(RMSE_M_Avg,3)],I_M_Avg);
Opt_Param_M(1)=C0(I1_MA);
Opt_Param_M(2)=a(I2_MA);
Opt_Param_M(3)=Ka(I3_MA);
Opt_Param_M(4)=RMSE_M_Avg(I1_MA,I2_MA,I3_MA);
%--------------------------------------------------------------------------
%           OPTIMUM ESTIMATES OF TRAVEL DISTANCES FOR EACH LOCATION
%--------------------------------------------------------------------------
for m=1:size(Locations_Analysis,2)
    Param=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
    Sample_Size=size(Param,1); % Size of the sample in location m
    if Sample_Size>0
        C{m}=sort(Param);
        C_Avg{m}=mean(Param);
    end
    mu{m}=mu_func(C{m},mu_max,C0(I1_AM),C_max,a(I2_AM));
    mu_Avg{m}=mu_func(C_Avg{m},mu_max,C0(I1_MA),C_max,a(I2_MA));
    Ls{m}=Length_func(mu{m},Ka(I3_AM));
    Ls_Avg{m}=Length_func(mu_Avg{m},Ka(I3_MA));
    TD_Avg(m)=mean(Ls{m}(:),'all');
    TD_M(m)=Ls_Avg{m};
end
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
Fig_Avg_Const_Param=figure (FigNum+1)
scatter(TD_E_Stats(1,:)./1000,TD_M./1000,100,'filled')
hold on;
scatter(TD_E_Stats(1,:)./1000,TD_Avg./1000,100,'+')
axeslim=max([xlim;ylim],[],'all');
xlim([0 axeslim]);
ylim([0 axeslim]);
plot(xlim, ylim, '--k')
legend('Model of the TD from the Average Circularity','Average of Modelled TD from Circularity')
title('Model results Hyp: const model parameters');
xlabel('Empirical Travel Distance from DTM [km]');
ylabel('Modelled Travel Distance [km]');
FigNameformat='Scatter_plot_for_Averages_of_%s_Path_Number_%d_Constant_Parameters.fig';
saveas(Fig_Avg_Const_Param,sprintf(FigNameformat,Lithology,PathNum))
%--------------------------------------------------------------------------

A_Avg(:,:)=RMSE_Avg_M(I1_AM,:,:);
[minimumA_Avg, IndMinA_Avg]=min(A_Avg,[],'all','linear');
[colj pagek]=ind2sub([size(A_Avg,1) size(A_Avg,2)],IndMinA_Avg);

FigNum=length(findobj('type','figure'));
Fig_RMSE_Avg_M_Const_Param=figure (FigNum+1)
surf(A_Avg)
hold on;
scatter3(pagek,colj,A_Avg(colj,pagek),100,'r','filled')
zlim([0 10000])
title('RMSE for Average of Model results; Hyp: const model parameters');
xlabel('C0');
ylabel('a');
FigNameRMSE_Avg_M_Const_Param='RMSE_Surf_for_Model_of_Averages_of_%s_Path_Number_%d_Constant_Parameters.fig';
saveas(Fig_RMSE_Avg_M_Const_Param,sprintf(FigNameRMSE_Avg_M_Const_Param,Lithology,PathNum))

%--------------------------------------------------------------------------

A_M(:,:)=RMSE_M_Avg(I1_MA,:,:);
[minimumA_M, IndMinA_M]=min(A_M,[],'all','linear');
[colj pagek]=ind2sub([size(A_M,1) size(A_M,2)],IndMinA_M);

FigNum=length(findobj('type','figure'));
Fig_RMSE_M_Avg_Const_Param=figure (FigNum+1)
surf(A_M)
hold on;
scatter3(pagek,colj,A_M(colj,pagek),100,'r','filled')
zlim([0 10000])
title('RMSE for Model results of Averages; Hyp: const model parameters');
xlabel('C0');
ylabel('a');
FigNameRMSE_M_Avg_Const_Param='RMSE_Surf_for_Model_of_Averages_of_%s_Path_Number_%d_Constant_Parameters.fig';
saveas(Fig_RMSE_M_Avg_Const_Param,sprintf(FigNameRMSE_M_Avg_Const_Param,Lithology,PathNum))
%--------------------------------------------------------------------------

save(matfile);
toc