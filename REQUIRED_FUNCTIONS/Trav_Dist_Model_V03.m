function Travi_Dist_Model_V03(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
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
disp('Running Travel Distance Model')
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
%                      DEFINE RANGES FOR MODEL PARAMETERS VALUES
mu_max=1.0;
C_max=1.0;
% C0=[0.48:0.01:0.55];
% a=[0.6:0.1:10];
% Ka_exp=[1.0:0.25:7];
% C0=0.48:0.02:0.50;
% a=0.8:0.2:4;
% Ka_exp=1.0:0.25:6;
Ka=10.^-(Ka_exp);
% ParamSpace=size(C0,2)*size(a,2)*size(Ka,2);
%--------------------------------------------------------------------------
%                      INITIALIZE VARIABLES FOR EACH LOCATION
TDe_Stats=zeros(2,size(Locations,2));       % Average and St.Dev. of empirical Travel Distances
beta1=zeros(1,size(Locations,2));           % First parameter of Beta distribution
beta2=zeros(1,size(Locations,2));           % Second parameter of Beta distribution
h_beta=zeros(1,size(Locations,2));          % Hypothesis of Beta distribution fit
p_beta=zeros(1,size(Locations,2));          % p_value of Beta distribution fit
TD_M_avg=zeros(1,size(Locations,2));        % Average of modelled Travel Distances
RMSE_avg=zeros(1,size(Locations,2));        % Root Mean Square Error for the Average fo Travel Distances
I=zeros(1,size(Locations,2));               % Linear Index for the minimum RMSA for the Average of Travel Distances
Opt_Param_avg=zeros(4,size(Locations,2));   % Optimum parameters for the Average of Travel Distances (the fourth row contains the minimum RMSE)
RMSE=zeros(1,size(Locations,2));            % Minimum RMSE for the distribution of Travel Distance for each Location
I_distr=zeros(1,size(Locations,2));         % Linear Index for the minimum RMSA for the distribution of Travel Distances
Opt_Param_distr=zeros(4,size(Locations,2)); % Optimum parameters for the distribution of Travel Distances (the fourth row contains the minimum RMSE)
%--------------------------------------------------------------------------
%                      INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE
%                      CYCLE FOR EACH LOCATION
TDm_avg=zeros(size(C0,2),size(a,2),size(Ka,2));     % Average Travel Distance modelled for each parameter combination
ERR_avg=zeros(size(C0,2),size(a,2),size(Ka,2));     % Error function for the Average of modelled Travel Distance
RMSE_iter=zeros(size(C0,2),size(a,2),size(Ka,2));   % RMSE for Travel Distance distributions for each combination of parameters
%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
% Cycle through all locations of the path considered
for m=1:size(Locations,2)
    disp(['Location: ',num2str(Locations(m))])
    %----------------------------------------------------------------------
    %           ERASE VARIABLES AT EACH CYCLE TO SAVE MEMORY
    %----------------------------------------------------------------------
    RasterData_TDtest=[];
    Emp_SampSz=[];
    mu=[];
    Ls=[];
    TDm_avg=[];
    ERR_avg=[];
    RMSE_iter=[];
    %----------------------------------------------------------------------
    % CONVERT TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM GENERATED) IN A
    % COLUMN VECTOR, REMOVE NaNs AND DEFINE THE SAMPLE SIZE
    RasterData_TDtest=RastersData{1,RastersLocations==Locations(m)}(:); % Select the RasterData of Travel Distances of Location m
    RasterData_TDtest(isnan(RasterData_TDtest))=[];                     % Remove NaN values
    RasterData_TDtest=sort(RasterData_TDtest);                          % Sort Travel Distances for Location m
    Emp_SampSz=size(RasterData_TDtest,1);                                   % Number of pixels in the RasterData of Location m
    %----------------------------------------------------------------------
    % COMPUTE THE AVERAGE AND STANDARD DEVIATION OF TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM
    % GENERATED) IN [m] 
    %----------------------------------------------------------------------
    TDe_Stats(1,m)=mean(RasterData_TDtest,'all','omitnan');             % Mean of empirical Travel Distance from RasterData for location m
    TDe_Stats(2,m)=std(RasterData_TDtest);                              % St.Dev. of empirical Travel Distance from RasterData for location m
    %----------------------------------------------------------------------    
    Param=T.(Parameter)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m)); % Parameter values measured in location m
    Sample_Size=size(Param,1);                                          % Size of the sample in location m
    if Sample_Size>0
        %------------------------------------------------------------------
        % EXTRACT THE SHAPE PARAMETERS VALUES TO BE USED BY THE MODEL
        %------------------------------------------------------------------
        if or(Fitting_Param=="Y",Fitting_Param=="y")
            % [beta1(m),beta2(m),h_beta(m),p_beta(m)]=BetaFit(T,Parameter,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
            % C=betarnd(beta1(m),beta2(m),Emp_SampSz,1);
            FreqDistr{PathNum,1}=Freq_Distr_Extract(T,Parameter,FitTest,alpha,Lithology,OrigPath{PathNum},LocPath{PathNum},'N','pdf');
            C=betarnd(FreqDistr{PathNum,1}{1,Locations(m)}(1,1),FreqDistr{PathNum,1}{1,Locations(m)}(1,2),Emp_SampSz,1);
        else
            C=Param;
        end
        C=sort(C);
        %------------------------------------------------------------------
        %           ACTUAL MODEL TEST FOR EACH PARAMETER
        %------------------------------------------------------------------
        Iter_N=0;
        for i=1:size(C0,2)
            for j=1:size(a,2)
                for k=1:size(Ka,2)
                    Iter_N=Iter_N+1
                    mu=mu_func(C,mu_max,C0(i),C_max,a(j));
                    Ls=Length_func(mu,Ka(k));
                    TDm_avg(i,j,k)=mean(Ls,'all');
                    ERR_avg(i,j,k)= sqrt((TDe_Stats(1,m)-TDm_avg(i,j,k))^2);
                    RMSE_iter(i,j,k)=sqrt(mean((Ls-RasterData_TDtest).^2));
                end
            end
        end
        %------------------------------------------------------------------
        %           OPTIMUM PARAMETERS FOR AVERAGES
        %------------------------------------------------------------------
        [RMSE_avg(m),I(m)]=min(ERR_avg,[],'all','linear');
        % Convert linear indices I(m) into row,col indices basedon the
        % dimensions of the matrix R2
        [I1,I2,I3]=ind2sub([size(ERR_avg,1) size(ERR_avg,2) size(ERR_avg,3)],I(m));
        TD_M_avg(m)=TDm_avg(I1,I2,I3);
        Opt_Param_avg(1,m)=C0(I1);
        Opt_Param_avg(2,m)=a(I2);
        Opt_Param_avg(3,m)=Ka(I3);
        Opt_Param_avg(4,m)=ERR_avg(I1,I2,I3);
        %------------------------------------------------------------------
        %           OPTIMUM PARAMETERS FOR DISTRIBUTIONS
        %------------------------------------------------------------------
        [RMSE(m),I_distr(m)]=min(RMSE_iter,[],'all','linear');
        [I_d1,I_d2,I_d3]=ind2sub([size(RMSE_iter,1) size(RMSE_iter,2) size(RMSE_iter,3)],I_distr(m));
        Opt_Param_distr(1,m)=C0(I_d1);
        Opt_Param_distr(2,m)=a(I_d2);
        Opt_Param_distr(3,m)=Ka(I_d3);
        Opt_Param_distr(4,m)=RMSE_iter(I_d1,I_d2,I_d3);
    end
    clear mu Ls TDm_avg ERR_avg
    % clear RMSE_iter
end
save(matfile);
toc