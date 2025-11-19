function Trav_Dist_Model_Averages_V1(T,Parameter,Distance_Units,Distance_Metric,...
    Fitting_Param,FitTest,CLowerC0,mu_max,C_max,C0,a,Ka,Lithology,Origins,...
    Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
% 
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    Distance_Units (1,:) string         % e.g. "km", or "m"
    Distance_Metric (1,:) string        % e.g. "SourceDist" or "TravDist"
    Fitting_Param (1,1) string
    FitTest (1,1) string
    CLowerC0 (1,1) string
    mu_max (:,:) double
    C_max (:,:) double
    C0 (:,:) double
    a (:,:) double
    Ka (:,:) double
    Lithology (:,:) string
    Origins (:,:) string
    Paths (:,:) double
    PathNum double
    RastersData cell
    RastersLocations (:,:) double
    alpha (1,1) double
    TD_Model_file (:,:) string
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
%                      INITIALIZE VARIABLES FOR EACH LOCATION
TD_E_Stats=zeros(2,size(Locations,2));      % Average and St.Dev. of empirical Travel Distances
TD_Avg=zeros(1,size(Locations,2));          % Average of the modelled Travel Distances
TD_M=zeros(1,size(Locations,2));            % Model of the Average of Travel Distances
RMSE_Avg_M=zeros(1,size(Locations,2));      % Root Mean Square Error for the Average of the Modelled Travel Distances
RMSE_M_Avg=zeros(1,size(Locations,2));      % Root Mean Square Error for the Model of the Average fo Travel Distances
I_Avg_M=zeros(1,size(Locations,2));         % Linear Index for the minimum RMSA for the Average of Modelled Travel Distances
I_M_Avg=zeros(1,size(Locations,2));         % Linear Index for the minimum RMSA for the Model of the Average Travel Distance
Opt_Param_Avg=zeros(4,size(Locations,2));   % Optimum parameters for the Average of Modelled Travel Distances (the fourth row contains the minimum RMSE)
Opt_Param_M=zeros(4,size(Locations,2));     % Optimum parameters for the Model of the Average Travel Distance (the fourth row contains the minimum RMSE)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
% Cycle through all locations of the path considered
for m=1:size(Locations,2)
    disp(['Location: ',num2str(Locations(m))])
    %----------------------------------------------------------------------
    % CONVERT TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM GENERATED) IN A
    % COLUMN VECTOR, REMOVE NaNs AND DEFINE THE SAMPLE SIZE
    RasterData_TDtest=RastersData{1,RastersLocations==Locations(m)}(:); % Select the RasterData of Travel Distances of Location m
    RasterData_TDtest(isnan(RasterData_TDtest))=[];                     % Remove NaN values
    RasterData_TDtest=sort(RasterData_TDtest);                          % Sort Travel Distances for Location m
    Emp_SampSz=size(RasterData_TDtest,1);                               % Number of pixels in the RasterData of Location m
    %----------------------------------------------------------------------
    % COMPUTE THE AVERAGE AND STANDARD DEVIATION OF TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM
    % GENERATED) IN [m] 
    %----------------------------------------------------------------------
    TD_E_Stats(1,m)=mean(RasterData_TDtest,'all','omitnan');             % Mean of empirical Travel Distance from RasterData for location m
    TD_E_Stats(2,m)=std(RasterData_TDtest);                              % St.Dev. of empirical Travel Distance from RasterData for location m
    %----------------------------------------------------------------------
    % SHAPE PARAMETER USED TO ESTIMATE TRAVEL DISTANCE
    Param=T.(Parameter)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m)); % Parameter values measured in location m
    Sample_Size=size(Param,1);                                          % Size of the sample in location m
    if Sample_Size>0
        %------------------------------------------------------------------
        %           ACTUAL MODEL TEST FOR EACH PARAMETER
        %------------------------------------------------------------------
        %                      INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE
        %                      CYCLE FOR EACH LOCATION
        TD_Avg_M=zeros(size(C0,2),size(a,2),size(Ka,2));     
        TD_M_Avg=zeros(size(C0,2),size(a,2),size(Ka,2));     
        ERR_Avg_M=zeros(size(C0,2),size(a,2),size(Ka,2));    
        ERR_M_Avg=zeros(size(C0,2),size(a,2),size(Ka,2));    
        %------------------------------------------------------------------
        Iter_N=0;
        for i=1:size(C0,2)
            for j=1:size(a,2)
                C{m}=sort(Param);
                %------------------------------------------------------------------
                % CHOOSE WHETHER C VALUES LOWER THAN C0 HAVE TO BE DELTED (Y) ON NOT (N)
                if or(CLowerC0=="Y",CLowerC0=="y")
                    C{1,m}(C{1,m}<C0(i))=[];
                end
                C_Avg{m}=mean(C{m});
                for k=1:size(Ka,2)
                    Iter_N=Iter_N+1;
                    % disp(['Iteration: ',num2str(Iter_N)])
                    mu=mu_func(C{m},mu_max,C0(i),C_max,a(j));
                    mu_Avg=mu_func(C_Avg{m},mu_max,C0(i),C_max,a(j));
                    Ls=Length_func(mu,Ka(k));
                    Ls_Avg=Length_func(mu_Avg,Ka(k));
                    TD_Avg_M(i,j,k)=mean(Ls,'all');
                    TD_M_Avg(i,j,k)=Ls_Avg;
                    ERR_Avg_M(i,j,k)= sqrt((TD_E_Stats(1,m)-TD_Avg_M(i,j,k))^2);
                    ERR_M_Avg(i,j,k)= sqrt((TD_E_Stats(1,m)-TD_M_Avg(i,j,k))^2);
                end
            end
        end
        %------------------------------------------------------------------
        %           OPTIMUM PARAMETERS FOR AVERAGES
        %------------------------------------------------------------------
        [RMSE_Avg_M(m),I_Avg_M(m)]=min(ERR_Avg_M,[],'all','linear');
        % Convert linear indices I(m) into row,col indices basedon the
        % dimensions of the matrix R2
        [I1_AM,I2_AM,I3_AM]=ind2sub([size(ERR_Avg_M,1) size(ERR_Avg_M,2) size(ERR_Avg_M,3)],I_Avg_M(m));
        TD_Avg(m)=TD_Avg_M(I1_AM,I2_AM,I3_AM);
        Opt_Param_Avg(1,m)=C0(I1_AM);
        Opt_Param_Avg(2,m)=a(I2_AM);
        Opt_Param_Avg(3,m)=Ka(I3_AM);
        Opt_Param_Avg(4,m)=ERR_Avg_M(I1_AM,I2_AM,I3_AM);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        [RMSE_M_Avg(m),I_M_Avg(m)]=min(ERR_M_Avg,[],'all','linear');
        % Convert linear indices I(m) into row,col indices basedon the
        % dimensions of the matrix R2
        [I1_MA,I2_MA,I3_MA]=ind2sub([size(ERR_M_Avg,1) size(ERR_M_Avg,2) size(ERR_M_Avg,3)],I_M_Avg(m));
        TD_M(m)=TD_M_Avg(I1_MA,I2_MA,I3_MA);
        Opt_Param_M(1,m)=C0(I1_MA);
        Opt_Param_M(2,m)=a(I2_MA);
        Opt_Param_M(3,m)=Ka(I3_MA);
        Opt_Param_M(4,m)=ERR_M_Avg(I1_MA,I2_MA,I3_MA);        
        %------------------------------------------------------------------
    end
    clear mu mu_Avg Ls Ls_Avg TD_Avg_M TD_M_Avg ERR_Avg_M ERR_M_Avg
end
%--------------------------------------------------------------------------
% FigNum=get(gcf,'Number');
FigNum=length(findobj('type','figure'));
Fig_Avg=figure (FigNum+1)
scatter(TD_E_Stats(1,:)./1000,TD_M./1000,100,'filled')
hold on;
scatter(TD_E_Stats(1,:)./1000,TD_Avg./1000,100,'+')
plot(xlim, ylim, '--k')
legend('Model of the TD from the Average Circularity','Average of Modelled TD from Circularity')
% title('Model results Hyp: C0 = 0.48');
xlabel('Empirical Travel Distance from DTM [km]');
ylabel('Modelled Travel Distance [km]');
hold off;
FigNameformat='Scatter_plot_for_Averages_of_%s_Path_Number_%d.fig';
saveas(Fig_Avg,sprintf(FigNameformat,Lithology,PathNum))
%--------------------------------------------------------------------------
save(matfile);
toc