function Trav_Dist_Model_Distributions(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
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
end
%--------------------------------------------------------------------------
tic
disp('Running Travel Distance Model Distributions Local Global')
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
TD_Rasters=cell(1,size(Locations,2));
Emp_SampSz=cell(1,size(Locations,2));
% TD_E_Stats=zeros(2,size(Locations,2));
for m=1:size(Locations,2)
    TD_Rasters{m}=RastersData{1,RastersLocations==Locations(m)}(:); % Select the RasterData of Travel Distances of Location m
    TD_Rasters{m}(isnan(TD_Rasters{m}))=[];                         % Remove NaN values
    TD_Rasters{m}=sort(TD_Rasters{m});                              % Sort Travel Distances for Location m
    Emp_SampSz{m}=size(TD_Rasters{m},1);                            % Number of pixels in the RasterData of Location m
    %----------------------------------------------------------------------
    % COMPUTE THE AVERAGE AND STANDARD DEVIATION OF TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM
    % GENERATED) IN [m]
    %----------------------------------------------------------------------
    TD_E_Stats(1,m)=mean(TD_Rasters{m}(:),'all','omitnan');         % Mean of empirical Travel Distance from RasterData for location m
    TD_E_Stats(2,m)=std(TD_Rasters{m}(:));                          % St.Dev. of empirical Travel Distance from RasterData for location m
end
%--------------------------------------------------------------------------
%                      MODEL PARAMETERS VALUES
mu_max=1.0;
C_max=1.0;
%--------------------------------------------------------------------------
%                      LOCAL PARAMETERS TO BE USED ONLY WHEN RUNNING THE
%                      FUNCTION ALONE
%--------------------------------------------------------------------------
Ka=10.^-(Ka_exp);
% LocalParamSpace=size(Ka,2);
% GlobalParamSpace=size(C0,2)*size(a,2);
%--------------------------------------------------------------------------
Locations_Analysis=Locations;
Locations_Analysis(end-1:end)=[];
Origin_Analysis=OrigPath{PathNum};
Origin_Analysis(end-1:end)=[];
% TD_E_Stats(:,end-1:end)=[];
TD_Rasters_Analysis=TD_Rasters;
TD_Rasters_Analysis{end}=[];
TD_Rasters_Analysis{end-1}=[];
%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
%           INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE CYCLE
%--------------------------------------------------------------------------
%                      INITIALIZE VARIABLES
RMSE_Global=zeros(size(C0,2),size(a,2));                                % Sum of Root Mean Square Errors of all Locations Analyzed

I_Ka_Local=zeros(size(C0,2),size(a,2),size(Locations_Analysis,2));      % Index of the value of Ka that produces the minimum RMSE for a given location and for a given couple of C0 and a parameters

TD_Local_ij=cell(size(C0,2),size(a,2),size(Locations_Analysis,2));      % Modelled Travel Distances for a given location and for a given couple of C0 and a parameters

QuantilesP=[0.01:0.01:1.00];                                            % Quantiles used for the discretization of the Empirical and Modelled Distributions of travel Distances

TD_E_Quant=zeros(size(QuantilesP,2),size(Locations_Analysis,2));        % Computation of Empirical quantiles of the Travel Distances for each Locations Analyzed
for m=1:size(Locations_Analysis,2)
    TD_E_Quant(:,m)=quantile(TD_Rasters_Analysis{1,m}(:),QuantilesP);
    Param{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
    Sample_Size(m)=size(Param{m}(:),1);                                 % Size of the sample in location m
end
%--------------------------------------------------------------------------
% Cycle through all parameters that are considered common for a given
% lithology (C0 and a)
for i=1:size(C0,2)
    for j=1:size(a,2)
        % Cycle through all locations of the path considered
        for m=1:size(Locations_Analysis,2)
            %--------------------------------------------------------------
            % Delete variables for each location to avoid rewriting on them
            % and indexing errors
            clear FreqDistr C mu Ls TD_Local TD_M_Quant TD_Local_ij_Quant ERR_Local
            %--------------------------------------------------------------
            % SHAPE PARAMETER USED TO ESTIMATE TRAVEL DISTANCE
            % Param{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
            % Sample_Size=size(Param{m}(:),1); % Size of the sample in location m
            if Sample_Size(m)>0
                %------------------------------------------------------------------
                % EXTRACT THE SHAPE PARAMETERS VALUES TO BE USED BY THE MODEL
                %------------------------------------------------------------------
                if or(Fitting_Param=="Y",Fitting_Param=="y")
                    % [beta1(m),beta2(m),h_beta(m),p_beta(m)]=BetaFit(T,Parameter,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
                    % C=betarnd(beta1(m),beta2(m),Emp_SampSz{m},1);
                    FreqDistr{PathNum,1}=Freq_Distr_Extract(T,Parameter,FitTest,alpha,Lithology,OrigPath{PathNum},LocPath{PathNum},'N','pdf');
                    C=betarnd(FreqDistr{PathNum,1}{1,Locations(m)}(1,1),FreqDistr{PathNum,1}{1,Locations(m)}(1,2),Emp_SampSz{m},1);
                else
                    C=Param{m}(:);
                end
                C=sort(C);
                %----------------------------------------------------------
                % CHOOSE WHETHER C VALUES LOWER THAN C0 HAVE TO BE DELTED (Y) ON NOT (N)
                if or(CLowerC0=="Y",CLowerC0=="y")
                    C(C<C0(i))=[];
                end
                %----------------------------------------------------------
                % TRANSFORMATION MODEL
                mu=mu_func(C,mu_max,C0(i),C_max,a(j)); %C,mu_max,C0,C_max,a
                %----------------------------------------------------------           
                % Ls=cell(1,size(Ka,2));
                TD_Local=cell(1,size(Ka,2));
                TD_M_Quant=zeros(size(QuantilesP,2),size(Ka,2));
                TD_Local_ij_Quant=zeros(size(QuantilesP,2),1);
                ERR_Local=zeros(1,size(Ka,2));
                for k=1:size(Ka,2)
                    % Ls{k}=Length_func(mu,Ka(k));
                    % TD_Local{k}=Ls{k}(:);
                    TD_Local{k}=Length_func(mu,Ka(k));
                    TD_M_Quant(:,k)=quantile(TD_Local{1,k}(:),QuantilesP);
                    ERR_Local(k)=rmse(TD_E_Quant(:,m),TD_M_Quant(:,k));
                end
                %----------------------------------------------------------
                % IDENTIFICATION OF THE POSITION INDEX OF THE Ka THAT
                % MINIMIZES THE ABSOLUTE ERROR FOR A GIVEN LOCATION m AND
                % COMBINATION OF GLOBAL PARAMETERS C0 AND a
                % AND RECORD OF THE ASSOCIETED TRAVEL DISTANCE 
                [min_ERR_Avg_M,I_Ka_Local(i,j,m)]=min((ERR_Local),[],'all','linear');
                TD_Local_ij{i,j,m}=TD_Local{I_Ka_Local(i,j,m)};
                % Convergence_Figure_check(i,j,k,m,QuantilesP,TD_E_Quant,TD_M_Quant,I_Ka_Local);
                TD_Local_ij_Quant(:,1)=quantile(TD_Local_ij{i,j,m}(:),QuantilesP);
                RMSE_Global(i,j)=RMSE_Global(i,j)+rmse(TD_Local_ij_Quant,TD_E_Quant(:,m));
            end
        end
        %------------------------------------------------------------------
    end
end
%--------------------------------------------------------------------------
%           OPTIMUM PARAMETERS FOR AVERAGES
%--------------------------------------------------------------------------
[min_RMSE_Global,I_Global]=min(RMSE_Global,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_Dist,I2_Dist]=ind2sub([size(RMSE_Global,1) size(RMSE_Global,2)],I_Global);
Opt_Param_Dist{1}=C0(I1_Dist);
Opt_Param_Dist{2}=a(I2_Dist);
for m=1:size(I_Ka_Local,3)
    I3_Dist(m)=I_Ka_Local(I1_Dist,I2_Dist,m);
end
Opt_Param_Dist{3}=Ka_exp(I3_Dist);
Opt_Param_Dist{4}=RMSE_Global(I1_Dist,I2_Dist);
% for m=1:size(TD_Local_ij,3)
%     TD_Avg(m)=TD_Local_ij(I1_Dist,I2_Dist,m);
% end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
FigNum=length(findobj('type','figure'));
Fig_Surf_RMSE_Dist_Local_Global=figure(FigNum+1)

[X_coor,Y_coor] = meshgrid(C0,a);
surf(X_coor,Y_coor,RMSE_Global');
hold on;
scatter3(C0(I1_Dist),a(I2_Dist),RMSE_Global(I1_Dist,I2_Dist),100,'r','filled');
% title('RMSE for the Distributions');
xlabel('C0 [-]');
ylabel('a [-]');
zlabel('Cumulative RMSE for all locations [m]');
% zlim([0 10000])
FigNameRMSE_Dist='RMSE_Surf_for_Distributions_of_Model_of_%s_Path_Number_%d_Local_Global.fig';
saveas(Fig_Surf_RMSE_Dist_Local_Global,sprintf(FigNameRMSE_Dist,Lithology,PathNum))

save(matfile);
toc