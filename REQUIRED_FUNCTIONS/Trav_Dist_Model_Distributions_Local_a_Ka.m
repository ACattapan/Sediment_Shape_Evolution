function Trav_Dist_Model_Distributions_Local_a_Ka(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
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
disp('Running Travel Distance Model Distributions Local a Ka')
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
RMSE_Global=zeros(size(C0,2),1);                                % Sum of Root Mean Square Errors of all Locations Analyzed

I_Local=zeros(size(C0,2),size(Locations_Analysis,2));      % Index of the value of Ka that produces the minimum RMSE for a given location and for a given couple of C0 and a parameters

TD_Local_i=cell(size(C0,2),size(Locations_Analysis,2));      % Modelled Travel Distances for a given location and for a given couple of C0 and a parameters

QuantilesP=[0.01:0.01:1.00];                                            % Quantiles used for the discretization of the Empirical and Modelled Distributions of travel Distances

TD_E_Quant=zeros(size(QuantilesP,2),size(Locations_Analysis,2));        % Computation of Empirical quantiles of the Travel Distances for each Locations Analyzed
for m=1:size(Locations_Analysis,2)
    TD_E_Quant(:,m)=quantile(TD_Rasters_Analysis{1,m}(:),QuantilesP);
    Param{m}=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
    Sample_Size(m)=size(Param{m}(:),1);                                 % Size of the sample in location m
end
%--------------------------------------------------------------------------
% Cycle through all parameters that are considered common for a given
% lithology (C0)
for i=1:size(C0,2)
        % Cycle through all locations of the path considered
        for m=1:size(Locations_Analysis,2)
            %--------------------------------------------------------------
            % Delete variables for each location to avoid rewriting on them
            % and indexing errors
            clear FreqDistr C mu Ls TD_Local TD_M_Quant TD_Local_i_Quant ERR_Local
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
                TD_Local=cell(size(a,2),size(Ka,2));
                TD_M_Quant=zeros(size(QuantilesP,2),1);
                TD_Local_i_Quant=zeros(size(QuantilesP,2),1);
                ERR_Local=zeros(size(a,2),size(Ka,2));
                %----------------------------------------------------------
                % TRANSFORMATION MODEL
                for j=1:size(a,2)
                    mu=mu_func(C,mu_max,C0(i),C_max,a(j)); %C,mu_max,C0,C_max,a
                    for k=1:size(Ka,2)
                        TD_Local{j,k}=Length_func(mu,Ka(k));
                        TD_M_Quant(:,1)=quantile(TD_Local{j,k}(:),QuantilesP);
                        ERR_Local(j,k)=rmse(TD_E_Quant(:,m),TD_M_Quant(:,1));
                    end
                end
                %----------------------------------------------------------
                % IDENTIFICATION OF THE POSITION INDEX OF THE Ka THAT
                % MINIMIZES THE ABSOLUTE ERROR FOR A GIVEN LOCATION m AND
                % COMBINATION OF GLOBAL PARAMETERS C0 AND a
                % AND RECORD OF THE ASSOCIETED TRAVEL DISTANCE 
                [min_ERR_Avg_M,I_Local(i,m)]=min((ERR_Local),[],'all','linear');
                [idx_2,idx_3]=ind2sub([size(ERR_Local,1) size(ERR_Local,2)],I_Local(i,m));
                TD_Local_i{i,m}=TD_Local{idx_2,idx_3};
                TD_Local_i_Quant(:,1)=quantile(TD_Local_i{i,m}(:),QuantilesP);
                RMSE_Global(i,1)=RMSE_Global(i)+rmse(TD_E_Quant(:,m),TD_Local_i_Quant(:,1));
            end
        end
        %------------------------------------------------------------------
end
%--------------------------------------------------------------------------
%           OPTIMUM PARAMETERS FOR AVERAGES
%--------------------------------------------------------------------------
[min_RMSE_Global,I1_Dist_Loc_aKa]=min(RMSE_Global,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
Opt_Param_Dist_Loc_aKa{1}=C0(I1_Dist_Loc_aKa);
for m=1:size(I_Local,2)
    [I2_Dist_Loc_aKa(m) I3_Dist_Loc_aKa(m)]=ind2sub([size(ERR_Local,1) size(ERR_Local,2)],I_Local(I1_Dist_Loc_aKa,m));
    Opt_Param_Dist_Loc_aKa{2}(m)=a(I2_Dist_Loc_aKa(m));
    Opt_Param_Dist_Loc_aKa{3}(m)=Ka_exp(I3_Dist_Loc_aKa(m));
end
Opt_Param_Dist_Loc_aKa{4}=RMSE_Global(I1_Dist_Loc_aKa,1);%min_RMSE_Global
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% FigNum=length(findobj('type','figure'));
% Fig_Surf_RMSE_Dist_Local_aKa=figure(FigNum+1)
% 
% [X_coor,Y_coor] = meshgrid(C0,a);
% surf(X_coor,Y_coor,RMSE_Global');
% hold on;
% scatter3(C0(I1_Dist),a(I2_Dist_Loc_aKa),RMSE_Global(I1_Dist,I2_Dist_Loc_aKa),100,'r','filled');
% % title('RMSE for the Distributions');
% xlabel('C0 [-]');
% ylabel('a [-]');
% zlabel('Cumulative RMSE for all locations [m]');
% % zlim([0 10000])
% FigNameRMSE_Dist='RMSE_Surf_for_Distributions_of_Model_of_%s_Path_Number_%d_Local_Global.fig';
% saveas(Fig_Surf_RMSE_Dist_Local_aKa,sprintf(FigNameRMSE_Dist,Lithology,PathNum))

save(matfile);
toc