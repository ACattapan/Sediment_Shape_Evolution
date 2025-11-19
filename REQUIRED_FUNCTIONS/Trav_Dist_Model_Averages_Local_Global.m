function Trav_Dist_Model_Averages_Local_Global(T,Parameter,Fitting_Param,FitTest,CLowerC0,C0,a,Ka_exp,Lithology,Origins,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
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
disp('Running Travel Distance Model Averages Local Global')
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
TD_E_Stats=zeros(2,size(Locations,2));
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
% ParamSpace=size(C0,2)*size(a,2)*size(Ka,2);
%--------------------------------------------------------------------------
Locations_Analysis=Locations;
Locations_Analysis(end-1:end)=[];
Origin_Analysis=OrigPath{PathNum};
Origin_Analysis(end-1:end)=[];
TD_E_Stats(:,end-1:end)=[];
%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
%           INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE CYCLE
%--------------------------------------------------------------------------
%                      INITIALIZE VARIABLES
RMSE_Avg_M=zeros(size(C0,2),size(a,2));          % Root Mean Square Error for the Average of the Modelled Travel Distances
RMSE_M_Avg=zeros(size(C0,2),size(a,2));          % Root Mean Square Error for the Model of the Average fo Travel Distances
I_ka_Avg_M=zeros(size(C0,2),size(a,2),size(Locations_Analysis,2));
I_ka_M_Avg=zeros(size(C0,2),size(a,2),size(Locations_Analysis,2));
TD_Avg_M_ij=zeros(size(C0,2),size(a,2),size(Locations_Analysis,2));
TD_M_Avg_ij=zeros(size(C0,2),size(a,2),size(Locations_Analysis,2));
%--------------------------------------------------------------------------
for i=1:size(C0,2)
    %----------------------------------------------------------------------
    for j=1:size(a,2)
        % Cycle through all locations of the path considered
        for m=1:size(Locations_Analysis,2)
            %--------------------------------------------------------------
            % Delete variables for each location to avoid rewriting on them
            % and indexing errors
            clear C C_Avg mu mu_Avg Ls Ls_Avg TD_Avg_M TD_M_Avg ERR_Avg_M ERR_M_Avg
            %--------------------------------------------------------------
            % SHAPE PARAMETER USED TO ESTIMATE TRAVEL DISTANCE
            Param=T.(Parameter)(T.Lithology==Lithology&T.Location==Locations_Analysis(m)&T.Origin==Origin_Analysis(m)); % Parameter values measured in location m
            Sample_Size=size(Param,1); % Size of the sample in location m
            if Sample_Size>0
                C=sort(Param);
                %----------------------------------------------------------
                % CHOOSE WHETHER C VALUES LOWER THAN C0 HAVE TO BE DELTED (Y) ON NOT (N)
                if or(CLowerC0=="Y",CLowerC0=="y")
                    C(C<C0(i))=[];
                end
                %----------------------------------------------------------
                C_Avg=mean(C);
                %----------------------------------------------------------
                mu=mu_func(C,mu_max,C0(i),C_max,a(j));
                mu_Avg=mu_func(C_Avg,mu_max,C0(i),C_max,a(j));
                %----------------------------------------------------------
                for k=1:size(Ka,2)
                    Ls{k}=Length_func(mu,Ka(k));
                    Ls_Avg{k}=Length_func(mu_Avg,Ka(k));
                    TD_Avg_M(k)=mean(Ls{k},'all');
                    TD_M_Avg(k)=Ls_Avg{k};
                    ERR_Avg_M(k)= (TD_E_Stats(1,m)-TD_Avg_M(k)).^2;
                    ERR_M_Avg(k)= (TD_E_Stats(1,m)-TD_M_Avg(k)).^2;
                    %------------------------------------------------------
                end
                %----------------------------------------------------------
                % IDENTIFICATION OF THE POSITION INDEX OF THE Ka THAT
                % MINIMIZES THE ABSOLUTE ERROR FOR A GIVEN LOCATION m AND
                % COMBINATION OF GLOBAL PARAMETERS C0 AND a
                % AND RECORD OF THE ASSOCIETED TRAVEL DISTANCE 
                [min_ERR_Avg_M,I_ka_Avg_M(i,j,m)]=min((ERR_Avg_M),[],'all','linear');
                [min_ERR_M_Avg,I_ka_M_Avg(i,j,m)]=min((ERR_M_Avg),[],'all','linear');
                TD_Avg_M_ij(i,j,m)=TD_Avg_M(I_ka_Avg_M(i,j,m));
                TD_M_Avg_ij(i,j,m)= TD_M_Avg(I_ka_M_Avg(i,j,m));
            end
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % DEFINE THE SET OF BEST SET OF Kas FOR THE WHOLE CATCHMENT AND FOR
        % COMBINATION OF GLOBAL PARAMETERS C0 AND a
        for m=1:size(Locations_Analysis,2)
            L_Avg_M(m)=TD_Avg_M_ij(i,j,m);
            L_M_Avg(m)=TD_M_Avg_ij(i,j,m);
        end
        RMSE_Avg_M(i,j)= rmse(L_Avg_M(1:m),TD_E_Stats(1,1:m));
        RMSE_M_Avg(i,j)= rmse(L_M_Avg(1:m),TD_E_Stats(1,1:m));
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%           OPTIMUM PARAMETERS FOR AVERAGES
%--------------------------------------------------------------------------
[min_RMSE_Avg_M,I_Avg_M]=min(RMSE_Avg_M,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_AM,I2_AM]=ind2sub([size(RMSE_Avg_M,1) size(RMSE_Avg_M,2)],I_Avg_M);
Opt_Param_Avg{1}=C0(I1_AM);
Opt_Param_Avg{2}=a(I2_AM);
for m=1:size(I_ka_Avg_M,3)
    I3_AM(m)=I_ka_Avg_M(I1_AM,I2_AM,m);
end
Opt_Param_Avg{3}=Ka_exp(I3_AM);
Opt_Param_Avg{4}=RMSE_Avg_M(I1_AM,I2_AM);
for m=1:size(TD_Avg_M_ij,3)
    TD_Avg(m)=TD_Avg_M_ij(I1_AM,I2_AM,m);
end
%--------------------------------------------------------------------------
[min_RMSE_M_Avg,I_M_Avg]=min(RMSE_M_Avg,[],'all','linear');
% Convert linear indices I(m) into row,col indices basedon the
% dimensions of the matrix R2
[I1_MA,I2_MA]=ind2sub([size(RMSE_M_Avg,1) size(RMSE_M_Avg,2)],I_M_Avg);
Opt_Param_M{1}=C0(I1_MA);
Opt_Param_M{2}=a(I2_MA);
for m=1:size(I_ka_M_Avg,3)
    I3_MA(m)=I_ka_M_Avg(I1_MA,I2_MA,m);
end
Opt_Param_M{3}=Ka_exp(I3_MA);
Opt_Param_M{4}=RMSE_M_Avg(I1_MA,I2_MA);
for m=1:size(TD_M_Avg_ij,3)
    TD_M(m)=TD_M_Avg_ij(I1_MA,I2_MA,m);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% FigNum=get(gcf,'Number');
% Fig_1=figure (FigNum+1)
FigNum=length(findobj('type','figure'));
Fig_Avg_Local_Global=figure (FigNum+1)
scatter(TD_E_Stats(1,:)./1000,TD_M./1000,100,'filled')
hold on;
scatter(TD_E_Stats(1,:)./1000,TD_Avg./1000,100,'+')
axeslim=max([xlim;ylim],[],'all');
xlim([0 axeslim]);
ylim([0 axeslim]);
plot(xlim, ylim, '--k')
legend('Model of the TD from the Average Circularity','Average of Modelled TD from Circularity', 'Location','southeast')
legend('boxoff')
title('Model results Hyp:'); 
subtitle('Lithology parameters: C0, a; Environmental parameter: ka');
xlabel('Empirical Travel Distance from DTM [km]');
ylabel('Modelled Travel Distance [km]');
FigNameformat='Scatter_plot_for_Averages_of_%s_Path_Number_%d_Local_Global.fig';
saveas(Fig_Avg_Local_Global,sprintf(FigNameformat,Lithology,PathNum))
%--------------------------------------------------------------------------

if size(C0,2)>1
    FigNum=get(gcf,'Number');
    Fig_Surf_RMSE_Avg_M=figure(FigNum+1)
    surf(RMSE_Avg_M)
    hold on;
    scatter3(I2_AM,I1_AM,RMSE_Avg_M(I1_AM,I2_AM),100,'r','filled')
    title('RMSE for the Average of Model results');
    ylabel('C0');
    xlabel('a');
    % zlim([0 10000])
    FigNameRMSE_Avg='RMSE_Surf_for_Averages_of_Model_of_%s_Path_Number_%d_Local_Global.fig';
    saveas(Fig_Surf_RMSE_Avg_M,sprintf(FigNameRMSE_Avg,Lithology,PathNum))

    %----------------------------------------------------------------------

    FigNum=get(gcf,'Number');
    Fig_Surf_RMSE_M_Avg=figure(FigNum+1)
    surf(RMSE_M_Avg)
    hold on;
    scatter3(I2_MA,I1_MA,RMSE_M_Avg(I1_MA,I2_MA),100,'r','filled')
    title('RMSE for Model results of the Average Circularity');
    ylabel('C0');
    xlabel('a');
    % zlim([0 10000])
    FigNameRMSE_M='RMSE_Surf_for_Model_of_Averages_of_%s_Path_Number_%d_Local_Global.fig';
    saveas(Fig_Surf_RMSE_M_Avg,sprintf(FigNameRMSE_M,Lithology,PathNum))
end
%--------------------------------------------------------------------------

save(matfile);
toc