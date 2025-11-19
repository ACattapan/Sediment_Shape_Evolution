function Trav_Dist_Model_V02(T,Param,Lithology,Paths,PathNum,RastersData,RastersLocations,alpha,TD_Model_file)
tic
disp('Running Travel Distance Model')
% Current_folder_path='D:\08_PhD\04_CASE_STUDIES\Piave\DATA\SURVEYS\DATA_ANALYSIS';
% TD_Model_file='TD_MODEL_p3.mat';
% matfile=fullfile(Current_folder_path,TD_Model_file);
matfile=TD_Model_file;
%--------------------------------------------------------------------------
%                      EXTRACT PATHS FOR A GIVEN LITHOLOGY
[LocPath, OrigPath]=Path_Extract(T,Paths);
%--------------------------------------------------------------------------
%                      CHOOSE WHETHER TO FIT SHAPE PARAMETERS
Fitting_Param="Y";
%--------------------------------------------------------------------------
%                      CHOOSE WHETHER TO FIT SHAPE PARAMETERS
% Possible choices are: "RMSE" or "KS_test"
Minimization_Criteria="RMSE";
%--------------------------------------------------------------------------
%                      DEFINE RANGES FOR MODEL PARAMETERS VALUES
mu_max=1.0;
C_max=1.0;
% C0=[0.2:0.2:0.8];
% a=[1.0 4 6];
% K_exp=[3 4 6.5];
C0=[0.2:0.05:0.8];
a=[0.6:0.2:8];
K_exp=[1.0:0.25:7];
K=10.^-(K_exp);
ParamSpace=size(C0,2)*size(a,2)*size(K,2);
%--------------------------------------------------------------------------
%                      INITIALIZE VARIABLES FOR EACH LOCATION
TDe_avg=zeros(2,size(LocPath{PathNum},2));
beta1=zeros(1,size(LocPath{PathNum},2));
beta2=zeros(1,size(LocPath{PathNum},2));
h_beta=zeros(1,size(LocPath{PathNum},2));
p_beta=zeros(1,size(LocPath{PathNum},2));
TD_M_avg=zeros(1,size(LocPath{PathNum},2));
RMSE_avg=zeros(1,size(LocPath{PathNum},2));
I=zeros(1,size(LocPath{PathNum},2));
Opt_Param_avg=zeros(4,size(LocPath{PathNum},2));
if Minimization_Criteria=="RMSE"
    RMSE=zeros(1,size(LocPath{PathNum},2));
    I_distr=zeros(1,size(LocPath{PathNum},2));
    Opt_Param_distr=zeros(4,size(LocPath{PathNum},2));
elseif Minimization_Criteria=="KS_test"
    p_max=NaN(1,size(LocPath{PathNum},2));
    hip_p_max=NaN(1,size(LocPath{PathNum},2));
    I_distr=zeros(1,size(LocPath{PathNum},2));
    Opt_Param_distr=zeros(4,size(LocPath{PathNum},2));
end
%--------------------------------------------------------------------------
%                      INTERNAL VARIABLES USED FOR OPTIMIZATION WITHIN THE
%                      CYCLE FOR EACH LOCATION
TDm_avg=zeros(size(C0,2),size(a,2),size(K,2));
ERR_avg=zeros(size(C0,2),size(a,2),size(K,2));
if Minimization_Criteria=="RMSE"
    RMSE_iter=zeros(size(C0,2),size(a,2),size(K,2));
elseif Minimization_Criteria=="KS_test"
    p_values=zeros(size(C0,2),size(a,2),size(K,2));
    hip=zeros(size(C0,2),size(a,2),size(K,2));
end
%--------------------------------------------------------------------------
%                      OPTIMIZATION MODEL BEGIN
%--------------------------------------------------------------------------
for m=1:size(LocPath{PathNum},2)
    m
    %----------------------------------------------------------------------
    %           ERASE VARIABLES AT EACH CYCLE TO SAVE MEMORY
    %----------------------------------------------------------------------
    RasterData_TDtest=[];
    SampSz=[];
    mu=[];
    Ls=[];
    TDm_avg=[];
    ERR_avg=[];
    if Minimization_Criteria=="RMSE"
        RMSE_iter=[];
    elseif Minimization_Criteria=="KS_test"
        hip=[];
        p_values=[];
    end
    %----------------------------------------------------------------------
    % CONVERT TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM GENERATED) IN A
    % COLUMN VECTOR, REMOVE NaNs AND DEFINE THE SAMPLE SIZE
%     RasterData_TDtest=RastersData{1,find(RastersLocations==LocPath{PathNum}(m))}(:);
    RasterData_TDtest=RastersData{1,RastersLocations==LocPath{PathNum}(m)}(:);
    RasterData_TDtest(isnan(RasterData_TDtest))=[];
    RasterData_TDtest=sort(RasterData_TDtest);
    SampSz=size(RasterData_TDtest,1);
    %----------------------------------------------------------------------
    % COMPUTE THE AVERAGE TRAVEL DISTANCE FROM EMPIRICAL DATA (DEM
    % GENERATED) IN [m] AND ITS STANDARD DEVIATION
    %----------------------------------------------------------------------
    TDe_avg(1,m)=mean(RasterData_TDtest,'all','omitnan');
    TDe_avg(2,m)=std(RasterData_TDtest);
    %----------------------------------------------------------------------    
    if size(T.Pebble_Id(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m)),1)>0
        %------------------------------------------------------------------
        % EXTRACT THE SHAPE PARAMETERS VALUES TO BE USED BY THE MODEL
        %------------------------------------------------------------------
        if or(Fitting_Param=="Y",Fitting_Param=="y")
            [beta1(m),beta2(m),h_beta(m),p_beta(m)]=BetaFit(T,Param,alpha,Lithology,OrigPath{PathNum}(m),LocPath{PathNum}(m),"N","pdf");
            C=betarnd(beta1(m),beta2(m),SampSz,1);            
        else
            C=T.(Param)(T.Lithology==Lithology&T.Location==LocPath{PathNum}(m)&T.Origin==OrigPath{PathNum}(m));
        end
        C=sort(C);
        %------------------------------------------------------------------
        %           ACTUAL MODEL TEST FOR EACH PARAMETER
        %------------------------------------------------------------------
        Iter_N=0;
        for i=1:size(C0,2)
            for j=1:size(a,2)
                for k=1:size(K,2)
                    Iter_N=Iter_N+1
                    mu=mu_func(C,mu_max,C0(i),C_max,a(j));
                    Ls=Length_func(mu,K(k));
                    TDm_avg(i,j,k)=mean(Ls,'all');
                    ERR_avg(i,j,k)= sqrt((TDe_avg(1,m)-TDm_avg(i,j,k))^2);
                    if Minimization_Criteria=="RMSE"
                        RMSE_iter(i,j,k)=sqrt(mean((Ls-RasterData_TDtest).^2));
                    elseif Minimization_Criteria=="KS_test"
                        [hip(i,j,k),p_values(i,j,k)]=kstest2(RasterData_TDtest,Ls,'Alpha',alpha);
                    end
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
        Opt_Param_avg(3,m)=K(I3);
        Opt_Param_avg(4,m)=ERR_avg(I1,I2,I3);
        %------------------------------------------------------------------
        %           OPTIMUM PARAMETERS FOR DISTRIBUTIONS
        %------------------------------------------------------------------
        if Minimization_Criteria=="RMSE"
            [RMSE(m),I_distr(m)]=min(RMSE_iter,[],'all','linear');
            [I_d1,I_d2,I_d3]=ind2sub([size(RMSE_iter,1) size(RMSE_iter,2) size(RMSE_iter,3)],I_distr(m));
            Opt_Param_distr(1,m)=C0(I_d1);
            Opt_Param_distr(2,m)=a(I_d2);
            Opt_Param_distr(3,m)=K(I_d3);
            Opt_Param_distr(4,m)=RMSE_iter(I_d1,I_d2,I_d3);            
        elseif Minimization_Criteria=="KS_test"
            [p_max(m),I_distr(m)]=max(p_values,[],'all','linear');
            [Ip1,Ip2,Ip3]=ind2sub([size(p_values,1) size(p_values,2) size(p_values,3)],I_distr(m));
            hip_p_max(m)=hip(I1,I2,I3);
            Opt_Param_distr(1,m)=C0(Ip1);
            Opt_Param_distr(2,m)=a(Ip2);
            Opt_Param_distr(3,m)=K(Ip3);
            Opt_Param_distr(4,m)=p_values(Ip1,Ip2,Ip3);
        end
    end
    clear mu Ls TDm_avg ERR_avg
%     if Minimization_Criteria=="RMSE"
%         clear RMSE_iter
%     elseif Minimization_Criteria=="KS_test"
%         clear hip p_values
%     end
end
save(matfile);
toc