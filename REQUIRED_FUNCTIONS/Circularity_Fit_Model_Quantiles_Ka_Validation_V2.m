tic
disp('Running Circularity Validation model')
Transf_Mode="AVG";                      % Possible values: "AVG", "QUANTILES"
Sampl_Num=5;
Quantiles_Edges= conv([Quantiles_Used./100], [0.5, 0.5], 'valid');
Quantiles_Edges=[0 Quantiles_Edges 1];
%--------------------------------------------------------------------------
% if or(Distance_Units=="km",Distance_Units=="Km")
    % RastersData=RastersData';
    % for i=1:numel(RastersData)
    %     RastersData{i,1}=RastersData{i,1}./1000;
    % end
% end
%--------------------------------------------------------------------------
% CYCLE THROUGH ALL LOCATIONS THAT COMPOSE A PATH
for i=1:size(LocPath{PathNum,1},2)
    % READ FROM THE RASTERS DATA, THE DISTANCES OF SUCH LOCATION
    Rasters_Data_path{i,1}=RastersData{find(RastersLocations==LocPath{PathNum,1}(i),1)};
    %----------------------------------------------------------------------
    % if or(Distance_Units=="km",Distance_Units=="Km")
    %     Rasters_Data_path{i}=Rasters_Data_path{i}./1000;
    % end       
    %----------------------------------------------------------------------
    % CREATE A VECTOR OF DISTANCES OF THE SAME SIZE AS THE NUMBER OF
    % PIXELS IN THE RASTER, TIMES THE Sampl_Num: 
    % THE NUMBER OF PEBBLES MODELLED FOR EACH PIXEL
    L_model{i,1}=repmat(Rasters_Data_path{i},1,Sampl_Num);
    Model_Size=numel(L_model{i,1});
    %----------------------------------------------------------------------
    % EXTRACT A SET OF INITIAL CIRCULARITY VALUES FROM THE DISTRIBUTION
    % AT THE SOURCE: THE FIRST LOCATION OF A PATH = {1,1}
    C0_M_x{i,1}=betarnd(FreqDistr{PathNum,1}{1,1}(1),FreqDistr{PathNum,1}{1,1}(2),[Model_Size,1]);
    %----------------------------------------------------------------------
    % histogram(C0_model{i,1},'normalization','pdf')
    % hold on
    % plot(0.5:0.01:1,betapdf(0.5:0.01:1,FreqDistr{PathNum,1}{1,1}(1),FreqDistr{PathNum,1}{1,1}(2)));
    %----------------------------------------------------------------------
    if Transf_Mode=="AVG"
        dummyC=[];
        for j=1:numel(C0_M_x{i})       %For each value of initial circularity randomply sampled
            C0_model=C0_M_x{i}(j);
            a_model=Opt_Param_C_L{1,2};
            Ka_model=Opt_Param_C_L{1,3};
            dummyC(j,1)=C_of_L_C0_func(L_model{i}(1,j),mu_max,C0_model,C_max,a_model,Ka_model);
        end
    elseif Transf_Mode=="QUANTILES"
        % FOR EACH INITIAL SHAPE VALUE, COMPUTE ITS COMULATIVE RELATIVE
        % FREQUENCY
        for j=1:numel(C0_M_x{i})
            C0_M_x_quantile{i}(j,1)=find(sort(C0_M_x{i})==C0_M_x{i}(j))/numel(C0_M_x{i});
        end
        % ASSIGN EACH SHAPE VALUE TO ITS QUANTILE USING THE EDGES OF
        % QUANTILES BINS PREVIOUSLY DEFINED
        C0_M_x_quantile_indx{i}=discretize(C0_M_x_quantile{i},Quantiles_Edges);
        % USE THE SELECTED TRANSFORMATION MODE TO ESTIMATE THE CIRCULARITY
        % THAT THE PEBBLE WILL HAVE AFTER HAVING TRAVELLED A DISTANCE
        % L_model
        dummyC=[];
        for j=1:numel(C0_M_x{i})       %For each value of initial circularity randomply sampled
            C0_model=C0_M_x{i}(j);
            % quantile_indx=C0_M_x_quantile_indx(k);
            a_model=Opt_Param_C_L_quant{C0_M_x_quantile_indx{i}(j),2};
            Ka_model=Opt_Param_C_L_quant{C0_M_x_quantile_indx{i}(j),3};
            dummyC(j,1)=C_of_L_C0_func(L_model{i}(1,j),mu_max,C0_model,C_max,a_model,Ka_model);
        end
    end
    C_M_x{i,1}=dummyC;
    C_M_x_AVG(i,1)=mean(C_M_x{i});
    %----------------------------------------------------------------------
    % IF SUCH LOCATION WAS SAMPLED AND THEREFORE IS PRESENT IN X, THEN
    % CONSIDER IT
    % if ismember(LocPath{PathNum,1}(i),X{PathNum,1}.Location)
    % m IS THE INDEX OF THE ROW OF SUCH LOCATION, INSIDE X
    % m=find(X{PathNum,1}.Location==LocPath{PathNum,1}(i));%the index m runs through the locations present in X
    % end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
toc
for i=1:NumPaths
    for j=1:numel(LocPath{i})
        LocDist_path{i,1}(j,1)=(LocDist(LocPath{i}(1),1)-LocDist(LocPath{i}(j),1))./1000;
    end
end

FigNum=length(findobj('type','figure'));
Fig_Validation=figure(FigNum+1);
hold on;
scatter(X{1,1}.SourceDist,X{1,1}.Average)
scatter(LocDist_path{PathNum},C_M_x_AVG);
%--------------------------------------------------------------------------
% save("D:\08_PhD\04_CASE_STUDIES\Piave\DATA\SURVEYS\DATA_ANALYSIS\RESULTS\Validation_Arenites_AVG_NormCirc_SourceDist_Path_1_V02")
%--------------------------------------------------------------------------