tic
disp('Running Circularity Validation model')
Sampl_Num=5;
Quantiles_Edges= conv([Quantiles_Used./100], [0.5, 0.5], 'valid');
Quantiles_Edges=[0 Quantiles_Edges 1];
% CYCLE THROUGH ALL LOCATIONS THAT COMPOSE A PATH
for i=1:size(LocPath{PathNum,1},2)
    % IF SUCH LOCATION WAS SAMPLED AND THEREFORE IS PRESENT IN X, THEN
    % CONSIDER IT
    % if ismember(LocPath{PathNum,1}(i),X{PathNum,1}.Location)
        % m IS THE INDEX OF THE ROW OF SUCH LOCATION, INSIDE X
        m=find(X{PathNum,1}.Location==LocPath{PathNum,1}(i));%the index m runs through the locations present in X
        % L_model=X{PathNum,1}.TravDist(m);
        % READ FROM THE RASTERS DATA, THE DISTANCES OF SUCH LOCATION
        Rasters_Data_path{m}=RastersData{1,find(RastersLocations==X{PathNum,1}.Location(m),1)};
        %------------------------------------------------------------------
        if or(Distance_Units=="km",Distance_Units=="Km")
            Rasters_Data_path{m}=Rasters_Data_path{m}./1000;
        end
        %------------------------------------------------------------------
        % CREATE A VECTOR OF DISTANCES OF THE SAME SIZE OF THE AVAILABLE
        % PIXELS, TIMES THE Sampl_Num: THE NUMBER OF PEBBLES MODELLED FOR EACH PIXEL 
        L_model{m}=repmat(Rasters_Data_path{m},1,Sampl_Num);
        dummyC0=[];
        for j=1:size(Rasters_Data_path{m},2)
            % betarnd(FreqDistr_Norm_Circ{PathNum,1}{1,i}(1),FreqDistr_Norm_Circ{PathNum,1}{1,i}(2),[Sampl_Num,1]);
            % betarnd(FreqDistr_Norm_Circ{PathNum,1}{1,i}(1),FreqDistr_Norm_Circ{PathNum,1}{1,i}(2),[Sampl_Num,1]);
            % EXTRACT A SET OF SIZE Sampl_Num OF INITIAL VALUES OF SHAPE PARAMETER FROM THE DISTRIBUTION
            % AT THE SOURCE: THE FIRST LOCATION IN THE PATH: {1,1}
            dummyC0=[dummyC0; betarnd(FreqDistr{PathNum,1}{1,1}(1),FreqDistr{PathNum,1}{1,1}(2),[Sampl_Num,1])];
            %--------------------------------------------------------------
            % histogram(dummyC0,'normalization','pdf')
            % hold on
            % plot(0.5:0.01:1,betapdf(0.5:0.01:1,FreqDistr{PathNum,1}{1,1}(1),FreqDistr{PathNum,1}{1,1}(2)));
            %--------------------------------------------------------------
        end
        C0_M_x{m}=dummyC0;
        % FOR EACH INITIAL SHAPE VALUE, COMPUTE ITS COMULATIVE RELATIVE
        % FREQUENCY 
        for j=1:size(C0_M_x{m},1)
            C0_M_x_quantile{m}(j,1)=find(sort(C0_M_x{m})==C0_M_x{m}(j))/size(C0_M_x{m},1);
        end
        % ASSIGN EACH SHAPE VALUE TO ITS QUANTILE USING THE EDGES OF
        % QUANTILES BINS PREVIOUSLY DEFINED
        C0_M_x_quantile_indx{m}=discretize(C0_M_x_quantile{m},Quantiles_Edges);
        
        dummyC=[];
        for j=1:size(C0_M_x{m},1)       %For each value of initial circularity randomply sampled
            C0_model=C0_M_x{m}(j);
            % quantile_indx=C0_M_x_quantile_indx(k);
            a_model=Opt_Param_C_L_quant{C0_M_x_quantile_indx{m}(j),2};
            Ka_model=Opt_Param_C_L_quant{C0_M_x_quantile_indx{m}(j),3};
            dummyC(j,1)=C_of_L_C0_func(L_model{m}(1,j),mu_max,C0_model,C_max,a_model,Ka_model);
        end
        C_M_x{m}=dummyC;
        C_M_x_AVG(m)=mean(C_M_x{m});
    % end
end


