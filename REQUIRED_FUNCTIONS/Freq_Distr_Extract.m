function FreqDistr=Freq_Distr_Extract(T,Parameter,FitTest,alpha,Lithology,OrigPath,LocPath,Plotting_flag,Plot_Type)
% FreqDistr=Freq_Distr_Extract(T,Parameter,alpha,Lithology,OrigPath{PathNum},LocPath{PathNum},"Y","cdf")
%--------------------------------------------------------------------------
% FOR EACH LOCATION BELONGING TO A PATH, FIT THE VALUES OF THE CHOSEN 
% PARAMETER WITH A PROBABILITY DISTRIBUTION.
% THE MODEL BY DEFAULT USES A BETA BISTRIBUTION FOR PARAMETERS WHOSE VALUES
% ARE BOUNDED UPWARDS BY 1 (e.g. Circularity or elongation) 
% OR A KERNEL SMOOTHING DISTRIBUTION FOR PARAMETERS WITH MAXIMUM VALUE 
% HIGHER THAN 1 (e.g. b or rP). 
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    FitTest (1,1) string
    alpha (1,1) double
    Lithology (:,:) string
    OrigPath (1,:) string
    LocPath (1,:) %{mustBeNumeric(LocPath)}
    Plotting_flag (1,1) string
    Plot_Type (1,1) string
end
%--------------------------------------------------------------------------
tic
disp('Running Frequency Distribution Extract')
%-----CHECK WHETHER THE LOCATIONS FOR THE SPECIFIC PATH WERE INTRODUCED----
if isempty(LocPath)
    Locations=unique(T.Location);
    LocDist=sort(unique(T.Distance),'descend');
else
    Locations=LocPath;
    for i=1:size(Locations,2)
        LocDist(i)=unique(T.Distance(T.Location==Locations(1,i)));
    end
    LocDist=sort(LocDist,'descend');
end
%--------------------------------------------------------------------------
%       CHOOSE WHICH METHOD TO USE FOR FITTING PROB. DISTRIBUTION
%--------------------------------------------------------------------------
% IF THE MAXIMUM VALUE OF THE PARAMETER IS LESS THEN 1, 
% USE A BETA DISTRIBUTION TO FIT HTE DATA
% IF THE MAXIMUM VALUE OF THE PARAMETER IS MORE THEN 1, 
% USE A KERNEL SMOOTHING DISTRIBUTION TO FIT HTE DATA
MaxDataLoc=[];
for k=1:size(Locations,2)
    sz=size(T.(Parameter)(T.Lithology==Lithology&T.Location==Locations(k)&T.Origin==OrigPath(k),1),1);
    if sz>0
        MaxDataLoc(k)=max(T.(Parameter)(T.Lithology==Lithology&T.Origin==OrigPath(k)&T.Location==Locations(k)));
    end
end
MaxData=max(MaxDataLoc);
%--------------------------------------------------------------------------
%                   USE A BETA DISTRIBUTION TO FIT HTE DATA
%--------------------------------------------------------------------------
if MaxData<=1
    for k=1:size(Locations,2)
        sz=size(T.(Parameter)(T.Lithology==Lithology&T.Location==Locations(k)&T.Origin==OrigPath(k),1),1);
        if sz>0
            [a,b,h,p]=BetaFit(T,Parameter,FitTest,alpha,Lithology,OrigPath(k),Locations(k),Plotting_flag,Plot_Type);
            FreqDistr{1,k}(1,:)=[a,b];
            FreqDistr{2,k}(1,:)=[h,p];
        else
            FreqDistr{1,k}(1,:)=[NaN,NaN];
            FreqDistr{2,k}(1,:)=[NaN,NaN];
        end
    end
end
%--------------------------------------------------------------------------
%             USE A KERNEL SMOOTHING DISTRIBUTION TO FIT HTE DATA
%--------------------------------------------------------------------------
if MaxData>1
    for k=1:size(Locations,2)
        sz=size(T.(Parameter)(T.Lithology==Lithology&T.Location==Locations(k)&T.Origin==OrigPath(k),1),1);
        if sz>0
            [a,b,h,p]=KSdensityFit(T,Parameter,FitTest,alpha,Lithology,OrigPath(k),Locations(k),Plotting_flag,Plot_Type);
            FreqDistr{1,k}(:,:)=[a,b];
            FreqDistr{2,k}(1,:)=[h,p];
        else
            FreqDistr{1,k}(:,:)=[NaN,NaN];
            FreqDistr{2,k}(1,:)=[NaN,NaN];
        end
    end
end
%--------------------------------------------------------------------------
% SampleSize=[];
% FreqDistr=[];
% for i=1:size(Origins,1)
%     for j=1:size(Lithologies,1)
%         for k=1:size(Locations,1)
%             SampleSize(i,j,k)=size(T.Pebble_Id(T.Lithology==Lithologies(j)&T.Location==Locations(k)&T.Origin==Origins(i)),1);
%             %----------------------------------------------------------------------
%             % COMPUTE THE GRAIN SIZE DISTRIBUTION FOR EACH LOCATION AND FOR EACH
%             % LITHOLOGY AND ORIGIN
%             if SampleSize(i,j,k)>0
%                 [f,x]=ecdf(T.b(T.Lithology==Lithologies(j)&T.Location==Locations(k)&T.Origin==Origins(i)));
%                 GSD=cat(2,f,x);
%                 FreqDistr{i,j,k}=GSD;
%                 clear f x GSD
%
%             end
%             %----------------------------------------------------------------------
%         end
%     end
% end
toc