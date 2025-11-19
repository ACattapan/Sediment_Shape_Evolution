function X=Stats_Extract_Path(T,Parameter,alpha,Quantiles_Used,Lith,OrigPath,LocPath)
% 
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    alpha (1,1) double
    Quantiles_Used (1,:) double
    Lith (:,:) string
    OrigPath (:,:) string
    LocPath (1,:) %{mustBeNumeric(LocPath)}
end
%--------------------------------------------------------------------------
tic
disp('Running Stats Extract Path')
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

Quantiles_Parameters=zeros(size(Locations,2),size(Quantiles_Used,2));
%--------------------------------------------------------------------------
Location=[];
Distance=[]; %N.B. Distance is produced in km here!!!
TravDist=[];
SampleSize=[];
Min=[];
Max=[];
Average=[];
Median=[];
Mode=[];
StDev=[];
CV=[];
SEM=[];
cil=[];
cih=[];
for m=1:size(Locations,2)
    sz=size(T(T.Lithology==Lith&T.Location==Locations(m)&T.Origin==OrigPath(m),1),1);
    if sz>0
        Sed_Param=T.(Parameter)(T.Lithology==Lith&T.Location==Locations(m)&T.Origin==OrigPath(m));
        Location=[Location; Locations(m)];
        Distance=[Distance; LocDist(m)/1000];
        TravDist=[TravDist; (LocDist(1)-LocDist(m))/1000];
        SampleSize=[SampleSize; sz];
        Minimum=min(Sed_Param);
        Min=[Min; Minimum];
        Maximum=max(Sed_Param);
        Max=[Max; Maximum];
        xbar=mean(Sed_Param);
        Average=[Average; xbar];
        x50=quantile(Sed_Param,0.5);
        Median=[Median; x50];
        Mode=[Mode; mode(Sed_Param)];
        StDeviation=std(Sed_Param);
        StDev=[StDev; StDeviation];
        CV=[CV; StDeviation/xbar];
        SEM=[SEM; StDeviation/sqrt(sz)];
        ts = tinv([alpha/2  1-alpha/2],sz-1);      % T-Score
        cil=[cil; xbar+ts(1)*StDeviation/sqrt(sz)];
        cih=[cih; xbar+ts(2)*StDeviation/sqrt(sz)];
        Quantiles_Parameters(m,:)=quantile(Sed_Param,Quantiles_Used./100);
        %scatter(Quantiles_Parameters(m,:),Quantiles_Used./100, 'filled')
    end
    %----------------------------------------------------------------------
end
Quantiles_Names="q_"+string(Quantiles_Used);
X_Quantiles=table('Size',size(Quantiles_Parameters),'VariableTypes',repmat("double",1, size(Quantiles_Parameters,2)),'VariableNames',Quantiles_Names);
for i=1:size(Quantiles_Parameters,2)
    X_Quantiles{:,i}=Quantiles_Parameters(:,i);
end
X=table(Location,Distance,TravDist,SampleSize,Min,Max,Average,Median,Mode,StDev,CV,SEM,cil,cih);
X=[X X_Quantiles];
toc
end
