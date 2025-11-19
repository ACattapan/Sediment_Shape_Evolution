function [LocPath, OrigPath, NumPaths]=Path_Extract(Origins,Paths)
%--------------------------------------------------------------------------
% EXTRACT THE LOCATIONS THAT COMPOSE A CERTAIN PATH AND THEIR ASSOCIATED
% ORIGINS FOR EACH LOCATION
% Possible Origines are    :    CREEK=1,    RIVER=2,        SOURCE=3;
%--------------------------------------------------------------------------
arguments
    % T (:,:) table
    Origins (:,:) string
    Paths (:,:) double {mustBeNumeric(Paths)}
end
%--------------------------------------------------------------------------
tic
disp('Running Path Extract')
% Origins=string(unique(T.Origin(:)));
k=0;
for i=1:size(Paths,1)                                   % For each row of the matrix Paths
    if Paths(i,i)~=0 %~isempty(find(Paths(i,:)~=0,1))   % If the diagonal element is not empty
        k=k+1;                                          % Increase the index
        LocPath{k,1}=find(Paths(i,:)~=0);               % Insert in positionk, the list of 
        OrigPath{k,1}=Origins(Paths(i,LocPath{k,1}),1)';
    end
end
NumPaths=k;
toc

