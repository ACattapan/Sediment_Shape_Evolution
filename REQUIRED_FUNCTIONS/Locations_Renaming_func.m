function LocationsNew=Locations_Renaming_func(Location,LocationsRenaming)
%--------------------------------------------------------------------------
% THIS CODE COMPUTES THE RELATIVE MASS LOSS AS A FUNCTION OF TRAVEL
% DISTANCE.
%--------------------------------------------------------------------------
arguments
    Location double
    LocationsRenaming double
end
%--------------------------------------------------------------------------
% if isrow(Locations)
%     dim=1;
% elseif iscolumn(Locations)
%     dim=2;
% end
% length(LocationsRenaming);
% len=max(size(LocationsRenaming));
LocationsNew=LocationsRenaming.LocationNew(LocationsRenaming.LocationOld==Location)
