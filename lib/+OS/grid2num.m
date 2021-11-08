function [ absEasting, absNorthing ] = grid2num( gridCode, easting, northing )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   grid2num.m  $
% $Revision:   1.3  $
% $Author:   andrew.berkeley  $
% $Date:   Sep 03 2014 11:46:52  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Converts an OS grid-letter based location reference to an OS all-numeric reference.
    %
    % Usage:
    %
    %    [ absEasting, absNorthing ] = OS.grid2num( gridCode, easting, northing );
    %
    % where gridCode, easting and northing represent a 3-part grid-letter
    % based OS reference. Note, all arguments must be strings. The reason
    % for this is that the easting and northings may contain leading zeros
    % which are sigificant to the meaning of the reference.
    %
    % OUTPUT:
    %    
    %   absEasting: the easting, in m from the absolute OS origin, of the
    %               point passed in
    %   absNorthing: the northing, in m from the absolute OS origin, of the
    %                point passed in
    %
    %EXAMPLES:
    %
    %     [east,north] = OS.grid2num('HU','123','456')
    %     east =
    %           412300
    %     north =
    %          1145600
    %
    % DEPENDENCIES:
    %
    %  - +OS/code2indx.m
    %  - +OS/gridDistanceMetres.m
    %
    
    [eastGridCell, northGridCell] = OS.code2indx(gridCode);
    
    absEasting  = (eastGridCell - 1)  * 100000 + OS.gridDistanceMetres(easting);
    absNorthing = (northGridCell - 1) * 100000 + OS.gridDistanceMetres(northing);
end

