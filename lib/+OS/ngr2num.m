function [ absEasting, absNorthing ] = grid2num(varargin)
    % Converts an OS grid-letter based location reference to an OS all-numeric reference.
    %
    % Usage:
    %
    %    [ absEasting, absNorthing ] = OS.grid2num( gridCode, relativeEasting, relativeNorthing );
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
    
    if nargin == 3
        gridCode         = varargin{1};
        relativeEasting  = varargin{2}; 
        relativeNorthing = varargin{3};
    elseif nargin == 1
        parts = strsplit(varargin{1}, ' ');
        
        gridCode         = parts{1};
        relativeEasting  = str2num(parts{2}; 
        relativeNorthing = str2num(parts{3});
    end
    
    [eastGridCell, northGridCell] = find(strcmp(gridCode, OS.nationalGrid));
    
    absEasting  = (eastGridCell - 1)  * 100000 + OS.gridDistanceMetres(relativeEasting);
    absNorthing = (northGridCell - 1) * 100000 + OS.gridDistanceMetres(relativeNorthing);
end

