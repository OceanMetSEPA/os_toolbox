function [ ngrCode, relativeEasting, relativeNorthing ] = num2ngr(absEasting, absNorthing)
    % Converts an OS all-numeric location reference to an OS grid-letter based location 
    % reference.
    %
    %
    % Usage:
    %
    %    [ ngrCode, relativeEasting, cellNorthing ] = OS.num2ngr(absEasting, absNorthing);
    %
    % where absEasting, absNorthing represent a 2-part all-numeric OS location reference
    %
    % OUTPUT:
    %    
    %   ngrCode    :
    %   relativeEasting : the easting, in m from the OS grid cell origin, of the
    %                 point passed in
    %   relativeNorthing: the northing, in m from the OS grid cell origin, of the
    %                 point passed in
    %
    %EXAMPLES:
    %
    %     [letter,east,north] = OS.num2ngr(412300,1145600)
    %     letter =
    %     HU
    %     east =
    %            12300
    %     north =
    %            45600
    %
    % DEPENDENCIES:
    %
    %  - +OS/nationalGrid.m
    %

    % Get reference grid
    grid = OS.nationalGrid();

    % Get indices for grid cell
    cellEasting  = ceil(absEasting/100000);
    cellNorthing = ceil(absNorthing/100000);
    
    % Set values
    ngrCode          = grid{cellNorthing, cellEasting};
    relativeEasting  = absEasting  - (cellEasting - 1) * 100000;
    relativeNorthing = absNorthing - (cellNorthing - 1) * 100000;
    
    if nargout < 3
        ngrCode = strjoin({ngrCode num2str(relativeEasting) num2str(relativeNorthing)}, ' ')
    end
end

