function [ gridCode, cellEasting, cellNorthing ] = num2grid(absEasting, absNorthing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   num2grid.m  $
% $Revision:   1.2  $
% $Author:   ted.schlicke  $
% $Date:   May 28 2014 13:04:08  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Converts an OS all-numeric location reference to an OS grid-letter based location 
    % reference.
    %
    %
    % Usage:
    %
    %    [ gridCode, cellEasting, cellNorthing ] = OS.num2grid(absEasting, absNorthing);
    %
    % where absEasting, absNorthing represent a 2-part all-numeric OS location reference
    %
    % OUTPUT:
    %    
    %   gridCode    :
    %   cellEasting : the easting, in m from the OS grid cell origin, of the
    %                 point passed in
    %   cellNorthing: the northing, in m from the OS grid cell origin, of the
    %                 point passed in
    %
    %EXAMPLES:
    %
    %     [letter,east,north] = OS.num2grid(412300,1145600)
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
    gridEasting  = ceil(absEasting/100000);
    gridNorthing = ceil(absNorthing/100000);
    
    % Set values
    gridCode     = grid{gridNorthing, gridEasting};
    cellEasting  = absEasting  - (gridEasting - 1) * 100000;
    cellNorthing = absNorthing - (gridNorthing - 1) * 100000;
end

