function [ grid ] = nationalGrid(varargin)
    % Returns a cell array representing the Ordnance Survey Nation Grid labelling system.
    % This funciton is used to drive the other conversion functions found
    % in +OS.
    %
    % Note, the grid is indexed from the btoom left (i.e. southwest) corner
    % and so appears upside down when viewing the cell array itself.
    %
    %
    % Usage:
    %
    %    [ gr ] = OS.nationalGrid();
    %
    %
    % OUTPUT:
    %    
    %   gr:  cell array describing OS grid labels
    %
    %EXAMPLES:
    %
    %    OS.nationalGrid
    %     ans = 
    %         'SV'    'SW'    'SX'    'SY'    'SZ'    'TV'      []
    %           []    'SR'    'SS'    'ST'    'SU'    'TQ'    'TR'
    %         ...
    %           []      []      []      []    'HP'      []      []
    %
    %
    % DEPENDENCIES:
    %
    %  - None
    %
    
    firstLetters  = {'S' 'T' ; 'N' 'O' ; 'H' 'J'};
    secondLetters = flipud(reshape(strrep(char(65:90),'I',''),5,5)');

    grid = cell(15, 10);

    for r = 1:15
        for c = 1:10
            firstLetter  = firstLetters{ceil(r/5), ceil(c/5)};
            secondLetter = secondLetters(r-(ceil(r/5)-1)*5, c-(ceil(c/5)-1)*5);
            grid{r,c} = [firstLetter secondLetter];
        end
    end
                
end

