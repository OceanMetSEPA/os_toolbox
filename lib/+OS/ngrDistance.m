function [ metres ] = ngrDistance(string)
    % Converts an OS grid reference distance into metres.
    %
    % The distance references used in the OS 3-part grid references
    % represent different distances depending on the number of digits
    % quoted and any leading zeros. For example, '7' represents 70,000 m,
    % '707' represents 70,700 m and '007' represents 7,000 m. 
    %
    % Distances should be strings of 5 digits or less.
    %
    %
    % Usage:
    %
    %    [ metres ] = OS.ngrDistance(distance);
    %
    %
    % OUTPUT:
    %    
    %   metres:  a distance in metres
    %
    % EXAMPLES:
    %
    %    OS.gridDistanceMetres('007')
    %     ans = 
    %         7000
    %
    %
    % DEPENDENCIES:
    %
    %  - None
    %
    
    metres = [];
    leadingZeros = 0;
    
    % Count the number of leading zeros
    for i = 1:length(string)
        if string(i) == '0'
            leadingZeros =+ 1;
        else
            break
        end
    end
    
    if length(string) < 6
        % Figure out how many digits were pased in and scale up accordingly
        metres = str2double(string)*(10^(5-ceil(log10(str2double(string)))-leadingZeros));        
    else
       error('Distance passed in is greater than 5 digits.');        
    end
end

