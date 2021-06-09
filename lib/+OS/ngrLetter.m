function [ varargout ] = ngrLetter( ngrString,varargin )
% This function returns the easting & northing numbers for 2 character NGR
% codes.
%
% It can also convert strings with NGR codes to Eastings, Northings
% *************************************************************************
%
% The OS map is divided into large squares, each of which is subdivided into 25
% (5x5) smaller squares
%
% The large square determines the first letter of the 2 character code
% Letters 'H','N' and 'S' correspond to the west of the country;
% 'T' corresponds to the east of England.
% (Actually there's a small point OV in England but we'll ignore that for
% now)
%
% The position of the small square within the big square determines the
% second character.
% Each of these small squares is identified by a 2 letter code.
% The 1st character is the same of each
% The large squares are identified by the
% NGR codes consist of 2 letters
%
% INPUT:
% 2 character char, or space-delimited char starting with NGR code
%
% OUTPUT: Depends on number of outputs -
%  Single output - either
%                       1) vector [X,Y]
%                       2) struct array (with fields 'String','X','Y')
% Two outputs - [X,Y] as separate variables
%
% OPTIONS:
% sepChar (' ') - gap between letters and numbers in NGR string
% autoScale (true) - if true, it'll use the size of the NGR string
% to scale the numeric part.
% numeric (true) - return vector rather than struct
%
% EXAMPLE:
% s=ngrLetter('HU') % return [4,11]
% s=ngrLetter('HU','numeric',0) % return struct('String','HU','X',4,'Y',11)
% [easting,northing]=ngrLetter('NM 74890 95390') % return two double values
% ngrLetter('NM 749 954') % returns 174800,795300
% ngrLetter('NM 749 954','autoScale',0) % returns [100749,700954]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   ngrLetter.m  $
% $Revision:   1.1  $
% $Author:   ted.schlicke  $
% $Date:   May 28 2014 13:04:06  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=struct;
options.sepChar=' ';
options.verbose=false;
options.autoScale=true; % uses number of digits in string to estimate order of magnitude
options.numeric=true;
options=checkArguments(options,varargin);

ngrLetter1={'H','N','S','T'}; % UK 2 letter grid references all start with one of these letters (apart from tiny 'O' section in England which we'll ignore!)
ngrE={0,0,0,5}; % H,N,S all correspond to left hand side of NGR map; % T is 5
ngrN=num2cell([3,2,1,1]*5-1);% Number of cells above bottom left of grid (SV) of top left subsquare (second character = 'A')
% Subsquares are labelled by letters A:Z except for letter 'I'.
ngrLookupStruct=struct('Letter',ngrLetter1,'E',ngrE,'N',ngrN);

if ischar(ngrString)
    ngrString=cellstr(ngrString); % convert to char
end

Ns=length(ngrString);
op=cell(Ns,1);
for i=1:Ns
    s=struct('String',ngrString{i});
    ngrStringUpperCase=upper(ngrString{i});
    subString=regexp(ngrStringUpperCase,options.sepChar,'split'); % Individual strings, separated by spaces in input
    Nsi=length(subString);
    if options.verbose
        fprintf('%d of %d and string = ''%s''; number of substrings = %d\n',i,Ns,ngrString{i},Nsi)
    end
    if Nsi==2 % 3 strings? Assume it's a string with easting northing
        x=str2double(subString{1});
        y=str2double(subString{2});
        s.X=x;
        s.Y=y;
    else % Assume there's a NGR code
        ngrStringUpperCase=subString{1};
        % Find numeric values of letters:
        c1=ngrStringUpperCase(1); % 1st character
        c2=ngrStringUpperCase(2); % 2nd character
        ucl=65:90; % ascii value of upper case letters
        n2=find(c2==char(ucl)); % how far through the alphabet second character is (A=1; Z=26)
        
        % These are the top left corners of the main squares, H,N,S,T:
        luti=stringFinder(ngrLetter1,c1,'output','index');
        if isempty(luti)
            warning('OH:DEAR','Character 1 ''%c'' doesn''t appear to be a valid NGR code',c1)
            Noutput=nargout;
            op=cell(Noutput,1);
            for j=1:Noutput
                op{j}=NaN;
            end
            varargout=op;
            return
        end
        lut=ngrLookupStruct(luti);
        if strcmp(c2,'I')
            error('Second character can''t be ''I''')
        end
        %find index of subsquare:
        n=n2-1-(c2>'H'); % need to compensate for missing 'I'
        % Convert subsquare index to x,y coordinates:
        x=lut.E+mod(n,5);
        y=lut.N-floor(n/5);
        s.X=x;
        s.Y=y;
    end
    
    if Nsi==3 % Assume strings 2 & 3 define easting / northing
        es=str2double(subString{2});
        ns=str2double(subString{3});
        % Determine order of magnitude based on string length:
        eOOM=5-length(subString{2});
        nOOM=5-length(subString{3});
        if options.autoScale
            ase=power(10,eOOM);
            asn=power(10,nOOM);
            if options.verbose
                fprintf('AUTOSCALE %d ; %d\n',ase,asn);
            end
            es=es*ase;
            ns=ns*asn;
        end
        es=x*1e5+es;
        ns=y*1e5+ns;
        s.('X')=es;
        s.('Y')=ns;
    end
    op{i}=s;
end
try
    op=vertcat(op{:}); % Try to convert cell array to struct array
catch
end
if nargout<2
    if options.numeric
        x=[op.X];
        y=[op.Y];
        varargout{1}=[x(:),y(:)];
    else
        varargout{1}=op;
    end
    return
elseif nargout==2
    varargout{1}=op.X;
    varargout{2}=op.Y;
end

end

