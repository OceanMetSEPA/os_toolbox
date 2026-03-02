function [ varargout ] = catCoordinates(x,y,varargin )
% This function will  convert and/or transform coordinates, including Eastings/Northings and Longitude/Latitude
%
% Usage:
% catCoordinates(x,y,options)%
%   where x,y are long/lats or easting/northings
%
% Function will guess whether inputs are lat/longs or eastings/northings based on
% size. By default, it will assume you want to change between them and that
% lat/longs are specified in ESTR89 (GPS) coordinate system and eastings/northingss
% are specified in OSGB36 (Airy 1830 (National) grid).
%
% All these parameters can be specified explicitly if required.
%
% Options:
%   'from'  : either 'LL' (for lat/longs) or 'EN' (Easting/Northings)
%   'to'    : either 'LL' or 'EN'
%    'verbose' : display messages plus results of intermediate calculations
%
%    'warning' : alert user if |x| > |y|. Unlikely to be the case for our
%             areas of interest (user may have got input arguments wrong way round)
%
% OUTPUT:
%    x,y coordinates
%
% EXAMPLES:
%    [x,y]=catCoordinates(-1,52,'lut',lut,'csin','GPS','csout','NG','from','LL','to','EN') %will convert lat/longs from GPS to Eastings/Northings on the National Grid
%    [x,y]=catCoordinates(-1,52) % x = easting, y = northing. Same results as above (assumptions made about input values)
%    x=catCoordinates(-1,52) % x is now a 1x2 vector where x(1) = easting,  x(2) = northing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminology:
%
% CONVERSION is changing from Lat/Long to Easting/Northings (or vice versa)
% in the SAME coordinate system (eg WGS84)
% TRANSFORMATION is changing from from Lat/Long (or E/N) in one coordinate
% system to Lat/Long (or E/N) in another coordinate system.
%
% The transformation defined in this function is for Eastings/Northings only.
% So to transform Lat/Longs, it is necessary to convert them first!
%
% To change from GPS Lat/Longs, therefore, the following steps are
% required:
% 1) Convert from Lat/Long to Eastings/Northings in ETRS89 (WGS84)
%    coordinate system
% 2) Transform from Easting/Northing in ETRS89 to E/N in OSGB36
%
% In other words, both a conversion and a transformation are required.
%
% This function was developed using information obtained in the following
% documents, referred to as Document 1 and Document 2 in comments below:
%
% 1) "A guide to coordinate systems in Great Britain". This contains useful background
% information about transformations/conversions.
%
% 2) "Transformations and OSGM02 User Guide". This contains details of the calculations
% and worked examples in the appendices:
%   * Appendix A - transformation between ETRS89 and OSGB36 (done using
%     eastings/northings).
%   * Appendix B - conversion from lat/longs to Eastings/Northings
%   * Appendix C - conversion between E/N and Lat/Longs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING:
% This function was tested in several ways to ensure its accuracy:
%
% 1) The results from Appendix A of document 2 were recreated to 10 significant figures:
%    %%% Coordinates of Caister Water tower:
%    lonCWT=degreeConverter([1,42,57.8663],1);
%    latCWT=degreeConverter([52,39,28.8282],1);
%    %Step 1: convert LL to EN in ETRS89 coordinate system:
%    [eCWT_ETRS89,nCWT_ETRS89]=catCoordinates(lonCWT,latCWT,'csin','GPS','csout','GPS','warning',false) % matches values on p29
%    % Step 2: transform between coordinate systems:
%    [eCWT_NG,nCWT_NG]=catCoordinates(eCWT_ETRS89,nCWT_ETRS89,'to','EN','csin','GPS','csout','NG','warning',false) % matches values on p31
%    % Above steps in single command:
%    [eCWT_NG,nCWT_NG]=catCoordinates(lonCWT,latCWT,'warning',false)
%    % Inverse: OSGB36 to ETRS89
%    [eCWT_ETRS89,nCWT_ETRS89]=catCoordinates(eCWT_NG,nCWT_NG,'to','EN','csin','NG','csout','GPS','warning',false) % values on p32 (to nearest mm)
%
% 2) Appendix B of document 2 were reproduced using:
%       lat=degreeConverter([52,39,27.2531],1);
%       lon=degreeConverter([1,43,4.5177],1);
%       catCoordinates(lon,lat,'csin','NG','csout','NG','warning',false) % matches values on p 36
%
% 3) Appendix C of document 2:
%       e=651409.903;
%       n=313177.270;
%       ll=catCoordinates(e,n,'csin','NG','csout','NG','warning',false);
%       degreeConverter(ll',3)   % matches values on p 39
%
% 4) Converting lat/longs (GPS) to eastings/northings was also tested by
%    comparing with GridInQuest results
%       catCoordinates(-3,55)
%    Function call above agrees with GridInQuest conversion from ETRS89
%    Geodetic reference system (GPS) to OSGB 1936 reference system to
%    given precision:
%    Eastings  335128.92 (m)
%    Northings 567726.48 (m)
%
% 5) Converting eastings/northings to GPS lat/longs
%       catCoordinates(300000,600000)
%    produces:
%           Longitude: -3.57596954759802
%           Latitude: 55.2839446908619
%   These also agree with GridInquest (2002 version)
%   With 2015 OSTN look-up table, values are 
%       Longitude: -3.57596973537939 (-3.575969735411)         
%       Latitude:  55.2839446614774 (55.283944661484)
%   The GridInquest results are in brackets. Agreement to ~9 decimal places
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   catCoordinates.m  $
% $Revision:   1.7  $
% $Author:   Ted.Schlicke  $
% $Date:   Sep 25 2018 15:53:02  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help OS.catCoordinates
    return
end

% Create a structure that contains our user options:
options=struct;
% ... and assign some default values...
options.from=''; % Valid = LL, EN
options.to=''; % As above
options.output='numeric';
options.verbose=false;
options.warning=true;
options.version=2015;
options=checkArguments(options,varargin{:});

global verbose; % Lots of sub-functions in this function. We use a global variable so we don't need to keep passing it around
verbose=options.verbose;
if verbose
    fprintf('_____________________________________________________________\n')
    fprintf('*************************************************************\n')
    fprintf('WELCOME TO catCoordinates function!\n')
    fprintf('_____________________________________________________________\n')
end

% Prepare OSTN look up table based on version
% --- Prepare LUT ---
persistent lut02 lut15
offshoreWarned=false;
switch options.version
    case 2002
        if isempty(lut02), lut02 = readmatrix('OSTN02_OSGM02_GB.txt'); end
        lut = lut02;
    case 2015
        if isempty(lut15), lut15 = readmatrix('OSTN15_OSGM15_DataFile.txt'); end
        lut = lut15;
    otherwise
        error('Invalid version');
end

% Store original size so we can revert at end (e.g. if user passes 2d values)
origSize=size(x);
x=x(:);
y=y(:);

% Some error checking now...
if(length(x)~=length(y))
    error('X and Y data sets are of different lengths! Aborting...')
end

if options.warning && any(abs(x)>abs(y)) % False for Scotland't latitude - check user didn't get inputs muddled up!
    warning('OH:DEAR','Input |X| value > |Y| value - please check input arguments are right way round!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is input Eastings/Northings or Lat/Longs?
% Is required output E/N or L/L?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If user didn't specify 'from' or 'to', we'll make a guess:
if isempty(options.from)
    if all(abs(x)<=180) && all(abs(y)<=180) % These look like Lat/Longs!
        options.from='LL';
    else % Some big numbers...
        options.from='EN'; % most likely Eastings/Northings
    end
end
if isempty(options.to) % no 'to' option specified
    if strcmp(options.from, 'LL') % if input is lat/long,
        options.to='EN'; % we'll assume we're converting to EN
    else
        options.to='LL'; % otherwise, vice-versa
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now check coordinate systems...
% INPUT:
if strcmp(options.from,'EN') % E/N
    csin='NG'; % assume E/N in National Grid
else
    csin='GPS'; % and L/L in GPS
end

% And OUTPUT:
if strcmp(options.to,'EN') % E/N
    csout='NG'; % assume E/N in National Grid
else
    csout='GPS'; % and L/L in GPS
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OK, option settings should now be fixed...
if verbose
    fprintf('OPTIONS:\n')
    disp(options)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the procedures required. We store them in a cell array
% whose label corresponds to the relevant appendix of
% "Transformations and OSGM02 User Guide".
%
% Appendix A : Transformation of ETRS89 to OSGB36. Label = 'A'
%              Inverse Transformation - OSBG36 to ETRS89 (requires
%              iteration). Label = 'A2'
%
% Appendix B : Conversion of lat/long to easting/northing. Label = 'B'
%
% Appendix C : Conversion of easting/northing to lat/long. Label = 'C'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(options.from,options.to)&&strcmp(csin,csout))
    error('Input and Outputs are the same - nothing to convert!')
end
procedureSteps={};
% Transformation is done using E/N. So if we start with L/L we need to
% convert...
if(strcmp(options.from,'LL'))
    procedureSteps(length(procedureSteps)+1)={'B'}; % Conversion
    if verbose
        fprintf('Conversion from Longitude/Latitude to Easting/Northing\n');
    end
end
% If coordinate systems change, we need to do transformation (A or A2)
if(~strcmp(csin,csout)) % input coordinate system is different from output?
    if(strcmp(csin,'GPS')) % Input is GPS (WGS84) then output must be OSGB36
        procedureSteps(length(procedureSteps)+1)={'A'};
        if verbose
            fprintf('Transforming from WGS84 to OSGB36\n');
        end
    else
        procedureSteps(length(procedureSteps)+1)={'A2'}; % Input is OSGB36 and output =GPS (WGS84)
        if verbose
            fprintf('Transforming from OSGB36 to WGS84\n');
        end
    end
end
% And if output required is lat/long, we need to convert again
if(strcmp(options.to,'LL'))
    procedureSteps(length(procedureSteps)+1)={'C'}; %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose
    fprintf('Transformation/Conversion requires procedures from appendices = \n');
    disp(procedureSteps)
end

% Now we have our required procedure (combination of conversions and
% transformations).

% Prepare structures which conveniently contain the ellipsoid we need.
% First we declare a function to set the various parameters:
    function [ellipsoidConstants]=setEllipsoidConstants(a,b)
        ellipsoidConstants=struct;
        ellipsoidConstants.a=a;
        ellipsoidConstants.b=b;
        ellipsoidConstants.e2=(a^2-b^2)/a^2;
        ellipsoidConstants.n=(a-b)/(a+b);
    end
% And now we define them.
% We have the Airy 1830 ellipsoid, used for the OSGB36 (OS) coordinate system
airy1830Ellipsoid=setEllipsoidConstants(6377563.396,6356256.910);
% And the GRS80 (WGS84) Ellipsoid used for the ETRS89 (GPS) coordinate system
grs80Ellipsoid=setEllipsoidConstants(6378137.000,6356752.3141);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We also need the projection constants. Note that SEPA probably doesn't
% need the ITM (Irish Transverse Mercator) projection, so we'll just use to
% National Grid one
nationalGridProjection=struct;
nationalGridProjection.F0=0.9996012717;
nationalGridProjection.phi0=49*pi/180;
nationalGridProjection.lambda0=-2*pi/180;
nationalGridProjection.E0=400000;
nationalGridProjection.N0=-100000;

% Store converted/transformed values in these 1d arrays
% Default val is NaN, in case input values are invalid
N=length(x);
newx=NaN(N,1);
newy=NaN(N,1);

% Finally, we're ready to start looping through the points we need to
% convert/transform
for pointIndex=1:N % for each point in data set...
    xi=x(pointIndex);
    yi=y(pointIndex);
    if verbose
        fprintf('Processing point %d of %d :::: (%f,%f)..........\n',pointIndex,N,xi,yi);
    end
    if(~isnan(xi)&&~isnan(yi))
        % We've found the conversions/transformations we need to do.
        % Loop through them and do them!
        for procno=1:length(procedureSteps)

            % CHECK for transformation:
            step=char(procedureSteps(procno));
            if verbose
                fprintf('Procedure %d of %d from OS guide Appendix %s\n',procno,length(procedureSteps),step);
            end
            if ~isempty(regexp(step,'A', 'once')) %Step contains 'A' - it's a transformation!
                if(strcmp('A',step)) % exactly 'A' ? It's forward transformation
                    direction=1;
                else % it's A2 - backward transformation
                    direction=-1;
                end
                if verbose
                    fprintf('Step %s is a transformation; Direction = %d\n',step,direction)
                end
                if verbose
                    fprintf('Doing OSTN transformation...\n')
                end
                [xi,yi]=OSTNTransformation(xi,yi,lut,direction);

            else % It's not a transformation, it's a conversion!
                if verbose
                    fprintf('Step %s is a conversion\n',step)
                end
                % Need to determine what ellipsoid we're using for conversions...
                ellipsoidName=csout;
                if(procno==1)
                    ellipsoidName=csin;
                end
                % We have the ellipsoid name- now set its parameters
                if(strcmp(ellipsoidName,'GPS'))
                    ellipsoid=grs80Ellipsoid;
                elseif(strcmp(ellipsoidName,'NG'))
                    ellipsoid=airy1830Ellipsoid;
                else
                    error('Cant find relevant ellipsoid! Please use ''GPS'' or ''NG''Aborting...')
                end
                if verbose
                    fprintf('Using ellipsoid ''%s'':\n',ellipsoidName)
                    disp(ellipsoid)
                end
                if(strcmp(step,'B'))
                    [xi,yi]=coordinateConversion(xi,yi,ellipsoid,nationalGridProjection,1);
                    if(xi<0||xi>2e6||yi<0||yi>2e6)
                        %                        warning('Easting/Northings appear out of range- please check conversion options appropriate')
                    end
                elseif(strcmp(step,'C'))
                    [xi,yi]=coordinateConversion(xi,yi,ellipsoid,nationalGridProjection,-1);
                    if(xi<-90||xi>90||yi<-90||yi>90)
                        error('Long/Lats appear out of range- please check conversion options appropriate')
                    end
                end
            end
        end % end of procedures- so we have our values!
        newx(pointIndex)=xi;
        newy(pointIndex)=yi;
        if options.warning && abs(xi)>abs(yi)
            %            warning('OH:DEAR','OUTPUT |X| value > |Y| value !') % OK for for parts of England!
        end
    end % end of 'isnan' test
end

newx=reshape(newx,origSize);
newy=reshape(newy,origSize);

% Prepare output arguments
if nargout<2
    varargout={[newx,newy]};
elseif nargout==2% numeric output
    %            op=varargout;
    varargout{1}=newx;
    varargout{2}=newy;
else
    error('Too many outputs!')
end
return % END OF CALCULATION!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONVERSION FUNCTIONS - convert between lat/longs and x/y values
%%%%% within a datum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ newx,newy ] = coordinateConversion(x,y,ellipsoid,projection,direction)
if(direction>0) % Lat/Long to East/North...
    [x,y]=convertLL2EN(x,y,ellipsoid,projection);
else
    [x,y]=convertEN2LL(x,y,ellipsoid,projection);
end
newx=x;
newy=y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are lots of intermediate, messy calculations required. We'll bundle
% these in a separate function
    function [intermediateValues]=calculateIntermediateValues(phi,ellipsoidConstants,projectionConstants)
        intermediateValues=struct; % Store values in this structure
        % Abbreviations to reduce amount of typing
        a=ellipsoidConstants.a;
        b=ellipsoidConstants.b;
        e2=ellipsoidConstants.e2; % B1 calculation - already done!
        n=ellipsoidConstants.n;   % B2 calculation - done too!
        F0=projectionConstants.F0;
        phi0=projectionConstants.phi0;
        N0=projectionConstants.N0;

        %%%%%% B3 calculation
        nu=a*F0*(1-e2*sin(phi)^2)^(-0.5);
        intermediateValues.nu=nu;
        %%%%% B4
        rho=a*F0*(1-e2)*(1-e2*sin(phi)^2)^(-1.5);
        intermediateValues.rho=rho;
        %%%% B5
        eta2=nu/rho-1;
        intermediateValues.eta2=eta2;
        intermediateValues.M=calculateM(n,phi,phi0,b,F0);
        %%%%%%%%% Next set of calcs...
        intermediateValues.I=intermediateValues.M+N0;
        intermediateValues.II=nu/2*sin(phi)*cos(phi);
        intermediateValues.III=nu/24*sin(phi)*cos(phi)^3*(5-tan(phi)^2+9*eta2);
        intermediateValues.IIIA=nu/720*sin(phi)*cos(phi)^5*(61-58*tan(phi)^2+tan(phi)^4);
        intermediateValues.IV=nu*cos(phi);
        intermediateValues.V=nu/6*cos(phi)^3*(nu/rho-tan(phi)^2);
        intermediateValues.VI=nu/120*cos(phi)^5*(5-18*tan(phi)^2+tan(phi)^4+14*eta2-58*tan(phi)^2*eta2);
        %%%% Calculations for Appendix C
        intermediateValues.VII=tan(phi)/(2*rho*nu);
        intermediateValues.VIII=tan(phi)/(24*rho*nu^3)*(5+3*tan(phi)^2+eta2-9*tan(phi)^2*eta2);
        intermediateValues.IX=tan(phi)/(720*rho*nu^5)*(61+90*tan(phi)^2+45*tan(phi)^4);
        intermediateValues.X=sec(phi)/nu;
        intermediateValues.XI=sec(phi)/(6*nu^3)*(nu/rho+2*tan(phi)^2);
        intermediateValues.XII=sec(phi)/(120*nu^5)*(5+28*tan(phi)^2+24*tan(phi)^4);
        intermediateValues.XIIA=sec(phi)/(5040*nu^7)*(61+662*tan(phi)^2+1320*tan(phi)^4+720*tan(phi)^6);
        if verbose
            fprintf('Intermediate values:\n')
            disp(intermediateValues)
        end
    end

%Separate little function to calculate M (calculation B6).
% We keep this separate from the calculations above as it needs to be
% called iteratively - so should be faster this way!
    function [M]=calculateM(n,phi,phi0,b,F0)
        %%%%% B6 - big calc! do in bits... NB this could be simplified but left as
        %%%%% in OS document to reduce risk of typos/errors
        t1=(1+n+5/4*n^2+5/4*n^3)*(phi-phi0);
        t2=-(3*n+3*n^2+21/8*n^3)*sin(phi-phi0)*cos(phi+phi0); % use degree version of trig functions
        t3=(15/8*n^2+15/8*n^3)*sin(2*(phi-phi0))*cos(2*(phi+phi0));
        t4=-35/24*n^3*sin(3*(phi-phi0))*cos(3*(phi+phi0));
        M=b*F0*(t1+t2+t3+t4);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FUNCTION TO CONVERT Lat/Longs to Easting/Northings (APPENDIX B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [easting,northing]=convertLL2EN(longitude,latitude,ellipsoid,projection)
        % Use conventional variable names:
        lambda=longitude*pi/180;
        phi=latitude*pi/180;
        iv=calculateIntermediateValues(phi,ellipsoid,projection);

        % abbreviation:
        dl=(lambda-projection.lambda0);
        % Now we can calculate Easting and Northing!
        easting=projection.E0+iv.IV*dl+iv.V*dl^3+iv.VI*dl^5;
        northing=iv.I+iv.II*dl^2+iv.III*dl^4+iv.IIIA*dl^6;
        if verbose
            fprintf('End of Appendix B function (''convertLL2EN'') and Easting, Northing are: %f %f\n',easting,northing);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF APPENDIX B algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Function to Convert Eastings/Northings to Lat/Longs (APPENDIX C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [longitude,latitude]=convertEN2LL(easting,northing,ellipsoid,projection)
        % Abbreviations to reduce amount of typing
        a=ellipsoid.a;
        b=ellipsoid.b;
        n=ellipsoid.n;   % B2 calculation - done too!
        F0=projection.F0;
        phi0=projection.phi0;
        lambda0=projection.lambda0;
        E0=projection.E0;
        N0=projection.N0;
        E=easting;
        N=northing;
        %%%%%%%% C1 calculation:
        phiDash=(N-N0)/(a*F0)+phi0;
        phi=phiDash;
        %%%%%%%%%% ITERATIVE BIT:
        iteration=0;
        M=calculateM(n,phi,phi0,b,F0);
        disc=abs(N-N0-M);
        while disc>0.01
            phi=phi+(N-N0-M)/(a*F0);
            M=calculateM(n,phi,phi0,b,F0);
            disc=abs(N-N0-M);
            iteration=iteration+1;
            if iteration>1000 % Check in case we don't get convergence-
                error('No CONVERGENCE!!!')
                %%% might need to do something than just aborting should
                %%% we ever get here. But we shouldn't!
            end
        end

        iv=calculateIntermediateValues(phi,ellipsoid,projection);
        dE=E-E0;
        phi=phi-iv.VII*dE^2+iv.VIII*dE^4-iv.IX*dE^6;
        lambda=lambda0+iv.X*dE-iv.XI*dE^3+iv.XII*dE^5-iv.XIIA*dE^7;

        % Now we can find our answer!
        longitude=lambda*180/pi;
        latitude=phi*180/pi;
        if verbose
            fprintf('End of Appendix C function (''convertEN2LL'') and Longitude/Latitude are: %f %f\n',longitude,latitude);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Functions to TRANSFORM coordinates between datums
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ newx,newy ] = OSTNTransformation(x,y,lut,direction)
if(length(x)~=length(y))
    error('X and Y data sets are of different lengths! Aborting...')
end
if(length(x)~=1)
    error('Can only do one point at a time!')
end
index=1;
N=length(x);
newx=zeros(N,1);
newy=zeros(N,1);

for osti=1:N
    if(direction>0)
        [tmpx,tmpy]=transformETRS89toOSGB36(x(osti),y(osti),lut);
    else
        [tmpx,tmpy]=transformOSGB36toETRS89(x(osti),y(osti),lut);
    end
    newx(index)=tmpx;
    newy(index)=tmpy;
    index=index+1;
end

% is point within lookup table?
    function ok=isInBounds(row,lut)
        ok=true;
        if(row<1||row>size(lut,1))
            error('POINT IS OUT OF TRANSFORMATION AREA! Aborting...')
%            ok=false;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FUNCTION To Apply Shifts TRANSFORM ETRS89 to OSGB36
% Uses variable names as used on p30 of Document 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [se,sn,sg]=findOSShifts(easting,northing,lut)
        % Find square km within which this point lies on OS grid:
        east_index=floor(easting/1000);
        north_index=floor(northing/1000);
        se=nan;
        sn=nan;
        sg=nan;

        %%%%%%%%%% FIND SHIFTS
        % Start with Bottom Left:
        i=east_index;
        j=north_index;
        row=i+j*701+1;
        if(~isInBounds(row,lut))
            return
        end
        se0=lut(row,4);
        sn0=lut(row,5);
        sg0=lut(row,6);

        % Bottom Right
        i=i+1;
        row=i+j*701+1;
        if(~isInBounds(row,lut))
            return
        end
        se1=lut(row,4);
        sn1=lut(row,5);
        sg1=lut(row,6);

        % Top Right
        j=j+1;
        row=i+j*701+1;
        if(~isInBounds(row,lut))
            return
        end
        se2=lut(row,4);
        sn2=lut(row,5);
        sg2=lut(row,6);

        % Top Left
        i=i-1;
        row=i+j*701+1;
        if(~isInBounds(row,lut))
            return
        end
        se3=lut(row,4);
        sn3=lut(row,5);
        sg3=lut(row,6);

        % isInBounds may return true even if all shifts are 0
        if se0==0 && se1==0 && se2==0 && se3==0 && sn0==0 && sn1==0 && sn2==0 && sn3==0 && ~offshoreWarned
            warning('Points outside OSTN coverage — zero shift applied.');
            offshoreWarned=true;
        end
        % Now we do bilinear interpolation to shift our E/N vals.
        % Some abbreviations first
        dx=easting-east_index*1000;
        dy=northing-north_index*1000;
        t=dx/1000;
        u=dy/1000;
        % now get shifts...
        se=(1-t)*(1-u)*se0+t*(1-u)*se1+t*u*se2+(1-t)*u*se3;
        sn=(1-t)*(1-u)*sn0+t*(1-u)*sn1+t*u*sn2+(1-t)*u*sn3;
        sg=(1-t)*(1-u)*sg0+t*(1-u)*sg1+t*u*sg2+(1-t)*u*sg3;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
%%%% THIS function transforms ETRS89 to OSGB using steps described in appendix A
%%%% It simply finds the shifts using the function above, and adds them on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [easting,northing,h]=transformETRS89toOSGB36(easting,northing,lut)
        h=0;
        [se,sn,sg]=findOSShifts(easting,northing,lut);
        % F apply these shifts to NG eastings and northings:
        easting=easting+se;
        northing=northing+sn;
        %            fprintf('Shifted by %f, %f\n',se,sn)
        h=h-sg;
        if verbose
            fprintf('End of Appendix A routine and Easting,Northing : %f %f\n',easting,northing);
            fprintf('height = %f\n',h)
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% This function does transformation from OSGB36 to ETRS89 using steps described in appendix A.
%%%%% Slightly more involved than reverse direction as we need to iterate...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [easting,northing,h]=transformOSGB36toETRS89(easting,northing,lut)
        h=0;
        seOld=0;
        snOld=0;
        disc=1;
        cutoff=1e-8;
        osEasting=easting;
        osNorthing=northing;
        % Iterative until change in shift is negligible (less than 'cutoff')
        counter=1;
        while disc>cutoff
            [se,sn,sg]=findOSShifts(easting,northing,lut);
            disce=abs(se-seOld);
            discn=abs(sn-snOld);
            if(disce>discn)
                disc=disce;
            else
                disc=discn;
            end
            seOld=se;
            snOld=sn;
            easting=osEasting-se;
            northing=osNorthing-sn;
            counter=counter+1;
        end
        h=h+sg;
        if verbose
            fprintf('End of Appendix A (inverse) routine and Easting,Northing : %f %f\n',easting,northing)
            fprintf('Orthometric height = %f\n',h)
            fprintf('Took %d iterations\n',counter)
        end
    end
end
end
