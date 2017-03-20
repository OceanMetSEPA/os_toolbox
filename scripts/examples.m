%% Convert OSGB easting/northing coordinates into lat/lon

[lon, lat] = OS.convertAndTransform(219055,656439)

% lon =
%          -4.88586469639923
% lat =
%            55.767345888062

%% Convert lat/lon coordinates into OSGB easting/northing
           
[easting, northing] = OS.convertAndTransform(-4.88586469639923, 55.767345888062)

% easting =
%           219055.000044897
% northing =
%           656438.999978426

%% Explicitly state the conversion types

% In many cases the use of lat/lon versus easting/northing can be detected automatically
% by the magnitude of the numbers used. This can be explicitly controlled though by using the 
% from and to optional arguments 

[lon, lat]          = OS.convertAndTransform(219055,656439, 'from', 'EN', 'to', 'LL')
[easting, northing] = OS.convertAndTransform(-4.88586469639923, 55.767345888062, 'from', 'LL', 'to', 'EN')

%% Convert an all-numeric OSGB easting and northing into a National Grid Reference

ngr = OS.num2ngr(219055,656439)
 
% ngrCode =
% NS 19055 56439

% The same thing can also be achieved with separate outputs

[ngrCode, ngrEasting, ngrNorthing] = OS.num2ngr(219055,656439)

% ngrCode =
% NS
% ngrEasting =
%        19055
% ngrNorthing =
%        56439

%% Convert a National Grid Reference string into an all-numeric OSGB easting and northing 

[e,n] = OS.ngr2num('NS', '19055', '56439')

% e =
%       219055
% n =
%       656439

% The same thing can be achieved by passing in a single string

[e,n] = OS.ngr2num('NS 19055 56439')

% e =
%       219055
% n =
%       656439

%%
