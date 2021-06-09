function varargout = getWaveExposureIndex(x,y,waveExposureStruct)
% Find the wave exposure index for specified x,y coordinates
%
% Wave Exposure values were extracted from .tif file here:
% \\sepa-fp-01\DIR SCIENCE\EAU\GIS\Data\Marine_Scotland\Wave_Exposure
%
% See also 
% http://sepa-app-gis01/interactivemapdrilldown/Metadata.aspx?strMetaDatasetID=1287
%
% INPUTS: 
%   Easting(s)
%   Northings(s) 
%
% OPTIONAL INPUT:
%   waveExposureStruct (see note below)
%
% OUTPUTS:
%   Wave Exposure Index(ices)
%   Distance(s) from query point(s) to nearest point in waveExposureStruct
%
% EXAMPLE:
% waveExposureIndex=getWaveExposureIndex(143190,780486)
% [waveExposureIndex,dist]=getWaveExposureIndex(143190,780486) % return distance as well)
%
% Note: waveExposureStruct generated in 'waveExposureStuff20200107' script.
% It should be in the Data folder on the VDrive. If it's not passed as an
% argument, the function will try to load it. Fractionally faster to pass
% it directly. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   getWaveExposureIndex.m  $
% $Revision:   1.0  $
% $Author:   Ted.Schlicke  $
% $Date:   Jan 07 2020 15:28:42  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    help OS.getWaveExposureIndex
    return
end

if ~exist('waveExposureStruct','var') % Not passed to function? 
    try % Then try to load it. 
        % Load then extract struct. That's more than twice as fast as:
        %waveExposureStruct=importdata('waveExposureStruct.mat')
        tmp=load('waveExposureStruct.mat');
        waveExposureStruct=tmp.waveExposureStruct;       
    catch
        error('Struct containing wave exposure data not found!')
    end
end
% Find closest point in waveExposureStruct to test points:
[dist,k]=distanceBetweenPoints(x,y,waveExposureStruct.x,waveExposureStruct.y,'min');
waveExposureIndex=waveExposureStruct.waveExposure(k);

if nargout==0
   fprintf('Checked %d points\n',length(x));
   fprintf('Wave index / distance:\n')
   disp([waveExposureIndex,dist])
   
else

if nargout>0
    varargout{1}=waveExposureIndex;
end
if nargout>1
    varargout{2}=dist; % Might want to know how far apart closest point is
end
if nargout>2
    error('too many outputs!')
end
end


end
