function op=closestTotalTidePortsBySea(longitude,latitude,varargin)
% Function to find closest TotalTide port to specfied point by seaway distance.
% 'As the crow flies' distance between two points may cross land.
% This function determines 'As the fish swims' distance.
%
% INPUTS:
% x0 - longitude
% y0 - latitude
%
% Optional Inputs:
% maxDistance (5e4) - distance (m) within which to search
% plot (true) - plot ports and paths
% mesh ([]) - model mesh used to determine paths
% closest (true) - only return closest port
% portOnly (false) - return TotalTide port(s)
% nameOnly (false) - just return port name (otherwise return name and distance in struct)
%
% OUPUT - either:
% struct with fields Name, Distance, Port
% port - (if portOnly==true) TotalTide port which can be passed to other function in TotalTide library
% name - (if portOnly==false and nameOnly==true) - name of TotalTide port
%
% This function uses Dijkstra's algorithm for finding shortest path between
% two points on a graph of nodes. In this case, the nodes are specified by
% the wider domain Scottish Shelf Model.
%
% We crop the model domain to our region of interest which speeds things
% up considerably. Algorithm can be quite slow with 10s of thousands of
% nodes!

options=struct;
options.maxDistance=5e4; % distance from point to look, metres
options.plot=0;
options.mesh=[];
options.closest=true;
options.nameOnly=false;
options.close=true;
options.portOnly=false;
options.heightPredicting=true;
options=checkArguments(options,varargin);

if latitude<90
    [x0,y0]=OS.catCoordinates(longitude,latitude);
else % assume input is easting/northing
    x0=longitude;
    y0=latitude;
end

%fprintf('Lat/long = %f,%f\n',longitude,latitude)
%fprintf('E/N = %f,%f\n',x0,y0)

% Check SSM Mesh
mikeMesh=options.mesh;
if isempty(mikeMesh)
    try
        f='\\asb-fp-mod01\AMMU\MarineModelling\MarineModelLibrary\ScottishShelfModel\SSM\ssmEN.mat';
        mikeMesh=importdata(f);
    catch
        error('Problem loading mike mesh :-(')
    end
end

%% Abbreviations
xMesh=mikeMesh.xMesh;
yMesh=mikeMesh.yMesh;
tri=mikeMesh.triMesh;
% Polyshape for plotting
ssmPolyshape=polyshape(mikeMesh.xMeshBoundary,mikeMesh.yMeshBoundary,'simplify',0);

% All total tide ports - stored in variable 'totalTidePortInfo' which is created by TotalTide.portFinder to speed things up
try
    totalTidePortInfo=evalin('base','totalTidePortInfo');
catch
    % oh well
end
if ~exist('totalTidePortInfo','var')
    %    fprintf('calling portFinder function\n')
    [~]=TotalTide.portFinder('dunbar');
    totalTidePortInfo=evalin('base','totalTidePortInfo');
end
listOfTotalTidePorts=evalin('base','listOfTotalTidePorts');
localPortInterfaces=arrayfun(@(x)listOfTotalTidePorts.Item(x),1:listOfTotalTidePorts.Count);

if options.heightPredicting
    % Remove ports where MeanSeaLevel is a char:
    % 'Invoke Error, Dispatch Exception:  Description: The data for this port does not allow height predictions'
    % Can't use these for water level predictions!
    % NB 20221216 This fails for Machrihanish- its MSL = 1.01, but other
    % heights are chars!
    %    k=cellfun(@ischar,{totalTidePortInfo.MeanSeaLevel});
    k=cellfun(@ischar,{totalTidePortInfo.MeanSeaLevel})|cellfun(@ischar,{totalTidePortInfo.HighestHighWater});

    totalTidePortInfo(k)=[];
    localPortInterfaces(k)=[];
end
xTotalTideAllPorts=[totalTidePortInfo.Longitude];
yTotalTideAllPorts=[totalTidePortInfo.Latitude];
% This causes warning due to points being outwith conversion domain:
[xTotalTideAllPorts,yTotalTideAllPorts]=OS.catCoordinates(xTotalTideAllPorts,yTotalTideAllPorts,'warn',0);

% Shortest paths within cropped model domain

% Define circle about specified coordinates
ang=linspace(0,2*pi,100);
r=options.maxDistance; % radius within which to check
xc=x0+r*cos(ang);
yc=y0+r*sin(ang);

% Find total tide ports within circle
k=inpolygon(xTotalTideAllPorts,yTotalTideAllPorts,xc,yc);
if ~any(k)
    error('No ports within range %.0fm - please increase maxDistance input!',r)
end
localTotalTidePortInfo=totalTidePortInfo(k);
localPortInterfaces=localPortInterfaces(k);
xTotalTideLocalPorts=[localTotalTidePortInfo.Longitude];
yTotalTideLocalPorts=[localTotalTidePortInfo.Latitude];
[xTotalTideLocalPorts,yTotalTideLocalPorts]=OS.catCoordinates(xTotalTideLocalPorts,yTotalTideLocalPorts);

% Find model mesh vertices within distance r of point of interest
k=inpolygon(xMesh,yMesh,xc,yc);

% % REINDEX TRIANGULATION
ktri=k(tri); % triangulation of points within circle
triSum=sum(ktri,2); % number of points of triangle within circle
k=triSum>2; % triangles to keep
cropTri=tri(k,:); % cropped triangulation
vk=unique(cropTri); % Ordered indices of cropped triangulation
ind=1:length(vk); % new indices, starting at 1
try
    cropTriChangem=changem(cropTri,ind,vk); % Renumbered triangulation
catch err
    disp(err)
    error('Oh dear, changem function from mapping toolbox required')
end
% Filter nodes, keeping those used in triangulation
xk=xMesh(vk);
yk=yMesh(vk);
xy=[xk,yk];

% Node connections of cropped triangulation
I=cropTriChangem; % Triangulation of interest
J = I(:,[2 3 1]); % Shifted triangulation
A = [I(:) J(:)];  % Bundled connections
Nnodes=length(xk);
% For each node, find all connections:
nodeConnections=cell(Nnodes,1); % space to store them
for nodeIndex=1:Nnodes % loop through nodes
    [rowIndex,~]=find(ismember(I,nodeIndex));
    iConnections=setdiff(unique(I(rowIndex,:)),nodeIndex);
    nodeConnections{nodeIndex}=iConnections(:);
end

% Find index of vertex closest to point of interest:
[~,ind0]=distanceBetweenPoints(x0,y0,xk,yk,'min');
% And indices of vertices closest to TotalTide ports:
[~,indTotalTide]=distanceBetweenPoints(xTotalTideLocalPorts,yTotalTideLocalPorts,xk,yk,'min');
%xy=[mikeMesh.xMesh(indTotalTide),mikeMesh.yMesh(indTotalTide)];

% Find paths to total tide node indices from single starting point:
[nodeDistance,nodeConnections]=dijkstra(xy,A,ind0,indTotalTide,false);

% Reorder connections and ports in order of increasing distance
[nodeDistance,costIndex]=sort(nodeDistance); % find indices defining order
nodeConnections=nodeConnections(costIndex); % reorder connections
localTotalTidePortInfo=localTotalTidePortInfo(costIndex); % and total tide ports
localPortInterfaces=localPortInterfaces(costIndex);

% Remove any with infinite distance (that can't be reached):
% 20250205 - Problem with location in River Carron - returns Clydebank?!?!?
% Cost == 0. So remove these as well. Not sure what issue is
k=nodeDistance==Inf | nodeDistance==0;
nodeConnections(k)=[];
nodeDistance(k)=[];
localTotalTidePortInfo(k)=[];
localPortInterfaces(k)=[];
if options.plot
    % Get ready to plot:
    N=length(nodeDistance);
    cm=rainbow(N);
    % Now plot model mesh within domain and connections
    str=sprintf('Height Predicting = %d',options.heightPredicting);
    prepareFigure('close',options.close,'title',str)
    plot(ssmPolyshape,'facecolor','w','facealpha',1)
    set(gca,'color',[94, 127, 57]/256)
    trimesh(cropTriChangem,xk,yk,xk*0,'edgecolor','k','facealpha',0)
    hPoint=scatter(x0,y0,100,'c','filled','displayname','Input coordinates');

    % Plot Total Tide coordinates within domain:
    %scatter(ttx,tty,30,'b','filled')
    % Plot closest vertices to Total Tide coordinates:
    scatter(xk(indTotalTide),yk(indTotalTide),30,'m','filled')
    ax=boundaryPolygon(xc,yc,'dx',0.1,'axis',1);
    axis(ax)

    N=length(nodeConnections);
    ok=true(N,1);
    hPath=nan(N,1);
    for pathIndex=1:N
        con=nodeConnections{pathIndex};
        try
            xi=xk(con);
            yi=yk(con);
            col=cm(pathIndex,:);
            tti=localTotalTidePortInfo(pathIndex);
            hPath(pathIndex)=plot(xi,yi,'Color',col,'linewidth',3,'DisplayName',tti.Name);
            %            fprintf('Port %25s and distance = %.0fm\n',tti.Name,nodeDistance(pathIndex))
        catch err
            disp(err)
            fprintf('Problem!\n')
            ok(pathIndex)=false;
            % oh well
        end
        %    pause
    end
    daspect([1,1,1])
    ax=gca;
    ax.XRuler.Exponent = 0;
    ax.YRuler.Exponent = 0;
    xtickformat('%8.f');
    ytickformat('%8.f');
    h4Legend=[hPoint;hPath];
    leg=legend(h4Legend);
    set(leg,'location','northeastoutside','color',[1,1,1])
    localTotalTidePortInfo(~ok)=[];
end

% Prepare output struct
%portNames={localTotalTidePortInfo(ok).Name};
%dist=nodeDistance(ok);
portNames={localTotalTidePortInfo.Name};
dist=nodeDistance;

%op=struct('Name',{transpose(portNames)},'Distance',dist','Port',localTotalTidePorts)'; % struct
op=struct('Name',portNames,'Distance',num2cell(dist),'Port',num2cell(localPortInterfaces),'PortInfo',num2cell(localTotalTidePortInfo)'); % struct array

if options.closest
    op=op(1);
end
if options.portOnly
    op=[op.Port];
    return
end
if options.nameOnly
    op={op.Name};
    if length(op)==1
        op=char(op);
    end
end

