%% This code generates ellipsoids at different locations in the imaging 
% domain. The 3D x matrix contains information about these ellipsoids as a
% binary value, 1 if that coordinate falls entirely inside ellipsoid, 0
% otherwise. 

function [objCenterGrid] = genLibraryTrueSizeMP(opts)
% -------------------------------------------------------------------------
% Input: 
% opts: Different input options in this structure. the default options are
% assigned if any particular field is not defined.
% -------------------------------------------------------------------------
% Output: 
% Alib: The new A for the generated library cases.
% X: Different images in column, with each column corresponding to a column
% of Alib.
% -------------------------------------------------------------------------
% The specifications and arrangement of library are shown as below: 
% ---------------------------SPECIFICATIONS--------------------------------
% Scaling is nearly 1:1. Quasi-2D simulation.
% Simulation size: (approx.) 3.6 m x 3.6 m x 0.6 m. 
% The coordinates are between: [0.2, 3.4; 0.2, 3.4; -0.25, 0.25] m for x, y 
% and z respectively. 
% Voxel size is [0.02;0.02;0.02] m; 
% Object1: Rectangular object.
% -----------------------------ARRANGEMENT---------------------------------
% Considering one quadrant in XY plane, 5 x 5 equally distributed grid
% points.
% -------------------------------------------------------------------------
% Pragya Sharma, ps847@cornell.edu, 
% Created on: 12/18/2018
% Last update on: 05/07/2019: For true scale library generation. Defaults
% will change. Object is simplified to be cuboidal.
% -------------------------------------------------------------------------

%% Input definitions, XYZ voxel coordinates and A matrix for a particular
% location of Tx-Rx antenna and frequencies.
if ~isfield(opts,'imSize')
    opts.imSize = [0.2, 3.4; 0.2, 3.4; -.25, 0.25]; % meters
end
if ~isfield(opts,'voxSz')
    opts.voxSz = [0.02;0.02;0.02]; % meters
end
if ~isfield(opts,'freq')
    opts.freq = [905500000;   910500000;    915500000;    920500000;...
                 925500000;   930500000];
end
if ~isfield(opts,'posRxTx')
    opts.posRxTx = 1;
end
if ~isfield(opts,'A')
    [~,~,rxPosition,tagPosition] = posRxTxTrueSize(opts.posRxTx);
    opts.A = genA(opts.imSize,opts.voxSz,tagPosition,rxPosition,opts.freq);
end
if ~isfield(opts,'savePath')
    opts.savePath = ['E:\ArpaE2018\3DImaging_Simulation\CST_Simulation',...
        'DataAnalysis\Algorithms\MP\LibATrueSize\'];
end
if ~isfield(opts,'fileName')
    fprintf('File Name is not input, keeping a new fileName as: lib_mmddyyyy_HHMM\n');
    opts.fileName = sprintf('lib_%s', datestr(now,'mmddyyyy_HHMM'));
end
if ~isfield(opts,'seeImg')
    opts.seeImg = 0; % Default visualization
end
if ~isfield(opts,'objSz')
    opts.objSz = [0.08 0.06 0.3];
end

%% Generating different possible object cases and corresponding image

[xyzVoxCoord,voxImCoord,xyzSlice,nVoxel] = genXYZ(opts.imSize,opts.voxSz);

% X = [];
Alib = [];
arrGrid = @arrGrid3; % Different grid arrangement functions
genImg = @genImgRect; % Different image generation

% First case:  grid with Object 1
objCenterGrid = arrGrid(opts.imSize,opts.objSz);
[opts.objSz,opts.objOrient] = arrObjPara(opts.objSz,1);
for i = 1:size(objCenterGrid,1)
    imgX = genImg(nVoxel,xyzVoxCoord,voxImCoord,opts.objSz,...
        objCenterGrid(i,:),opts.objOrient,xyzSlice,opts.imSize,opts.seeImg);
    Atemp = opts.A*imgX(:);
    Alib = [Alib, Atemp]; %#ok<*AGROW>
%     X = [X, imgX];
end

save([opts.savePath,opts.fileName,'.mat'],'Alib','objCenterGrid');

end

%% Generate images with different shapes
% Generating filled ellipsoidal objects
function imgX = genImgEllip(nVoxel,xyzVoxCoord,voxImCoord,objSz,objCenter,objOrient,xyzSlice,roomSize,seeImg)
% This function generates binary image with 1, if any voxel coordinate is
% inside object with given size, center and orientation
nVoxTot = nVoxel(1)*nVoxel(2)*nVoxel(3);
imgX = zeros(nVoxTot,1);
for i = 1:nVoxTot
    imgX(i) = isInEllipse(xyzVoxCoord(i,:),objSz,objCenter,objOrient);
end

if seeImg == 1
    visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize);
end

end

% Generate only an ellipsoidal shell
function imgX = genImgEllipShell(nVoxel,xyzVoxCoord,voxImCoord,objSz,objCenter,objOrient,xyzSlice,roomSize,seeImg)
% This function generates binary image with 1, if any voxel coordinate is
% inside object with given size, center and orientation
nVoxTot = nVoxel(1)*nVoxel(2)*nVoxel(3);
imgX = zeros(nVoxTot,1);
innerThick = [0.65 0.6 0.8];
for i = 1:nVoxTot
    inOuterEllipse = isInEllipse(xyzVoxCoord(i,:),objSz,objCenter,objOrient);
    inInnerEllipse = isInEllipse(xyzVoxCoord(i,:),objSz.*innerThick,objCenter,objOrient);
    imgX(i) = (inOuterEllipse)&&(~inInnerEllipse);
end

if seeImg == 1
    visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize);
end
end

% Generating cuboidal objects
function imgX = genImgRect(nVoxel,xyzVoxCoord,voxImCoord,objSz,objCenter,~,xyzSlice,roomSize,seeImg)
objSz = objSz.*0.9;
nVoxTot = nVoxel(1)*nVoxel(2)*nVoxel(3);
imgX = zeros(nVoxTot,1);
for i = 1:nVoxTot
    imgX(i) = isInCuboid(xyzVoxCoord(i,:),objSz,objCenter);
end

if seeImg == 1
    visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize);
end
end

% Is point inside ellipsoid
function outPt = isInEllipse(xyzCoord,ellSz,ellCent,ellAng)
% ellSz indicates axes of the ellipsoid, NOT semi-axes, hence
% multiplication with 0.5
outPt = ((cos(ellAng)*(xyzCoord(1)-ellCent(1)) + ...
    sin(ellAng)*(xyzCoord(2)-ellCent(2)))/(0.5*ellSz(1)))^2 + ...
    ((sin(ellAng)*(xyzCoord(1)-ellCent(1)) - ...
    cos(ellAng)*(xyzCoord(2)-ellCent(2)))/(0.5*ellSz(2)))^2 + ...
    ((xyzCoord(3)-ellCent(3))/(0.5*ellSz(3)))^2;
if outPt <= 1
    outPt = 1;
else
    outPt = 0;
end

end

% Is point inside cuboid
function outPt = isInCuboid(xyzCoord,objSz,objCent)

xMinMax = [objCent(1)-(objSz(1)/2), objCent(1)+(objSz(1)/2)];
yMinMax = [objCent(2)-(objSz(2)/2), objCent(2)+(objSz(2)/2)];
zMinMax = [objCent(3)-(objSz(3)/2), objCent(3)+(objSz(3)/2)];

testx = (xyzCoord(1) <= xMinMax(2)) && (xyzCoord(1) >= xMinMax(1));
testy = (xyzCoord(2) <= yMinMax(2)) && (xyzCoord(2) >= yMinMax(1));
testz = (xyzCoord(3) <= zMinMax(2)) && (xyzCoord(3) >= zMinMax(1));

outPt = testx && testy && testz;
end

%%
function visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize)

p = [2 1 3];
imgX = reshape(imgX,nVoxel);
imgX = permute(imgX,p);

x = reshape(voxImCoord(:,1),nVoxel(2),nVoxel(1),nVoxel(3));
y = reshape(voxImCoord(:,2),nVoxel(2),nVoxel(1),nVoxel(3));
z = reshape(voxImCoord(:,3),nVoxel(2),nVoxel(1),nVoxel(3));


figure
h = slice(x,y,z,imgX,...
    xyzSlice{1},xyzSlice{2},xyzSlice{3});
xlabel('x-axis','FontSize',14)
ylabel('y-axis','FontSize',14)
zlabel('z-axis','FontSize',14)
xlim([roomSize(1,1), roomSize(1,2)])
ylim([roomSize(2,1), roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')

end

%% Arranging object grid
function objCenter = arrGrid1(imSz,objSz)
% This function generates 5 x 5 grid points in first quadrant of the
% room, with z = zMin + (objHeight/2);
objSz = sort(objSz,'ascend');
zCoord = imSz(3,1) + (objSz(3)/2);
xLim = [(imSz(1,1)+imSz(1,2))/2, imSz(1,2)-objSz(2)];
yLim = [(imSz(2,1)+imSz(2,2))/2, imSz(2,2)-objSz(2)];
numPtsXY = [5, 5];
xCoord = linspace(xLim(1),xLim(2),numPtsXY(1));
yCoord = linspace(yLim(1),yLim(2),numPtsXY(2));

objCenter = combvec(xCoord,yCoord,zCoord)';

end

function objCenter = arrGrid2(imSz,objSz)
% This function generates 4 x 4 grid points in first quadrant of the
% room, with z = zMin + (objHeight/2);
objSz = sort(objSz,'ascend');
zCoord = imSz(3,1) + (objSz(3)/2);
xLim = [(imSz(1,1)+imSz(1,2))/2, imSz(1,2)-objSz(2)];
yLim = [(imSz(2,1)+imSz(2,2))/2, imSz(2,2)-objSz(2)];
numPtsXY = [4, 4];
xCoord = linspace(xLim(1),xLim(2),numPtsXY(1));
yCoord = linspace(yLim(1),yLim(2),numPtsXY(2));

objCenter = combvec(xCoord,yCoord,zCoord)';

end

function objCenter = arrGrid3(imSz,objSz)
% This function generates 20 x 20 grid points in entire room, with z = 0
%objSz = sort(objSz,'ascend');
zCoord = 0;
xLim = [(imSz(1,1)+objSz(1)), imSz(1,2)-objSz(1)];
yLim = [(imSz(2,1)+objSz(2)), imSz(2,2)-objSz(2)];
numPtsXY = [20, 20];
xCoord = linspace(xLim(1),xLim(2),numPtsXY(1));
yCoord = linspace(yLim(1),yLim(2),numPtsXY(2));

objCenter = combvec(xCoord,yCoord,zCoord)';

end

%%
function [objSz,objOrient] = arrObjPara(objSz,caseNum)
objSzReorder = sort(objSz);
objOrient = 0;
switch(caseNum)
    case 1
        % Same as input
        objSz = [objSz(1), objSz(2), objSz(3)];
    case 2
        % Arrange with major axis = z, next axis = x, smallest, y.
        objSz = [objSzReorder(2), objSzReorder(1), objSzReorder(3)];
        objOrient = pi/4;
    otherwise
        fprintf('Please enter valid object arrangement number. \n');
end
end
