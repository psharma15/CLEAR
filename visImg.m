function imgBrightness = visImg(imgComplex,roomSize,voxelSize)

xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
% Combination of all of these to get coordinates for all the voxels
xyzVoxelCoord = combvec(xVoxel,yVoxel,zVoxel)';

nx = length(xVoxel);
ny = length(yVoxel);
nz = length(zVoxel);

imgBrightness = (abs(imgComplex).^2);
maxBrightness = max(imgBrightness(:)); 
minBrightness = min(imgBrightness(:));
imgBrightness = (imgBrightness-minBrightness)/(maxBrightness-minBrightness);

% Visualizing reconstructed image
x = reshape(xyzVoxelCoord(:,1),nx,ny,nz);
y = reshape(xyzVoxelCoord(:,2),nx,ny,nz);
z = reshape(xyzVoxelCoord(:,3),nx,ny,nz);

xslice = roomSize(1,1):voxelSize(1):roomSize(1,2);
yslice = roomSize(2,1):voxelSize(2):roomSize(2,2);
zslice = roomSize(3,1):voxelSize(3):roomSize(3,2);

p = [2 1 3];
x = permute(x, p);
y = permute(y, p);
z = permute(z, p);
imgBrightness = reshape(imgBrightness,[nx ny nz]);
imgBrightness = permute(imgBrightness, p);

% figure('Position',[100,100,180,145]);
figure('Position',[400,400,450,350]);

fSz=7;
h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
xlabel('x (m)','FontSize',fSz)
ylabel('y (m)','FontSize',fSz)
zlabel('z (m)','FontSize',fSz)
xlim([roomSize(1,1),  roomSize(1,2)])
ylim([roomSize(2,1),  roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')

end