function [xyzVoxelCoord,xyzImCoord,xyzSlice,nVoxel] = genXYZ(roomSize,voxelSize)
%% Generating (x,y,z) coordinates
% The first two outputs may be confusing, refer to other sections of the
% code how they are used. The first one is used in generating A matrix.
% After image is generated and reshaped, it's organization changes for the
% plotting purposes ONLY, and the second definition is used.

xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
% Combination of all of these to get coordinates for all the voxels
xyzVoxelCoord = combvec(xVoxel,yVoxel,zVoxel)';

nVoxel = [length(xVoxel), length(yVoxel), length(zVoxel)];

x = reshape(xyzVoxelCoord(:,1),nVoxel(1),nVoxel(2),nVoxel(3));
y = reshape(xyzVoxelCoord(:,2),nVoxel(1),nVoxel(2),nVoxel(3));
z = reshape(xyzVoxelCoord(:,3),nVoxel(1),nVoxel(2),nVoxel(3));

xslice = roomSize(1,1):voxelSize(1):roomSize(1,2);
yslice = roomSize(2,1):voxelSize(2):roomSize(2,2);
zslice = roomSize(3,1):voxelSize(3):roomSize(3,2);

xyzSlice = {xslice,yslice,zslice};

p = [2 1 3];
x = permute(x, p);
y = permute(y, p);
z = permute(z, p);

xyzImCoord = [x(:), y(:), z(:)];

end