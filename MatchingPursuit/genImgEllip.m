function imgX = genImgEllip(nVoxel,xyzVoxCoord,voxImCoord,objSz,objCenter,objOrient,xyzSlice,roomSize)
% This function generates binary image with 1, if any voxel coordinate is
% inside object with given size, center and orientation
nVoxTot = nVoxel(1)*nVoxel(2)*nVoxel(3);
imgX = zeros(nVoxTot,1);
for i = 1:nVoxTot
    imgX(i) = isInEllipse(xyzVoxCoord(i,:),objSz,objCenter,objOrient);
end

visImg1(imgX,nVoxel,voxImCoord,xyzSlice,roomSize);

end

