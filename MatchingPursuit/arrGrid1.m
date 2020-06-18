function objCenter = arrGrid1(imSz,objSz)
% This function generates 5 x 5 grid points in first quadrant of the
% room, with z = zMin + (objHeight/2);
objSz = sort(objSz,'ascend');
zCoord = imSz(3,1) + (objSz(3)/2);
xLim = [(imSz(1,1)+imSz(1,2))/2, imSz(1,2)-objSz(2)];
yLim = [(imSz(2,1)+imSz(2,2))/2, imSz(2,2)-objSz(2)];

xCoord = linspace(xLim(1),xLim(2),5);
yCoord = linspace(yLim(1),yLim(2),5);

objCenter = combvec(xCoord,yCoord,zCoord)';

end
