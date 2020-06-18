function outPt = isInEllipse(xyzCoord,ellSz,ellCent,ellAng)

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