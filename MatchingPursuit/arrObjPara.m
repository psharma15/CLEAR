function [objSz,objOrient] = arrObjPara(objSz,caseNum)
objSzReorder = sort(objSz);
objOrient = 0;
switch(caseNum)
    case 1
        % Arrange with major axis = z, next axis = x, smallest, y.
        objSz = [objSzReorder(2), objSzReorder(1), objSzReorder(3)];
    case 2
        % Arrange with major axis = z, next axis = x, smallest, y.
        objSz = [objSzReorder(2), objSzReorder(1), objSzReorder(3)];
        objOrient = pi/4;
    otherwise
        fprintf('Please enter valid object arrangement number. \n');
end
end
