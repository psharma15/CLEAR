% Processing centroids, combining them, if distance is small

function clusterOut = clusterProcess(clusters,opts)

if ~isfield(opts,'distTh')
    opts.distTh = 0;
end
if ~isfield(opts,'distTh')
    opts.XYdistTh = 0; % Only cares if x-y are close
end

if ~isfield(opts,'elemNumTh')
    % Percent value of the maximum number of elements. Only centroids
    % greater than (opts.elemNumTh * max_elem_num) will be kept. 
    % 0 = no centroid rejected, 1 = all except once centroid are rejected.
    opts.elemNumTh = 0;
end

centroid = clusters.centroid;
elemNum = clusters.elemNum;
minmaxX = clusters.minmaxX;
minmaxY = clusters.minmaxY;
minmaxZ = clusters.minmaxZ;

% Implementing distance condition
for i = 1:size(centroid,1)
    centroidi = centroid(i,:);
    if elemNum(i) ~= 0
        for j = 1:size(centroid,1)
            if j~=i
                centroidj = centroid(j,:);
                    if norm(centroidi-centroidj) <= opts.distTh 
                        centroid(i,:) = (centroidi.*elemNum(i) + centroidj.*elemNum(j))/...
                            (elemNum(i)+elemNum(j));
                        minmaxX(i,:) = [min(minmaxX(i,1),minmaxX(j,1)), max(minmaxX(i,2),minmaxX(j,2))];
                        minmaxY(i,:) = [min(minmaxY(i,1),minmaxY(j,1)), max(minmaxY(i,2),minmaxY(j,2))];
                        minmaxZ(i,:) = [min(minmaxZ(i,1),minmaxZ(j,1)), max(minmaxZ(i,2),minmaxZ(j,2))];
                        elemNum(i) = elemNum(i) + elemNum(j);
                        elemNum(j) = 0;
                    end
            end
        end
    end
end
idxNonZeroElem = elemNum ~=0;
centroid = centroid(idxNonZeroElem,:);
elemNum = elemNum(idxNonZeroElem,:);
minmaxX = minmaxX(idxNonZeroElem,:);
minmaxY = minmaxY(idxNonZeroElem,:);
minmaxZ = minmaxZ(idxNonZeroElem,:);
% Implementing XY distance condition
% Implementing distance condition
for i = 1:size(centroid,1)
    centroidi = centroid(i,1:2);
    if elemNum(i)~=0
        for j = 1:size(centroid,1)
            if j~=i
                centroidj = centroid(j,1:2);
                    if norm(centroidi-centroidj) <= opts.XYdistTh 
                        centroid(i,:) = (centroid(i,:).*elemNum(i) + centroid(j,:).*elemNum(j))/...
                            (elemNum(i)+elemNum(j));
                        minmaxX(i,:) = [min(minmaxX(i,1),minmaxX(j,1)), max(minmaxX(i,2),minmaxX(j,2))];
                        minmaxY(i,:) = [min(minmaxY(i,1),minmaxY(j,1)), max(minmaxY(i,2),minmaxY(j,2))];
                        minmaxZ(i,:) = [min(minmaxZ(i,1),minmaxZ(j,1)), max(minmaxZ(i,2),minmaxZ(j,2))];
                        elemNum(i) = elemNum(i) + elemNum(j);
                        elemNum(j) = 0;
                    end
            end
        end
    end
end

objHt = minmaxZ(:,2) - minmaxZ(:,1); % Heights of each cluster

% Implementing minimum height condition
for i = 1:size(centroid,1)
    if elemNum(i)~=0
        if objHt(i) < opts.minHeightRatio*max(objHt)
            % If height is less than minimum, reject it
            elemNum(i) = 0;
        end
    end
end

idxNonZeroElem = elemNum ~=0;
centroid = centroid(idxNonZeroElem,:);
elemNum = elemNum(idxNonZeroElem,:);
minmaxX = minmaxX(idxNonZeroElem,:);
minmaxY = minmaxY(idxNonZeroElem,:);
minmaxZ = minmaxZ(idxNonZeroElem,:);

% Implementing elemNum condition
maxElemNum = max(elemNum);
idxLessElem = find(elemNum >= (maxElemNum .* opts.elemNumTh));
centroid = centroid(idxLessElem,:);
elemNum = elemNum(idxLessElem);
minmaxX = minmaxX(idxLessElem,:);
minmaxY = minmaxY(idxLessElem,:);
minmaxZ = minmaxZ(idxLessElem,:);


clusterOut.centroid = centroid;
clusterOut.elemNum = elemNum;
clusterOut.minmaxX = minmaxX;
clusterOut.minmaxY = minmaxY;
clusterOut.minmaxZ = minmaxZ;
