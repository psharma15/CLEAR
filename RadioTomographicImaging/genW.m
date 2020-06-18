% This function generates W in y=Wx, based on different models (elliptical
% is the first one) for linear RTI.

function W = genW(tagPosition,rxPosition,freq, roomSize,voxelSize,opts)

if ~isfield(opts,'model')
    opts.model = 'Ellipse';
end

nTag = length(tagPosition);
nRx = length(rxPosition);
nFreq = length(freq);

% Size of W is nData x [nVoxel_x*nVoxel_y*nVoxel_z]
nData = nTag*nRx*nFreq;
[xyzVoxelCoord,~,~,nVoxel] = genXYZ(roomSize,voxelSize); % Different folder
% Finding distances between tag and each voxel, and receiver and each voxel
% SIZE: [nTag, nVoxel]
distTagVoxel = ((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nTag 1]) ...
- repmat(permute(tagPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');

% SIZE: [nRx, nVoxel]
distRxVoxel = ((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nRx 1]) ...
- repmat(permute(rxPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');

distTagRx = ((sqrt(sum(abs( repmat(permute(rxPosition, [1 3 2]), [1 nTag 1]) ...
- repmat(permute(tagPosition, [3 1 2]), [nRx 1 1]) ).^2, 3)))');

W = zeros(nData,nVoxel(1)*nVoxel(2)*nVoxel(3));

switch(opts.model)
    case 'Ellipse'
        % Within this ellipse, weight is 1, otherwise 0. It could be
        % converted to Gaussian or something.
        if ~isfield(opts,'lambda')
            % Width parameter of ellipsoid, generally set low, such that 
            % essentially LoS model
            opts.lambda = 0.15; % In meters
        end
        % It is goind to compare each column of lhs with rhs. Tested.
        % W is a logical matrix in first step
        W = (repmat(distTagVoxel,nRx,1)+repelem(distRxVoxel,nTag,1)) < ...
            (distTagRx(:) + opts.lambda);
        W = (W./sqrt(distTagRx(:))); % This is in Patwari's paper, not Xu.
                
    otherwise
        fprintf('Enter correct genW model.\n');
end
   


end