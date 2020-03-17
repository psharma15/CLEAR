% Generating A given xyzVoxelCoord, Tx-Rx placement and frequencies
% Pragya Sharma, ps847@cornell.edu, 12/18/2018
% -------------------------------------------------------------------------
function expConst = genA(roomSize,voxelSize,tagPosition,rxPosition,freq)

[xyzVoxelCoord,~,~,nVoxel] = genXYZ(roomSize,voxelSize);

nTag = size(tagPosition,1);
nRecv = size(rxPosition,1);
nFreq = length(freq);

c = 3e8;

% Finding distances between tag and each voxel, and receiver and each voxel
% distTagVoxel = single((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nTag 1]) ...
% - repmat(permute(tagPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');
distTagVoxel = ((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nTag 1]) ...
- repmat(permute(tagPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');

% distRxVoxel = single((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nRecv 1]) ...
% - repmat(permute(rxPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');
distRxVoxel = ((sqrt(sum(abs( repmat(permute(xyzVoxelCoord, [1 3 2]), [1 nRecv 1]) ...
- repmat(permute(rxPosition, [3 1 2]), [nVoxel(1)*nVoxel(2)*nVoxel(3) 1 1]) ).^2, 3)))');

sumDistTagRxVoxel = repmat(distTagVoxel,nRecv,1)+repelem(distRxVoxel,nTag,1);
% freqConst = single(repelem((-1j*2*pi*freq/c),nTag*nRecv,1));
freqConst = (repelem((-1j*2*pi*freq/c),nTag*nRecv,1));
sumDistTagRxVoxel = (repmat(sumDistTagRxVoxel,nFreq,1));

mulDistTagRxVoxel = repmat(distTagVoxel,nRecv,1).*repelem(distRxVoxel,nTag,1);
mulDistTagRxVoxel = repmat(mulDistTagRxVoxel,nFreq,1);
mulDistTagRxVoxelFreq = 2*pi*repelem(freq,nTag*nRecv,1).*(mulDistTagRxVoxel)/c;

clearvars distTagVoxel distRxVoxel

% expConst = single(exp(freqConst .* distTagRxVoxel));
% expConst = (exp(freqConst .* sumDistTagRxVoxel));

expConst = (mulDistTagRxVoxelFreq).*(exp(freqConst .* sumDistTagRxVoxel));

clearvars distTagRxVoxel freqConst

end