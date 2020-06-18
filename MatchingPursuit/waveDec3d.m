% 3D wavelet decomposition 
codePathMP = 'D:\Research\SummerFall17Spring18\CnC\ArpaeProject\Algorithms\MP';
addpath(codePathMP);

opts = [];
% location of Tx-Rx antenna and frequencies.
if ~isfield(opts,'imSz')
    opts.imSz = [0.1, 0.8; 0.1, 0.8; 0, 0.3]; % meters
end
if ~isfield(opts,'voxSz')
    opts.voxSz = [0.01;0.01;0.01]; % meters
end
if ~isfield(opts,'freq')
    opts.freq = [1.7006;1.7591;1.8185;1.8779;1.9364;1.9958;2.0543;2.1137;2.1722;2.2316;2.2901]*1e9;
end
if ~isfield(opts,'posRxTx')
    opts.posRxTx = 1;
end

%% Generating image
[xyzVoxCoord,voxImCoord,xyzSlice,nVoxel] = genXYZ(opts.imSz,opts.voxSz);
obj1.Sz = [0.08, 0.06, 0.3];
objSzReorder = sort(obj1.Sz);
obj1.Sz = [objSzReorder(2), objSzReorder(1), objSzReorder(3)];
objCenterGrid = arrGrid1(opts.imSz,obj1.Sz);

caseNum = 10;
obj1.Center = objCenterGrid(caseNum,:);
obj1.Orient = 0;

imgX = genImgEllip(nVoxel,xyzVoxCoord,voxImCoord,obj1.Sz,obj1.Center,obj1.Orient,xyzSlice,opts.imSz);
p = [2 1 3];
imgX = reshape(imgX,nVoxel); % Reshaping image to 3D
imgX = permute(imgX,p);

%% Performing 1D Wavelet Transform

% % Taking a vector of the image along z axis
% a = imgX(42,65,:);
% a = a(:);
% 
% dydLen = nextpow2(length(a));
% a = [a;zeros(2^dydLen - length(a),1)]; % Padding with zeros
% 
% nMax = fix(log2(length(a))); % Number of level should be less than equal to this.
% wName = 'db4';
% wavLevel = 1;
% [c,l] = wavedec(a,wavLevel,wName); % Daubechies D4 (db2) Wavelet
% 
% figure
% subplot(3,1,1)
% plot(a); title('Original vector')
% 
% aNew = waverec(c,l,wName); % Near perfect reconstruction
% 
% subplot(3,1,2)
% plot(aNew); 
% title(['After perfect reconstruction at decomposition level ',num2str(wavLevel)])
% 
% fprintf('Error in perfect wavelet reconstruction: %f\n',norm(a-aNew));
% 
% % Approximation
% detCoeffLevel = 1:wavLevel;
% percentRej = 0.05;
% NC = wthcoef('t',c,l,detCoeffLevel,percentRej*ones(1,wavLevel),'s');
% aNew2 = waverec(NC,l,wName); 
% subplot(3,1,3)
% plot(aNew2); 
% title(['Reconstruction after thresholding ',num2str(percentRej),'% detail coefficients at each level'])
% 
% fprintf('Error in lossy wavelet reconstruction: %f\n',norm(a-aNew2));

%% Performing 2D Wavelet Transform

% Taking a YZ slice of the image.
a = imgX(:,66,:);
a = permute(a,[1,3,2]);

% Pad the image
[m,n] = size(a);
a = [a;zeros((2^nextpow2(m))-m,n)];
a = [a,zeros(2^nextpow2(m),(2^nextpow2(n))-n)];


figure
imagesc(a); title('Original Image after zero-padding')

wavLevel = 1;
wName = 'db4';
L = wmaxlev([m,n],wName);
fprintf('Maximum possible level of decomposition: %d\n',L);


[c,s] = wavedec2(a,wavLevel,wName);
[H1,V1,D1] = detcoef2('all',c,s,1);
A1 = appcoef2(c,s,wName,1);

figure
subplot(2,2,1)
imagesc(A1);
title('Approximation Coef. of Level 1');

subplot(2,2,2);
imagesc(H1);
title('Horizontal detail Coef. of Level 1');

subplot(2,2,3);
imagesc(V1);
title('Vertical detail Coef. of Level 1');

subplot(2,2,4);
imagesc(D1);
title('Diagonal detail Coef. of Level 1');

%%