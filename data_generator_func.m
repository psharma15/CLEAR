function [C,rssiStdDev,phStdDev] = data_generator_func(data, opts)

% System parameters

ChaNum = opts.nFreq;

TagNum = opts.nTag;

AntNum = opts.nRx;

%% Extract the data.
% The six columns of the data are (from left to right): channel index, tag
% index, antenna ID, tag RSSI/dBm, tag RSSI/nW, and tag phase.

% [~, text, raw] = xlsread(filename, 1, 'A2:F4000');
% data = str2double(raw);


%% Data extraction and grouping for each channel
C = zeros(TagNum, AntNum, ChaNum);
rssiStdDev = zeros(TagNum, AntNum, ChaNum);
phStdDev = zeros(TagNum, AntNum, ChaNum);

for i = 1:ChaNum
    [~, ~, C(:, :, i),rssiStdDev(:,:,i),phStdDev(:,:,i)] = data_mat(data, i, TagNum, AntNum, opts);
end

% C1 = C(1, :, :);
% C2 = C(2, :, :);
% C3 = C(3, :, :);
% C4 = C(4, :, :);
% C5 = C(5, :, :);
% C6 = C(6, :, :);
% 
end
