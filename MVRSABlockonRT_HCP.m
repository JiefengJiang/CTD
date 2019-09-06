% Modulation of block-onset reinstatement of CTD on subsequent trial RT
% This script can test the modulations at navigation cue onset (first row 
% of p) and block-onset (second row of p), but we only tested the 
% block-onset results for the 4 ROIs showing significant reinstatement of
% CTD
function [coef, p, q] = MVRSABlockonRT_HCP(matName, contrast)

if size(contrast, 1) == 1
    contrast = contrast';
end

ids = [1 3:20 22:24 26:32 34:36 38];
load(matName);
nSub = length(ids);
nROI = size(data{1}, 1);

coef = zeros(nSub, nROI, 2);
q = zeros(nSub, nROI, 2, 10); 
p = ones(2, nROI);

extIdx = [1:8 17:24 33:40 49:56 65:72 81:88];
intIdx = [9:16 25:32 41:48 57:64 73:80 89:96];

extIdx1 = [1:8 1:8 17:24 33:40 49:56 65:72];
intIdx1 = [9:16 9:16 25:32 41:48 57:64 73:80];

idx1 = zeros(4, 12);

for i = 1:length(ids)
    data1 = data{i}(:, extIdx, :);
    data2 = data{i}(:, intIdx, :);
    data3 = univar{i}(extIdx1, :);
    data4 = univar{i}(intIdx1, :);
    
    load(sprintf('trialIdx_%02d.mat', ids(i)));
    
    %temporal PE
    modelName = sprintf('%03d.mat', ids(i));
    load(modelName);
    x4 = [res.T.models.pe_t' res.T.models.pe_t'];
    x4(runs.task == 1, 1) = 0;
    x4(runs.task == 0, 2) = 0;
    x4 = res.T.models.pe_t';
    
    %mask out run 1
    runs.idx(1:64) = 0;
    
    unexpected = zeros(length(runs.idx), 1);
    idx2 = (runs.room > 1) & (runs.task == 0);
    unexpected(idx2) = 1;
    idx2 = (runs.room < 2) & (runs.task == 1);
    unexpected(idx2) = 1;
    
    y = log(runs.RT');
    y = y(runs.idx);
    nTrials = length(runs.idx);
    
    x1 = zeros(nTrials, 6);
    for j = 0 : 3
        x1(runs.room == j, j + 1) = 1;
    end
    x1(runs.task == 0, 5) = 1;
    x1(runs.task == 1, 6) = 1;
    
    %reconstruct room order
    rooms = reshape(runs.room, [8 48]);
    rooms = rooms(1, :);
    roomToTime = zeros(1, 48);
    for j = 0 : 3
        idx1(j + 1, :) = find(rooms == j);
    end
    
    for j = 1 : 6
        x3 = idx1(:, (j * 2 - 1 : j * 2))';
        roomToTime((j - 1) * 8 + 1 : j * 8) = reshape(x3, [1 8]);
    end
    
    m = zeros(length(runs.idx), 4);
    
    for j = 1 : nROI
        cur_ext = squeeze(data1(j, :, :)) * contrast;
        cur_int = squeeze(data2(j, :, :)) * contrast;
        
        cur_ext1 = data3(:, j);
        cur_int1 = data4(:, j);
        
        %run 2-6 only
        for k = 9:48
            trialIdx = (roomToTime(k) - 1) * 8 + 1 : roomToTime(k) * 8;
            m(trialIdx, 1) = cur_ext(k, 1);
            m(trialIdx, 2) = cur_int(k, 1);
            m(trialIdx, 3) = cur_ext1(k, 1);
            m(trialIdx, 4) = cur_int1(k, 1);
        end
        
        m(unexpected > 0.5, :) = -m(unexpected > 0.5, :);
        
        for k = 1 : 2
            x2 = zeros(nTrials, 4);
            x2(runs.task == 0, 1:2) = m(runs.task == 0, [k k + 2]);
            x2(runs.task == 1, 3:4) = m(runs.task == 1, [k k + 2]);
            
            x = [x2 x1(:, 5:6) x4];
            x = x(runs.idx, :);
            
            %get quintines
            q(i, j, k, 1:5) = quantileData(x(x(:, 1) ~= 0, 1), y(x(:, 1)~=0), 5);
            q(i, j, k, 6:10) = quantileData(x(x(:, 2) ~= 0, 2), y(x(:, 2)~=0), 5);
            

            x = [x2(:, 1:2) + x2(:, 3:4)];
            x = x(runs.idx, :);
            x = (x - repmat(mean(x), [size(x, 1) 1])) ./ repmat(std(x), [size(x, 1) 1]);
            x = [x x1(runs.idx, 5:6) x4(runs.idx, :)];
            
            b = pinv(x) * y;
            coef(i, j, k) = b(1);
        end
    end
end

for i = 1 : 2
    xx = squeeze(coef(:, :, i));
    [~, p(i, :), ~, stats] = ttest(xx);
end
    
    
    