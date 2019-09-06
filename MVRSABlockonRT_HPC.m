%For hippocampal only
%Modulation of hippocampal univariate activity and pattern separation
%measures on RT
function [coef, p, q] = MVRSABlockonRT_HPC()

matName = 'HPC_BB_Trial_RSA_run2-6.mat';

% load BT as a convaiate
BT = load('HCP_BT_Trial_RSA_run2-6.mat');
intIdx_BT = [25:32 41:48 57:64 73:80 89:96];

ids = [1 3:20 22:24 26:32 34:36 38];
load(matName);
nSub = length(ids);
nROI = 1;

coef = zeros(nSub, 5);
q = zeros(nSub, 20); 
p = 1;

intIdx = [1:8 9:16 17:24 25:32 33:40];

idx1 = zeros(4, 12);

for i = 1:length(ids)
    data2 = data{i}(intIdx, :);
    uni_data = univar{i};
    cur_BT = squeeze(BT.data{i}(266, intIdx_BT, :));
    cur_BT = cur_BT * [1 -1 1 -1 0 0]';    
    
    %this encodes information of which trials were excluded from fMRI
    %analysis
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
        cur_int = data2;
        
        %run 2-6 only
        for k = 9:48
            trialIdx = (roomToTime(k) - 1) * 8 + 1: roomToTime(k) * 8;
            for i0 = 1 : 3
                m(trialIdx, i0) = cur_int(k - 8, i0);
            end
            m(trialIdx, 4) = uni_data(k - 8);
            m(trialIdx, 5) = cur_BT(k - 8);
        end
        
        %align the direction of reinstatement
        m(unexpected > 0.5, :) = -m(unexpected > 0.5, :);

        m1 = m(runs.idx, :);
        
        %get quintile numbers
        q(i, 1:5) = quantileData(m1(x1(runs.idx, 5) > 0, 4), y(x1(runs.idx, 5) > 0), 5);
        q(i, 6:10) = quantileData(m1(x1(runs.idx, 6) > 0, 4), y(x1(runs.idx, 6) > 0), 5);
        q(i, 11:15) = quantileData(m(runs.idx, 3), y, 5);
        q(i, 16:20) = quantileData(m(runs.idx, 4), y, 5);

        %with CTD reinstatement in dlPFC as covariate
        x = m;
        %without CTD reinstatement in dlPFC as covariate
        x = m(:, 1:4);
        x = x(runs.idx, :);
        x = (x - repmat(mean(x), [size(x, 1) 1])) ./ repmat(std(x), [size(x, 1) 1]);
        x = [x  x1(runs.idx, 5:6) x4(runs.idx)];
        b = pinv(x) * y;
        coef(i, 1:4) = b(1:4);
    end
end

[~, p, ~, stats] = ttest(coef);

    
    