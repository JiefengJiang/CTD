%Testing modulation of hippocampal activity (pattern) on CTD reinstatement
%in HCP ROIs, here we tested right dlPFC/frontopolar (266)
%first 5 regressors are in the order shown in Fig. 5
function [coef, yy, q] = MVRSACorrelationAnalysis(BTDataName, BBDataName, BT_roi, contrast)

BT = load(BTDataName);
BB = load(BBDataName);
intIdx = [25:32 41:48 57:64 73:80 89:96];
intIdx1 = [25:32 41:48 57:64 73:80 89:96] - 16;
n_sub = length(BB.data);

if size(contrast, 1) == 1
    contrast = contrast';
end

coef = [];
yy = [];
q = [];

run_regressors = zeros(40, 5);
for i = 1 : 5
    run_regressors((i - 1) * 8 + 1 : i * 8, i) = 1;
end

item_regressors = zeros(40, 4);
for i = 1 : 4
    item_regressors((i - 1) * 2 + 1 : 8 : 40, i) = 1;
    item_regressors(i * 2 : 8 : 40, i) = 1;
end

for i = 1 : n_sub
    x1 = [BB.data{i}];
    x1 = (x1 - repmat(mean(x1), [40 1])) ./ repmat(std(x1), [40 1]);
    y = squeeze(BT.data{i}(BT_roi, intIdx, :));
    yy = [yy; mean(y)];
    y = y * contrast;
    y1 = y;
    y = (y - mean(y)) / std(y);
    x2 = [BT.univar{i}(intIdx1, BT_roi) BB.univar{i}];
    x2 = (x2 - repmat(mean(x2), [40 1])) ./ repmat(std(x2), [40 1]);
    meanX2 = mean(x2);
    stdX2 = std(x2);
    
    %exclude outliers
    uni_idx = x2(:, 1) >= -5 & x2(:, 1) <= 5 & x2(:, 2) >= -5 & x2(:, 2) <= 5;
    
    x = [x1 x2 (1:40)' ones(size(x1, 1), 1)];
    x = x(uni_idx, :);
    x = [x ones(size(x, 1), 1)];
    y = y(uni_idx);
    q = [q; quantileData(x1(uni_idx, 2), y1(uni_idx), 5)];
    b = pinv(x) * y;
    coef = [coef; b'];
end

