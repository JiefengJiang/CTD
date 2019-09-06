%ROI-level analysis of CTD reinstatement. Even though it gives results for
%all HCP ROIs, only frontoparietal ones were tested.
function [coef, data] = MVRSABTContrastAnalysis(BTDataName, contrast)

BT = load(BTDataName);
intIdx = [25:32 41:48 57:64 73:80 89:96];
n_sub = length(BT.data);

if size(contrast, 1) == 1
    contrast = contrast';
end

coef = [];
data = cell([1, size(BT.data{1}, 1)]);
for i = 1 : length(data)
    data{i} = [];
end

for i = 1 : n_sub
    y = squeeze(mean(BT.data{i}(:, intIdx, :), 2));
    coef = [coef; (y * contrast)'];
    for j = 1 : length(data)
        data{j} = [data{j}; y(j, :)];
    end
end

