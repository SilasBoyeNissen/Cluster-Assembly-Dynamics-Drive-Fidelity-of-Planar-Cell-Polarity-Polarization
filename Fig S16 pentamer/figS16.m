% Figure S16: Pentamer
clear; rng(1); tic;
DX = 3;
DT = 4;
BACK = 2;
BINW = 1;
FMIN = 80;
FMAX = 160;

RES = [];
FILE = dir('data/d*.tif');
for n = size(FILE, 1)-2:size(FILE, 1)
    F = FILE(n).name(1:end-4); disp(F);
    RAW = single(tiffreadVolume(['data/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [1 1 inf]}));
    if size(RAW, 3) == 1
        RAW = tiffreadVolume(['data/' F '.tif']);
        RAW = single(RAW(:, :, 1:1:end));
    end
    SMOO = movmedian(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), [DT 0], 3);
    MASK = (blockproc(SMOO(:, :, 1), [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2]));
    [r3, c3] = find(MASK); % puncta centers
    smoo = reshape(permute(SMOO, [3 1 2]), size(RAW, 3), []);
    smor = smoo(:, (MASK)); % puncta centers
    SMOO = SMOO(:, :, 1);
    res = nan(size(r3, 1), 9);
    tank = logical(triu(ones(size(smoo, 1), size(smoo, 1)), 1));
    parfor i = 1:size(smor, 2)
        wi = smor(:, i) - smor(:, i)';
        [pwr, f] = pspectrum(histcounts(wi(tank), 'BinWidth', BINW, 'Normalization', 'probability'), 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/FMIN]);
        I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
        if ~isempty(I)
            res(i, :) = [n 0 0 BINW*2*pi/f(I) smor(1, i)/(BINW*2*pi/f(I)) c3(i) r3(i) 4 1];
        end
    end
    res(isnan(res(:, 4)), :) = [];
    RES = [RES; res];
end

figure(1); clf;
set(gcf, 'color', 'w')
histogram(RES(:, 5), -0.5:1:15.5)
title(['N = ' num2str(size(RES, 1)) '; \mu = ' num2str(round(mean(RES(:, 5)), 1)) ' +- ' num2str(round(std(RES(:, 5)), 1))]);
xlabel('Cluster size', 'FontWeight', 'bold');
ylabel('Number of clusters', 'FontWeight', 'bold');
set(gca, 'LineWidth', 2, 'FontSize', 20);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('figS16', '-dpdf', '-r300');
