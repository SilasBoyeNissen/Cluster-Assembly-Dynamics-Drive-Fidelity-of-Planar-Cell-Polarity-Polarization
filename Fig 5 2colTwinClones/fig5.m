%% Figure 5: 2-color analyze
clear; rng(1); tic;
DX = 5;
BACK = 5;
DT = 5;

STAT = zeros(0, 5);
FILE = dir('1-TIF/c*.tif'); % Deposited on Stanford Digital Repository (https://doi.org/10.25740/sp776gx3909)
for n = 345:440
    F = FILE(n).name(1:end-4);
    disp(['n=' num2str(n) ': ' F]);
    load(['3-MAP/' F '.mat']);
    MAP = repelem(MAP, 2, 2, 1);
    CH1 = imgNEW2(F, 1, DT, BACK, DX);
    CH2 = imgNEW2(F, 2, DT, BACK, DX);
    for pd = unique(MAP(MAP>0))'
        ch1 = CH1;
        ch2 = CH2;
        ch1((blockproc(ch1, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP == pd) == 0) = 0;
        ch2((blockproc(ch2, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP == pd) == 0) = 0;
        ch1(~bwmorph(ch1, 'shrink')) = 0;
        ch2(~bwmorph(ch2, 'shrink')) = 0;
        [r1, c1, v1] = find(ch1);
        [r2, c2, v2] = find(ch2);
        C = intersect([knnsearch([r1 c1], [r2 c2]) (1:numel(r2))'], [(1:numel(r1))' knnsearch([r2 c2], [r1 c1])], 'rows');
        n1 = v1(setdiff(setdiff((1:numel(v1)), C(:, 1)), C(:, 1))');
        n2 = v2(setdiff(setdiff((1:numel(v2)), C(:, 2)), C(:, 2))');
        n12 = [v1(C(:, 1)) v2(C(:, 2)); n1 zeros(numel(n1), 1); zeros(numel(n2), 1) n2];
        STAT = [STAT; n*ones(size(n12, 1), 1) str2double(F(2:3))*ones(size(n12, 1), 1) pd*ones(size(n12, 1), 1) n12];
    end
    disp(['img' num2str(n) ' took ' num2str(round(toc)) ' sec']);
end
save('data.mat', 'STAT');

function I = imgNEW2(F, ch, DT, BACK, DX)
RAW = single(repelem(tiffreadVolume(['1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch 2 2*DT]}), 2, 2, 1));
if size(RAW, 3) == 1
    RAW = tiffreadVolume(['1-TIF/' F '.tif']);
    RAW = single(repelem(RAW(:, :, ch:2:2*DT), 2, 2, 1));
end
I = uint32(sum(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), 3));
end

%% Figure 5: 2-color visualize
figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [300 0 240 1100]);
tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'none');
fig2col(78, 1, 'Distal Vang', 'Proximal Fmi');
fig2col(76, 1, 'Fz', 'Proximal Vang');
fig2col(79, 1, 'Distal Vang', 'Proximal Vang');
fig2col(77, 1, 'Distal Vang', 'Proximal Pk');
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig5a', '-dpdf', '-r300');

figure(3); clf;
set(gcf, 'Color', 'w', 'Position', [540 0 240 1100]);
tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'none');
fig2col(78, 3, 'Proximal Vang', 'Distal Fmi');
fig2col(76, 3, 'Fz', 'Distal Vang');
fig2col(79, 3, 'Proximal Vang', 'Distal Vang');
fig2col(77, 3, 'Proximal Vang', 'Distal Pk');
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig5b', '-dpdf', '-r300');

function fig2col(v, pd, xlab, ylab)
nexttile;
load('data.mat');
n = numel(unique(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 1)));
n1 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 4))/10;
n2 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 5))/30;
n2(n1 == 0) = [];
n1(n1 == 0) = [];
n1(n2 == 0) = [];
n2(n2 == 0) = [];
scatter(n1, n2, 30, 'filled', 'MarkerFaceAlpha', 0.2);
[R, P, RL, RU] = corrcoef(n1, n2);
set(gca, 'Box', 'on', 'FontSize', 14, 'LineWidth', 1, 'XScale', 'log', 'YScale', 'log', 'XTick', [1 100 10000], 'YTick', [1 100 10000]);
title(['{\bf r=' num2str(R(2, 1), 1) '}, 95%CI [' num2str(RL(2, 1), 1) ', ' num2str(RU(2, 1), 1) ']' newline ...
    '{\bf p=' num2str(P(2, 1), 1) '}, n=' num2str(n) ', N=' num2str(numel(n1))], 'FontSize', 12, 'FontWeight', 'normal')
axis([1 10000 1 10000], 'square');
xlabel(xlab, 'FontWeight', 'bold', 'Color', '#EE220C', 'FontSize', 16);
ylabel(ylab, 'FontWeight', 'bold', 'Color', '#1DB100', 'FontSize', 16);
end