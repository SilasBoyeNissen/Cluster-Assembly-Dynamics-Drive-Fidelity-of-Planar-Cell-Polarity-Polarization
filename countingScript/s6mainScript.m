%% MAIN script. Script #6 counts molecules present at all locations within the map automatically or manually.
close all; clear; rng(1); tic;
DX = 3;
DT = 20;
BACK = 2;
BINW = 0.3;
FMAX = 110;

RES = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 1:size(FILE, 1)
    F = FILE(n).name(1:end-4); disp(F);
    CH = str2double(F(22));
    load(['../../3-MAP/' F '.mat']);
    DIM = [find(sum(MAP), 1) find(sum(MAP, 2), 1) find(sum(MAP), 1, 'last') find(sum(MAP, 2), 1, 'last')];
    MAP = MAP(DIM(2):DIM(4), DIM(1):DIM(3));
    for ch = 1:CH
        RAW = single(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[DIM(2) 1 DIM(4)], [DIM(1) 1 DIM(3)], [ch CH inf]}));
        if size(RAW, 3) == 1
            RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
            RAW = single(RAW(DIM(2):DIM(4), DIM(1):DIM(3), ch:CH:end));
        end
        SMOO = movmedian(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), [DT 0], 3);
        MASK = (blockproc(SMOO(:, :, 1), [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2]) & MAP);
        [r3, c3] = find(MASK); % puncta centers
        smoo = reshape(permute(SMOO, [3 1 2]), size(RAW, 3), []);
        smor = smoo(:, (MASK)); % puncta centers
        SMOO = SMOO(:, :, 1);
        res = nan(size(r3, 1), 9);
        tank = logical(triu(ones(size(smoo, 1), size(smoo, 1)), 1));
        for i = 1:size(smor, 2) % parfor
            wi = smor(:, i) - smor(:, i)';
            [pwr, f] = pspectrum(histcounts(wi(tank), 'BinWidth', BINW), 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/40]);
            I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
            if ~isempty(I)
                res(i, :) = [n str2double(F(17:18)) str2double(F(2:3)) BINW*2*pi/f(I) smor(1, i)/(BINW*2*pi/f(I)) c3(i) r3(i) MAP(r3(i), c3(i)) ch];
            else
                [pddf, edges] = histcounts(wi(tank), 'BinWidth', BINW, 'Normalization', 'probability');
                [pwr, f] = pspectrum(pddf, 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/40]);
                I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
                stepSize(F, i, smor, edges, pddf, f, pwr, BINW*2*pi/f(I), ch);
                figure(2);
                figure(1);
                w = waitforbuttonpress;
                if w
                    cs = get(gcf, 'CurrentCharacter');
                    disp([ num2str(n) ': ' num2str(i) ' out of ' num2str(size(smor, 2))])
                end
                res(i, :) = [n str2double(F(17:18)) str2double(F(2:3)) 0 str2double(cs)-1 c3(i) r3(i) MAP(r3(i), c3(i)) ch];
            end
        end
        save(['../../7-res/' F '.mat'], 'res');
        res(isnan(res(:, 4)), :) = [];
        RES = [RES; res];
    end
    disp(['img' num2str(n) ' took ' num2str(round(toc)) ' sec']);
end

%% SAVE step size, used for manual visualization of the bleaching traces and identifying the fluorophore step sizes.
function stepSize(F, i, smor, edges, pddf, f, pwr, ss, ch)
figure(1); clf;
tiledlayout(1, 3);
set(gcf, 'color', 'w', 'Position', [0 0 1800 600]);

nexttile;
%plot(squeeze(valr(:, i)), 'color', [0 0 1 0.15], 'LineWidth', 0.1); hold on;
plot([0 size(smor, 1)], [0 0], 'color', [0 0 0 0.5], 'LineWidth', 3); hold on;
plot(squeeze(smor(:, i)), 'color', [0, 0.4470, 0.7410 1], 'LineWidth', 4);
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:20:1000, 'YTick', 0:100:500);
xlabel('Frames', 'Fontweight', 'bold');
ylabel('Intensity', 'Fontweight', 'bold');
axis([0 100 -200 500]);

nexttile;
plot(edges(1:end-1), pddf, 'LineWidth', 6);
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:200:800, 'YTick', [1e-6 1e-5 1e-4 1e-3 1e-2], 'YScale', 'log');
xlabel('Pairwise difference', 'Fontweight', 'bold');
ylabel('Probability', 'Fontweight', 'bold');
title(['Channel ' num2str(ch)])
axis([-200 800 1e-6 0.0101]);

nexttile;
INTS = 2*(edges(2)-edges(1))*pi./f;
plot(INTS, pwr, 'LineWidth', 6); hold on;
scatter(round(ss), max(pwr(round(INTS) == round(ss))), 400, 'r', 'filled');
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:40:200);
xlabel('Intensity', 'Fontweight', 'bold');
ylabel('Power', 'Fontweight', 'bold');
axis([0 160 0 1.5e-7], 'auto y');

mkdir(['../../6-StepSize/' F '-' num2str(ch) '/']);
warning('off', 'MATLAB:MKDIR:DirectoryExists');
saveas(gcf, ['../../6-StepSize/' F '-' num2str(ch) '/' num2str(i) '.png']);
end