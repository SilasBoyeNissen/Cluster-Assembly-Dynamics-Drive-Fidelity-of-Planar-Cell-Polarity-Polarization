%% LOAD & Prepare
close all; clear;
RES = [];
FILE = dir('7-res/c*');
for n = 1:size(FILE, 1)
    load(['7-res/' FILE(n).name]);
    RES = [RES; res];
end
RES(:, 5) = floor(RES(:, 5));

PROT = {'Dgo', 'Dsh', 'Fmi', 'Fz', 'Pk', 'Vang'};
PD = {['{\color[rgb]{' num2str([1 113 1]/255) '}Distal}'], '', ['{\color[rgb]{' num2str([108 83 142]/255) '}Proximal}']};
MUT = {'', '', 'fz^{null}', 'vang^{null}', 'fz^{null}', 'fz^{null}', '', 'pk^{null}',};
CO4 = [0 0 0; 8 29 88; 34 94 168; 65 182 196]/255; % Colorbrewer (9 colors, sequential, color-blind)
CO6 = [90 174 97; 27 120 55; 186 135 45; 0 68 27; 153 112 171; 64 0 75]/255;
FMI = 1;
YMI = 0.01;
x = 0:1:200;
RES(RES(:, 5) < 0, :) = [];
%RES(RES(:, 8) == 5, 8) = 4; %%
[~, ia] = unique(RES(:, 1));
N = RES(ia, 1:3);
Nl = RES(ia, 1:3); % lambda
Nle = RES(ia, 1:3); % lambda error
FILE = dir('7-res/c*.mat');
for i = 1:size(FILE, 1)
    F = FILE(i).name(1:end-4); 
    for ch = 1:2
        for pd = 1:6
            temp = RES(RES(:, 1) == i & RES(:, 8) == pd & RES(:, 9) == ch, 6:7);
            [~, d] = knnsearch(temp, temp, 'K', 2);
            if size(d, 2) > 1
                RES(RES(:, 1) == i & RES(:, 8) == pd & RES(:, 9) == ch, 10) = d(:, 2)*65/1000;
            end
            N(N(:, 1) == i, (ch-1)*6+pd+3) = sum(RES(:, 1) == i & RES(:, 8) == pd & RES(:, 9) == ch);
            y1 = 1-histcounts(RES(RES(:, 1) == i & RES(:, 8) == pd & RES(:, 9) == ch, 5), x, 'Normalization', 'cdf');
            t1 = find(y1==mink(unique(y1), 1), 1); %
            if sum(RES(:, 1) == i & RES(:, 8) == pd & RES(:, 9) == ch) > 3
                MO1 = fitnlm(x(FMI:t1), y1(FMI:t1), @(b,x) b(1).*exp(-x./b(2)), [1 1]);
                Nl(Nl(:, 1) == i, (ch-1)*6+pd+3) = table2array(MO1.Coefficients(2,1));
                Nle(Nle(:, 1) == i, (ch-1)*6+pd+3) = table2array(MO1.Coefficients(2,2));
            end
        end
    end
end
N(isnan(N)) = 0;
Nl(Nl == 0) = NaN;
Nle(Nle == 0) = NaN;
XT = 15:5:30;

%% Figure S6: Mobility analyze
rng(1);
DX = 3;
BACK = 2;
DT = 20;

STAT = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 1:size(FILE, 1)
    F = FILE(n).name(1:end-4);
    if str2double(F(2)) == 2
        tic;
        ch = 1;
        CH = str2double(F(22));
        disp(['n=' num2str(n) ': ' F]);
        load(['../../3-MAP/' F '.mat']);
        DIM = [find(sum(MAP), 1) find(sum(MAP, 2), 1) find(sum(MAP), 1, 'last') find(sum(MAP, 2), 1, 'last')];
        RAW = single(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[DIM(2) 1 DIM(4)], [DIM(1) 1 DIM(3)], [ch CH inf]}));
        if size(RAW, 3) == 1
            RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
            RAW = single(RAW(DIM(2):DIM(4), DIM(1):DIM(3), ch:CH:end));
        end
        MAP = MAP(DIM(2):DIM(4), DIM(1):DIM(3));
        SMOO = movmedian(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), [DT 0], 3);
        MASK = false(size(SMOO));
        parfor i = 1:size(SMOO, 3)
            MASK(:, :, i) = (blockproc(SMOO(:, :, i), [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP);
        end
        CC = bwconncomp(MASK);
        L = labelmatrix(CC);
        RESs = RES(RES(:, 1) == n, :);
        for i = 1:size(RESs, 1)
            [x, y, z] = ind2sub(size(MASK), CC.PixelIdxList{L(RESs(i, 7), RESs(i, 6), 1)});
            RESs(i, 9:14) = [y(end) x(end) z(end) sqrt((x(end)-RESs(i, 7))^2 + (y(end)-RESs(i, 6))^2)/z(end) sum(sqrt((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2)) mean((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2)/4];
        end
        STAT = [STAT; RESs];
        toc;
    end
end
save('data.mat', 'STAT');

%% Figure S6: Mobility visualize
figure(1); clf;
load('data');
TH = 2;
set(gcf, 'Color', 'w', 'Position', [0 0 900 900]);
tl = tiledlayout(4, 2, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 15);
ylabel(tl, 'Log diffusion coefficient [\mum^2/s]', 'Fontweight', 'bold', 'FontSize', 15);
for v = 23:26
    tak = [];
    for pd = [3 1]
        nt = nexttile;
        nn = N(N(:, 3) == v, :);
        x = STAT(STAT(:, 3) == v & STAT(:, 11) > TH & STAT(:, 8) == pd, 5);
        y = STAT(STAT(:, 3) == v & STAT(:, 11) > TH & STAT(:, 8) == pd, 14)./(STAT(STAT(:, 3) == v & STAT(:, 11) > TH & STAT(:, 8) == pd, 11)-1)*(0.065)^2/0.050;
        sta = zeros(39, 3);
        x(y == 0) = [];
        y(y == 0) = [];
        y = log(y);
        tak = [tak; pd*ones(size(x)) x y];
        if pd == 1
            tbl = table(tak(:, 1), tak(:, 2), tak(:, 3), 'VariableNames', {'Side', 'ClusterSize', 'DiffCoef'});
            lme = fitlme(tbl,'DiffCoef~ClusterSize+Side+(ClusterSize|Side)');
        end
        for i = 1:39
            sta(i, 1:3) = [i mean(y(x == i)) std(y(x == i))];
        end
        sta(isnan(sta(:, 2)), :) = [];
        y(x==0) = [];
        x(x==0) = [];
        scatter(x, y, 15, [CO6(mod(v, 10), :)], 'filled', 'MarkerFaceAlpha', 0.3); hold on;
        plot(sta(:, 1), sta(:, 2), '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
        fill([sta(:, 1)', fliplr(sta(:, 1)')], [(sta(:, 2)+sta(:, 3))' fliplr((sta(:, 2)-sta(:, 3))')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xlabel('');
        ylabel('');
        set(gca, 'Box', 'on', 'FontSize', 15, 'LineWidth', 1, 'XTick', 0:5:20, 'YTick', -13:2:-5, 'YTickLabel', '');
        axis([0 25 -14 -3]);
        text(max(xlim)-diff(xlim)/30, max(ylim)-diff(ylim)/20, ['n=' num2str(size(nn, 1)) ', N=' num2str(size(x, 1))], 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        if mod(tilenum(nt), 2) == 0
            text(max(xlim)-diff(xlim)/30, max(ylim)-diff(ylim)/6, ['P'' = ' mat2str(round(double(lme.Coefficients(2, 6)), 1, 'significant'))], 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        end
    end
end
set(gca, 'XTick', 0:5:30);
set(nexttile(1), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(3), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(5), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(7), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5});
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('figS6', '-dpdf', '-r300');