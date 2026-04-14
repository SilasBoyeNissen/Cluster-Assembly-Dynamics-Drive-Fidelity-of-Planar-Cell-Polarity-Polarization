% LOAD & Prepare
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

% Figure S5: Density
tic;
for i = 1:max(RES(:, 1))
    RESs = RES(RES(:, 1) == i, :);
    [~, D] = knnsearch(RESs(:, 6:7), RESs(:, 6:7), 'K', 30);
    RES(RES(:, 1) == i, 10) = sum(D*0.065 < sqrt(1/pi), 2);
end

pd = 4;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 20);
ylabel(tl, 'Density [clusters/\mum^2]', 'Fontweight', 'bold', 'FontSize', 20);
for v = [11 12 14 13 15 16]
    nt = nexttile;
    RESs = RES(RES(:, 3) == v & RES(:, 8) > pd, :);
    j = 1;
    sta = zeros(numel(unique(RESs(:, 1))), 3);
    for i = unique(RESs(:, 1))'
        errorbar(unique(RESs(RESs(:, 1) == i, 2))+(rand()/3-1/6), mean(RESs(RESs(:, 1) == i, 10)), std(RESs(RESs(:, 1) == i, 10)), 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 2); hold on;
        sta(j, :) = [unique(RESs(RESs(:, 1) == i, 2)) mean(RESs(RESs(:, 1) == i, 10)) numel(RESs(RESs(:, 1) == i, 10))];
        j = j + 1;
    end
    mdl = fitlm(sta(:, 1), sta(:, 2), 'Weights', sta(:, 3));
    [Ypred, Yci] = predict(mdl, sta(:, 1));
    plot(sta(:, 1), Ypred, '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
    fill([sta(:, 1)', fliplr(sta(:, 1)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    xlabel('');
    ylabel('');
    set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1);
    axis([14 33 0 9]);
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, ['{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
        '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
        '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
        'FontSize', 16, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
set(nexttile(1), 'YTick', 0:2:8, 'XTick', '');
set(nexttile(2), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(3), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(4), 'YTick', 0:2:8);
set(nexttile(5), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(6), 'YTick', 0:2:8, 'YTickLabel', '');
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('figS5', '-dpdf', '-r300');