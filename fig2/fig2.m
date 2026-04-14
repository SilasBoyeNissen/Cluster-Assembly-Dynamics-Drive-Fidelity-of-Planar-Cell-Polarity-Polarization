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
RES(RES(:, 8) == 5, 8) = 4; %%
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

% Figure 2: DIX (VANG)
pd = 4;
figure(1); clf;
YT = [0:10:60; 0:5:30; 0:1:6];
set(gcf, 'Color', 'w', 'Position', [10 10 1200 460]);
tl = tiledlayout(1, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 18);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 18);
k = 0;
vv = [61 64 62];
for v = vv
    k = k + 1;
    nt = nexttile;
    if v > 61
        tbl = table(RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 3), RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 2), RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 5),'VariableNames',{'GeneticID','Age','ClusterSize'});
        lme = fitlme(tbl,'ClusterSize~Age+GeneticID+(Age|GeneticID)')
    end
    errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', [CO4(k, :) 1], 'LineWidth', 2, 'Marker', 'none'); hold on;
    if v > 61
        text(14, 16, ['P'' = ' num2str(round(double(lme.Coefficients(2, 6)), 1, 'significant'))], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        errorbar(14, 15, 3, 'horizontal', 'Color', 'k', 'LineWidth', 2);
        ax = gca;
        ax.Clipping = 'off';
    end
    mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
    N1 = sum(N(N(:, 3) == v, 7));
    [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
    plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', [CO4(k, :) 1], 'LineWidth', 2);
    fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO4(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    xlabel('');
    ylabel('');
    if v ~= 61
        set(gca, 'Box', 'on', 'FontSize', 18, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
    else
        set(gca, 'Box', 'on', 'FontSize', 18, 'LineWidth', 1, 'XTick', XT, 'YTick', '');
    end
    axis([14 33 0 24.99], 'square');
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
        '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
        '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
        '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
        'FontSize', 15, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
set(gca, 'XTick', XT);
set(nexttile(1), 'YTick', YT(2, :), 'XTick', XT);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig2', '-dpdf', '-r300');