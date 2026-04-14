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

% Figure 1: P or D single exponentials
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 20);
ylabel(tl, '1 - Cumulative Distribution Function', 'Fontweight', 'bold', 'FontSize', 20);
for pd = [1 3]
    for v = 23:26
        TXT = {};
        nt = nexttile;
        for apf = 1:2
            if apf == 1
                nn = N(N(:, 3) == v & N(:, 2) < 24, :);
                txt = ' 15-23 hr APF';
                ft1 = '--';
                ft2 = ':';
                ti = 1;
            else
                nn = N(N(:, 3) == v & N(:, 2) > 23, :);
                txt = ' 24-32 hr APF';
                ft1 = '-';
                ft2 = '-';
                ti = 0.5;
            end
            stat = nan(size(nn, 1), 2);
            j = 1;
            for i = nn(:, 1)'
                y1 = 1-histcounts(RES(RES(:, 1) == i & RES(:, 8) == pd, 5), x, 'Normalization', 'cdf');
                plot(x(1:end-1), y1, ft1, 'Color', [CO6(v-20, :) 1], 'LineWidth', ti); hold on;
                TXT{end+1} = '';
                try
                    t1 = find(y1==mink(unique(y1), 1), 1);
                    G1 = fit(x(FMI:t1)', y1(FMI:t1)', fittype('exp1'));
                    stat(j, :) = [G1.a -1/G1.b];
                    j = j + 1;
                catch
                    disp([num2str(pd) '-' num2str(v) '-' num2str(apf)]);
                end
            end
            plot(x, median(stat(:, 1),'omitnan')*exp(-x/median(stat(:, 2),'omitnan')), ft2, 'Color', [CO6(v-20, :) 1], 'LineWidth', ti*4); %
            TXT{end+1} = [txt ', n=' num2str(size(nn, 1)) ', N=' num2str(sum(nn(:, pd+3))) ', ' num2str(round(median(stat(:, 1),'omitnan'), 1)) 'e^{-x/\bf{' num2str(round(median(stat(:, 2),'omitnan'))) '}}']; %
        end
        xlabel('');
        ylabel('');
        if mod(v, 20) > 3
            set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1, 'XTick', 0:20:40, 'YScale', 'log', 'YTickLabel', '');
        else
            set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1, 'XTick', 0:20:40, 'YScale', 'log');
        end
        legend(TXT, 'Box', 'off', 'FontSize', 14);
        axis([0 60 YMI 5]);
    end
end
set(gca, 'XTick', 0:20:60);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTick', 20:20:40);
set(nexttile(5), 'YTick', [1e-3 1e-2 0.1 1]);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig1', '-dpdf', '-r300');