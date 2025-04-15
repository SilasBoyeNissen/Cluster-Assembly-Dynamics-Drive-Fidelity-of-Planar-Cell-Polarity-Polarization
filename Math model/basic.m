%% PARAMETERS
clear; close all;
DD = ['{\bf{\color[rgb]{' num2str([27 120 55]/255) '}F}_1}'];
DP = ['{\bf{\color[rgb]{' num2str([27 120 55]/255) '}F}_2}'];
PD = ['{\bf{\color[rgb]{' num2str([153 112 171]/255) '}V}_1}'];
PP = ['{\bf{\color[rgb]{' num2str([153 112 171]/255) '}V}_2}'];
CC = repmat([0.635 0.078 0.184; 0.466 0.674 0.188; 0 0.447 0.741], 2, 1);
PR = {DD, DP, PD, PP};

%% Panel B: Time series
clc; tic;
T = 6000;
A = 0.95;
B = 0.8;
C = 0.02;
D = 0.02;
rng(2);
S = zeros(T, 5);
for i = 2:T
    S(i, :) = S(i-1, :);
    p = rand(1, 2); % find two random numbers between 0 and 1
    G = [S(i,3)*D+1             S(i,4)*D+1              S(i,1)*D+1              S(i,2)*D+1 ...
        A+(S(i,4)-S(i,1))*C     B+(S(i,3)-S(i,2))*C     A+(S(i,2)-S(i,3))*C     A+(S(i,1)-S(i,4))*C];
    G(1:4) = G(1:4) - G(5:8).*(G(5:8) < 0);
    G(5:8) = G(5:8).*(G(5:8) >= 0);
    G(S(i, :) == 0) = 0;
    S(i, 5) = S(i, 5) + 1/sum(G)*log(1/p(1)); % calculate the new time point
    if p(2) < G(1)/sum(G)
        S(i, 1) = S(i, 1) - 1;
    elseif p(2) < sum(G(1:2))/sum(G)
        S(i, 2) = S(i, 2) - 1;
    elseif p(2) < sum(G(1:3))/sum(G)
        S(i, 3) = S(i, 3) - 1;
    elseif p(2) < sum(G(1:4))/sum(G)
        S(i, 4) = S(i, 4) - 1;
    elseif p(2) < sum(G(1:5))/sum(G)
        S(i, 1) = S(i, 1) + 1;
    elseif p(2) < sum(G(1:6))/sum(G)
        S(i, 2) = S(i, 2) + 1;
    elseif p(2) < sum(G(1:7))/sum(G)
        S(i, 3) = S(i, 3) + 1;
    else
        S(i, 4) = S(i, 4) + 1;
    end
end

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 750 433]);
tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');
nexttile;
S(:, 3:4) = -S(:, 3:4);
plot(S(:, 5), S(:, 1), 'Color', '#1b7837', 'LineWidth', 1); hold on;
plot(S(:, 5), S(:, 2), 'Color', '#9970ab', 'LineWidth', 1);
plot(S(:, 5), S(:, 3), 'Color', '#5AAE61', 'LineWidth', 1);
plot(S(:, 5), S(:, 4), 'Color', '#762a83', 'LineWidth', 1);
plot([0 1e5], [0 0], 'Color', 'k', 'LineWidth', 4);
set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1, 'XTick', 0:200:1600, 'YTick', -40:20:40, 'YTickLabel', {'40', '20', '0', '20', '40'});
text(410, 40, ['\gamma = ' num2str(C) '; \delta = ' num2str(D)], 'FontSize', 22, 'HorizontalAlignment', 'center')
xlabel('Time ({\itt})', 'FontSize', 24);
ylabel('Cluster size', 'FontSize', 24);
axis([0 820 -45 45]);
saveas(gcf, 'Fig6b.png');
toc;

%% Panel C: Cluster size distribution
N = 1e4;
T = 1e5;
Av = 0.99*ones(1,6);
Bv = 0.8*ones(1,6);
Cv = [0 0 0.1 0 0 0.1];
Dv = [0 1e-3 1e-3 0 1e-3 1e-3];
generateData(N, T, Av, Bv, Cv, Dv);

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 650 400]);
tl = tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'none');
ylabel(tl, '1 - CDF', 'FontSize', 24);
for i = 1:numel(Dv)
    if mod(i, 3) == 1
        nexttile;
    end
    load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(i)*100) '-B' num2str(Bv(i)*100) '-C'  num2str(Cv(i)*10000) '-D' num2str(Dv(i)*10000)]);
    plot(0:199, 1-histcounts(RES(:, floor((i-1)/3)+3), 0:200, 'Normalization', 'cdf'), '-', 'Color', CC(i, :), 'LineWidth', 3, ...
        'DisplayName', ['\gamma=' num2str(Cv(i)) ', \delta=' num2str(Dv(i))]); hold on;
    xlabel(['Cluster size ' PR{floor((i-1)/3)+3}], 'FontSize', 24);
    set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5, 'YScale', 'log', 'XTick', 0:20:80, 'YTick', '');
    axis([0 80 0.01 5]);
end
set(nexttile(1), 'XDir', 'reverse', 'YTick', [1e-3 1e-2 0.1 1]);
legend('Box', 'off', 'Location', 'northwest', 'FontSize', 18);        
saveas(gcf, 'Fig6c.png');

%% Panel D: Polarity vs. average cluster size
N = 1e4;
T = 1e5;
Av = [0.80  0.92  0.96  0.98  0.99  0.80  0.92  0.96  0.98  0.99  0.80  0.92  0.96  0.98  0.99];
Bv = [0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80  0.80];
Cv = [0     0     0     0     0     0     0     0     0     0     0.1   0.1   0.1   0.1   0.1]*5;
Dv = [0     0     0     0     0     0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
generateData(N, T, Av, Bv, Cv, Dv);

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 400 400]);
tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');
xlabel(tl, ['Average cluster size ' PR{3} '+' PR{4}], 'FontSize', 20);
ylabel(tl, ['Polarity ' PR{4} ' / ' PR{3}], 'FontSize', 20);
nexttile;
S = zeros(numel(Dv), 2);
for i = numel(Dv):-1:1
    load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(i)*100) '-B' num2str(Bv(i)*100) '-C' num2str(Cv(i)*10000) '-D' num2str(Dv(i)*10000)]);
    S(i, :) = [mean(RES(:, 4))/mean(RES(:, 3)) mean(RES(:, 6))];
    if mod(i, 5) == 1
        plot(S(i:i+4, 2), S(i:i+4, 1), '.-', 'Color', CC((i+4)/5, :), 'LineWidth', 2.5, 'MarkerSize', 40, ...
            'DisplayName', ['\gamma=' num2str(Cv(i)) ', \delta='  num2str(Dv(i))]); hold on;
    end
end
set(gca, 'Box', 'on', 'FontSize', 19, 'LineWidth', 1.5, 'XTick', 0:20:80, 'YTick', 1:2:9);
legend('Box', 'off', 'Location', 'northwest', 'FontSize', 18);
axis([0 90 0 10], 'square')
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('Fig6d.pdf', '-dpdf');

%% Panel E: Phase diagram
N = 1e4;
T = 1e5;
Av = 0.99*ones(1, 81);
Bv = 0.8*ones(1, 81);
Cv = [zeros(1,9) 1e-4*ones(1,9) 5e-4*ones(1,9) 1e-3*ones(1,9) 5e-3*ones(1,9) 1e-2*ones(1,9) 5e-2*ones(1,9) ones(1,9)/10 ones(1,9)/2];
Dv = repmat([0 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 0.1 0.5], 1, 9);
generateData(N, T, Av, Bv, Cv, Dv);

% dx = 5;
% figure(1); clf;
% set(gcf, 'color', 'w', 'Position', [0 0 480 400]);
% S = zeros(dx, dx);
% for C = Cv(1:dx:end)
%     y = 1;
%     for D = Dv(1:5)
%         load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(1)*100) '-B' num2str(Bv(1)*100) '-C' num2str(C*10000) '-D' num2str(D*10000)]);
%         S(dx, y) = mean(RES(:, 4))/mean(RES(:, 3));
%         y = y + 1;
%     end
%     dx = dx - 1;
% end
% imagesc(S);
% colormap(flipud(jet));
% c = colorbar('Ticks', 1:2:9, 'FontSize', 20, 'Limits', [1 9]);
% c.Label.FontSize = 22;
% c.Label.String = ['Polarity ' PR{4} '/' PR{3}];
% set(gca, 'FontSize', 22, 'XTick', 1:5, 'YTick', 1:5, 'XTickLabel', Dv(1:5), 'YTickLabel', Cv(end:-5:1), 'dataAspectRatio',[1 1 1]);
% ylabel('\gamma', 'FontSize', 24, 'FontWeight', 'bold');
% xlabel('\delta', 'FontSize', 24, 'FontWeight', 'bold');
% saveas(gcf, 'Fig6e.png');

%% Panel F: Misoriented clusters
N = 1e4;
T = 1e5;
Av = [0.99 0.99  0.99];
Bv = [0.8  0.8   0.8];
Cv = [0    0     0.1];
Dv = [0    0.001 0.001];
generateData(N, T, Av, Bv, Cv, Dv);

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 400 400]);
tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');
ylabel(tl, ['Fraction of ' PR{3} '≥' PR{4}], 'FontSize', 24);
xlabel(tl, ['Cluster size ' PR{3} '+' PR{4}], 'FontSize', 24);
TX = {'''Fz[null]''', '''FmiΔcad''', '''WT'''};
nexttile;
for i = 1:numel(Dv)
    load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(i)*100) '-B' num2str(Bv(i)*100) '-C' num2str(Cv(i)*10000) '-D' num2str(Dv(i)*10000)]);
    RES(:, 7) = RES(:, 3) + RES(:, 4);
    F = (RES(:, 3) >= RES(:, 4));
    G = zeros(80, 3);
    for j = 1:80
        G(j, :) = [j sum(F(RES(:, 7) == j)) sum(RES(:, 7) == j)];
    end
    plot(G(1:end, 1), G(1:end, 2)./G(1:end, 3), '-', 'Color', CC(i, :), 'LineWidth', i, ...
        'DisplayName', ['{\bf{' TX{i}, '}} \gamma=' num2str(Cv(i)) ', \delta='  num2str(Dv(i))]); hold on;
end
set(gca, 'Box', 'on', 'FontSize', 22, 'LineWidth', 1.5, 'XTick', 0:20:80, 'YTick', 0:0.5:1);
legend('Box', 'off', 'FontSize', 17);
axis([0 80 0 1], 'square');
saveas(gcf, 'Fig6f.png');

%% Panels GH: 2-color predictions
N = 1e4;
T = 1e5;
Av = 0.99;
Bv = 0.8;
Cv = 0.1;
Dv = 0.01;
generateData(N, T, Av, Bv, Cv, Dv);

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 800 400]);
tiledlayout(1, 2, 'TileSpacing', 'default', 'Padding', 'tight');
panelGH(N, T, Av, Bv, Cv, Dv, 6, 5, 80, 80, [PD '+' PP], [DD '+' DP]);
legend(['\gamma=' num2str(Cv) ', \delta='  num2str(Dv)], 'Location', 'northwest', 'Box', 'off', 'FontSize', 19);
panelGH(N, T, Av, Bv, Cv, Dv, 4, 3, 80, 25, PP, PD);
legend(['\gamma=' num2str(Cv) ', \delta=' num2str(Dv)], 'Location', 'northeast', 'Box', 'off', 'FontSize', 19);
saveas(gcf, 'Fig6gh.png');

function panelGH(N, T, Av, Bv, Cv, Dv, x, y, xm, ym, xt, yt)
nexttile;
load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av*100) '-B' num2str(Bv*100) '-C' num2str(Cv*10000) '-D' num2str(Dv*10000)]);
scatter(RES(:, x), RES(:, y), 200, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.2);
set(gca, 'Box', 'on', 'FontSize', 22, 'LineWidth', 1.5);
xlabel(['Cluster size ' xt], 'FontSize', 24);
ylabel(['Cluster size ' yt], 'FontSize', 24);
axis([0 xm 0 ym], 'square');
end

%% Figure S9: Parameter scan
N = 1e4;
T = 1e5;
Av = [0.99*ones(1,10) 0.9*ones(1,10) 0.99*ones(1,10) 0.9*ones(1,10) 0.99*ones(1,10) 0.9*ones(1,10)];
Bv = [repmat([0 0.4 0.8 0.9 0.99], 1, 4) 0.8*ones(1,40)];
Cv = [0.1*ones(1,20) repmat([0 1e-4 1e-3 1e-2 1e-1], 1, 4) 0.1*ones(1,20)];
Dv = [1e-3*ones(1,40) repmat([0 1e-4 1e-3 1e-2 1e-1], 1, 4)];
generateData(N, T, Av, Bv, Cv, Dv);

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 1470 1100]);
tl = tiledlayout(3, 4, 'TileSpacing', 'none', 'Padding', 'none');
title(tl, ' ', 'FontSize', 24, 'FontWeight', 'bold');
xlabel(tl, 'Cluster size', 'FontSize', 24, 'FontWeight', 'bold');
ylabel(tl, '1 - CDF', 'FontSize', 24, 'FontWeight', 'bold');
TT =  repmat([4 4 4 4 4 3 3 3 3 3], 1, 6);
CG = [repmat([0 68 27; 0 109 44; 35 139 69; 65 171 93; 116 196 118], 4, 1); ...
      repmat([63 0 125; 84 39 143; 106 81 163; 128 125 186; 158 154 200], 4, 1); ...
      repmat([103 0 13; 165 15 21; 203 24 29; 239 59 44; 251 106 74], 4, 1)]/255;
for i = 1:numel(Dv)
    if mod(i, 5) == 1
        nexttile;
    end
    load(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(i)*100) '-B' num2str(Bv(i)*100) '-C'  num2str(Cv(i)*10000) '-D' num2str(Dv(i)*10000)]);
    plot(0:199, 1-histcounts(RES(:, TT(i)), 0:200, 'Normalization', 'cdf'), '-', 'LineWidth', 3, ...
        'Color', CG(i, :), 'DisplayName', [PR{TT(i)} ', \beta=' num2str(Bv(i)) ', \gamma=' num2str(Cv(i)) ...
        ', \delta=' num2str(Dv(i)) ', \mu=' num2str(round(mean(RES(:, TT(i))), 0))]); hold on;
    set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5, 'YScale', 'log', 'XTick', 0:20:40, 'YTick', '');
    legend('Box', 'off', 'FontSize', 16);
    axis([0 60 0.01 30]);
end
set(gca, 'XTick', 0:20:60);
text(-120, 4.2e8, '\alpha = 0.99', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(0, 4.2e8, '\alpha = 0.9', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTickLabel', '');
set(nexttile(5), 'YTick', [1e-3 1e-2 0.1 1], 'XTickLabel', '');
set(nexttile(9), 'YTick', [1e-3 1e-2 0.1 1]);
saveas(gcf, 'FigS9.png');

%% GENERATE DATA
function generateData(N, T, Av, Bv, Cv, Dv)
for i = 1:numel(Dv)
    tic;
    A = Av(i); B = Bv(i); C = Cv(i); D = Dv(i);
    if ~isfile(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(A*100) '-B' num2str(B*100) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat'])
        rng(0);
        RES = zeros(N, 5);
        parfor c = 1:N
            S = zeros(1, 5);
            while S(5) < T
                p = rand(1, 2);
                G = [S(3)*D+1       S(4)*D+1        S(1)*D+1        S(2)*D+1 ...
                    A+(S(4)-S(1))*C B+(S(3)-S(2))*C A+(S(2)-S(3))*C A+(S(1)-S(4))*C];
                G(1:4) = G(1:4) - G(5:8).*(G(5:8) < 0);
                G(5:8) = G(5:8).*(G(5:8) >= 0);
                G(S == 0) = 0;
                S(5) = S(5) + 1/sum(G)*log(1/p(1)); % calculate the new time point
                if p(2) < G(1)/sum(G)
                    S(1) = S(1) - 1;
                elseif p(2) < sum(G(1:2))/sum(G)
                    S(2) = S(2) - 1;
                elseif p(2) < sum(G(1:3))/sum(G)
                    S(3) = S(3) - 1;
                elseif p(2) < sum(G(1:4))/sum(G)
                    S(4) = S(4) - 1;
                elseif p(2) < sum(G(1:5))/sum(G)
                    S(1) = S(1) + 1;
                elseif p(2) < sum(G(1:6))/sum(G)
                    S(2) = S(2) + 1;
                elseif p(2) < sum(G(1:7))/sum(G)
                    S(3) = S(3) + 1;
                else
                    S(4) = S(4) + 1;
                end
            end
            RES(c, :) = S;
        end
        RES(:, 5) = sum(RES(:, 1:2), 2);
        RES(:, 6) = sum(RES(:, 3:4), 2);
        save(['_data/N' num2str(N) '-T' num2str(T) '-A' num2str(A*100) '-B' num2str(B*100) '-C' num2str(C*10000) '-D' num2str(D*10000)], 'RES');
    end
    toc;
end
end