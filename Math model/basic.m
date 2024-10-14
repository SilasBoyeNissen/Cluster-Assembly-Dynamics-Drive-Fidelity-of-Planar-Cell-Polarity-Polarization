%% PARAMETERS
clear; close all;
C2 = [27 120 55; 153 112 171]/255;
CG = [0 68 27; 0 109 44; 35 139 69; 65 171 93; 116 196 118]/255;
CP = [63 0 125; 84 39 143; 106 81 163; 128 125 186; 158 154 200]/255;
CR = [103 0 13; 165 15 21; 203 24 29; 239 59 44; 251 106 74]/255;
DT = ['{\color[rgb]{' num2str(C2(1, :)) '}'];
PT = ['{\color[rgb]{' num2str(C2(2, :)) '}'];
DD = ['{\bf' DT 'F}_1}'];
DP = ['{\bf' DT 'F}_2}'];
PD = ['{\bf' PT 'V}_1}'];
PP = ['{\bf' PT 'V}_2}'];
Av = [0 0.2 0.4 0.45 0.495];
Bv = [0 0.2 0.4 0.45 0.495];
Cv = [0 1e-4 1e-3 1e-2 1e-1]; % 1
Dv = [0 1e-4 1e-3 1e-2 1e-1]; % 1
Tv = [1 1e1 1e2 1e3 1e4 1e5 1e6];

Av = [0.41 0.42 0.43 0.44 0.46 0.47 0.48 0.49];
Bv = 0.4;
Cv = 0.1; % 1
Dv = 0.001; % 1
Tv = 1e6;

PR = {DD, DP, PD, PP};

%% GENERATE data set
%parpool('local', str2double(getenv('SLURM_CPUS_PER_TASK'))); maxNumCompThreads(16);
clc; tic;
N = 1e3;

for m = 1%:2
    for T = Tv
        for A = Av
            for B = Bv
                toc; disp(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) ' started']); tic;
                for C = Cv
                    for D = Dv
                        if ~isfile(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat'])
                            rng(0);
                            RES = zeros(N, 5);
                            parfor c = 1:N
                                S = zeros(1, 5);
                                while S(5) < T
                                    p = rand(1, 2);
                                    if m == 1
                                        G = [S(3)*D+1-A         S(4)*D+1-A          S(1)*D+1-A          S(2)*D+1-A ...
                                            A+(S(4)-S(1))*C     B+(S(3)-S(2))*C     A+(S(2)-S(3))*C     A+(S(1)-S(4))*C];
                                    elseif m == 2
                                        G = [1-A                1-A                 1-A                 1-A ...
                                            A+(S(4)-S(1))*C     B+(S(3)-S(2))*C     A+(S(2)-S(3))*C     A+(S(1)-S(4))*C];
                                        G = G + [S(3:4)>0 S(1:2)>0 0 0 0 0]*D;
                                    end
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
                            save(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000)], 'RES');
                        end
                    end
                end
            end
        end
    end
end

%% GENERATE time series
clc; tic;
T = 30000;
A = [0.495 0.495 0.495 0.495];
B = [0.495 0.4 0.4 0.4];
C = [0 0 0.01 0.01];
D = [0 0 0 0.01];

for m = 1:2
    for c = 1:numel(A)
        rng(2);
        S = zeros(T, 5);
        for i = 2:T
            S(i, :) = S(i-1, :);
            p = rand(1, 2); % find two random numbers between 0 and 1
            if m == 1
                G = [S(i,3)*D(c)+1-A(c)       S(i,4)*D(c)+1-A(c)        S(i,1)*D(c)+1-A(c)        S(i,2)*D(c)+1-A(c) ...
                    A(c)+(S(i,4)-S(i,1))*C(c) B(c)+(S(i,3)-S(i,2))*C(c) A(c)+(S(i,2)-S(i,3))*C(c) A(c)+(S(i,1)-S(i,4))*C(c)];
            elseif m == 2
                G = [1-A(c)                   1-A(c)                    1-A(c)                    1-A(c) ...
                    A(c)+(S(i,4)-S(i,1))*C(c) B(c)+(S(i,3)-S(i,2))*C(c) A(c)+(S(i,2)-S(i,3))*C(c) A(c)+(S(i,1)-S(i,4))*C(c)];
                G = G + [S(i, 3:4)>0 S(i, 1:2)>0 0 0 0 0]*D(c);
            end
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
        save(['_time' num2str(m) '/T' num2str(T) '-A' num2str(A(c)*1000) '-B' num2str(B(c)*1000) '-C' num2str(C(c)*10000) '-D' num2str(D(c)*10000)], 'S');
    end
end

%% Panel B: Time series
m = 1;
c = 4;
figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 600 433]);
tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');
xlabel(tl, 'Reduced time ({\itt})', 'FontSize', 24);
ylabel(tl, 'Cluster size', 'FontSize', 24);
nexttile;
load(['_time' num2str(m) '/T' num2str(T) '-A' num2str(A(c)*1000) '-B' num2str(B(c)*1000) '-C' num2str(C(c)*10000) '-D' num2str(D(c)*10000) '.mat']);
S(:, 3:4) = -S(:, 3:4);
plot([0 1e5], [0 0], 'Color', 'k', 'LineWidth', 3); hold on;
plot(S(:, 5), S(:, 1), 'Color', '#5AAE61', 'LineWidth', 1);
plot(S(:, 5), S(:, 2), 'Color', '#A16EAF', 'LineWidth', 1);
plot(S(:, 5), S(:, 3), 'Color', '#5AAE61', 'LineWidth', 1);
plot(S(:, 5), S(:, 4), 'Color', '#A16EAF', 'LineWidth', 1);
set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1, 'YTick', '');
title(['\gamma = ' num2str(C(c)) '; \delta = ' num2str(D(c))], ...
    'Units', 'normalized', 'Position', [0.5, 0.9, 0], 'HorizontalAlignment', 'center', 'FontWeight', 'normal');
%[~, hobj] = legend(PR, 'Box', 'off', 'Location', 'southwest', 'Orientation', 'horizontal');
%set(findobj(hobj, 'type', 'line'), 'LineWidth', 2);
axis([0 2000 -40 40]);
set(gca, 'XTick', 0:1000:5000);
set(nexttile(1), 'XTick', 0:500:4000, 'YTick', -40:20:40, 'YTickLabel', {'40', '20', '0', '20', '40'});
saveas(gcf, ['Fig5b-m' num2str(m) '.png']);

%% Panel C: Cluster size distribution of P+D and P/D
PS = 2;
N = 1e3;
T = 1e6;
A = 0.495;
B = 0.4;
Ct = [0 0     0.1];
Dt = [0 0.01 0.01];
CL = [0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0 0.4470 0.7410];
TX = {'''Fz[null]''', '''FmiΔcad''', '''WT'''};
m = 1;
figure(4); clf;
set(gcf, 'color', 'w', 'Position', [0 0 650 400]);
tl = tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'none');
ylabel(tl, '1 - CDF', 'FontSize', 24);
LT = {[DD '+' DP], DD, DP, [PP '+' PD], PP, PD};
TT = [5 1 2 6 4 3];
XX = 0:1:200;
for c = 3:-1:2
    nexttile;
    for i = 1:numel(Ct)
        C = Ct(i); D = Dt(i);
        load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C'  num2str(C*10000) '-D' num2str(D*10000) '.mat']);
        plot(XX(1:end-1), 1-histcounts(RES(:, TT(c+(PS-1)*3)), XX, 'Normalization', 'cdf'), '-', 'Color', CL(i, :), 'LineWidth', 3, ...
            'DisplayName', ['\gamma=' num2str(C) ', \delta=' num2str(D) ', \mu=' num2str(round(mean(RES(:, TT(c+(PS-1)*3))), 0))]); hold on;
    end
    xlabel(['Cluster size ' LT{c+(PS-1)*3}], 'FontSize', 24);
    set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5, 'YScale', 'log', 'XTick', 0:20:40, 'YTick', '');
    axis([0 60 0.01 5]);
    legend('Box', 'off', 'FontSize', 20);
end
set(gca, 'XTick', 0:20:60);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1]);
saveas(gcf, ['Fig5c-m' num2str(m) '.png']);

%% Panel D: 2-color predictions, D vs P (same protein) and D vs P (other protein)
N = 1e3;
T = 1e6;
A = 0.495;
B = 0.4;
C = 0.1;
D = 0.01;

m = 1;
figure(5); clf;
set(gcf, 'color', 'w', 'Position', [0 0 750 400]);
tiledlayout(1, 2, 'TileSpacing', 'default', 'Padding', 'tight');
calcE(N, T, A, B, C, D, m, 5, 6, 90, 90, [DD '+' DP], [PD '+' PP]);
legend(['\gamma=' num2str(C) ', \delta=' num2str(D)], 'Location', 'northwest', 'Box', 'off')
calcE(N, T, A, B, C, D, m, 4, 3, 90, 25, PP, PD);
legend(['\gamma=' num2str(C) ', \delta=' num2str(D)], 'Location', 'northeast', 'Box', 'off')
saveas(gcf, ['Fig5d-m' num2str(m) '.png']);

function calcE(N, T, A, B, C, D, m, x, y, xm, ym, xt, yt)
nexttile;
load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat']);
scatter(RES(:, x), RES(:, y), 160, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.2);
set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5);
xlabel(['Cluster size ' xt], 'FontSize', 24);
ylabel(['Cluster size ' yt], 'FontSize', 24);
axis([0 xm 0 ym], 'square');
end

%% Panel E: Phase diagram of gamma and beta at the mean cluster size (Fmi delta cad in one corner, complete null in lower left corner)
N = 1e3;
T = 1e6;
A = 0.495;
B = 0.4;
PS = [4 3];
m = 1;
x = 5;
figure(7); clf;
set(gcf, 'color', 'w', 'Position', [0 0 500 500]);
S = zeros(numel(Cv), numel(Dv));
for C = Cv
    y = 1;
    for D = Dv
        load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat']);
        S(x, y) = mean(RES(:, PS(1)))/mean(RES(:, PS(2)));
        y = y + 1;
    end
    x = x - 1;
end
imagesc(S);
colormap(flipud(jet));
c = colorbar('Ticks', 1:5, 'FontSize', 20);
c.Label.FontSize = 24;
c.Label.String = ['{\bf Polarity (' PR{PS(1)} '/' PR{PS(2)} ')}']; % [PR{PS(1)} '{\bf mean cluster size (\mu)}'];
set(gca, 'FontSize', 20, 'XTick', 1:5, 'YTick', 1:5, 'XTickLabel', Dv, 'YTickLabel', Cv(end:-1:1))
ylabel('\gamma', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('\delta', 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, ['Fig5e-m' num2str(m) '.png']);

%% Panel F: Polarity, error rate as a function of beta (fraction of misoriented clusters)
PS = [3 4];
N = 1e3;
T = 1e6;
m = 1;
A = 0.495;
B = 0.4;
C = [0     0   0.1];
D = [0     0.001     0.001];
figure(6); clf;
set(gcf, 'color', 'w', 'Position', [0 0 465 465]);
tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');
ylabel(tl, ['Fraction of ' PR{PS(1)} '≥' PR{PS(2)}], 'FontSize', 24, 'FontWeight', 'bold');
xlabel(tl, ['Cluster size (' PR{PS(1)} '+' PR{PS(2)} ')'], 'FontSize', 24, 'FontWeight', 'bold');
title(tl, ['\alpha = ' num2str(A)], 'FontSize', 24, 'FontWeight', 'bold');
nexttile;

for c = 1:3
    load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C(c)*10000) '-D' num2str(D(c)*10000) '.mat']);
    RES(:, 7) = RES(:, PS(1)) + RES(:, PS(2));
    F = (RES(:, PS(1)) >= RES(:, PS(2)));
    G = zeros(60, 3);
    for i = 1:60
        G(i, :) = [i sum(F(RES(:, 7) == i)) sum(RES(:, 7) == i)];
    end
    plot(G(1:1:end, 1), G(1:1:end, 2)./G(1:1:end, 3), '-', 'LineWidth', 1.5, ...
        'DisplayName', ['\beta=' num2str(B) ', \gamma=' num2str(C(c)) ', \delta='  num2str(D(c))]); hold on;
end

set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5, 'XTick', 0:10:40);
legend('Box', 'off');
axis([0 60 0 1]);
set(gca, 'XTick', 0:10:100);
saveas(gcf, ['Fig5f-m' num2str(m) '.png']);

%% Figure S4: Parameter scan
PS = 2;
N = 1e3;
T = 1e6;
A = 0.45;
m = 1;

figure(4); clf;
set(gcf, 'color', 'w', 'Position', [0 0 750 1100]);
tl = tiledlayout(3, 2, 'TileSpacing', 'none', 'Padding', 'none');
title(tl, ['\alpha = ' num2str(A)], 'FontSize', 24, 'FontWeight', 'bold')
xlabel(tl, 'Cluster size', 'FontSize', 24, 'FontWeight', 'bold');
ylabel(tl, '1 - CDF', 'FontSize', 24, 'FontWeight', 'bold');
LS = {'-', '-'};
LT = {DD, DP, PP, PD};
TT = [1 2 4 3];
XX = 0:1:200;
for i = 1:6
    nexttile;
    if i < 3
        j = 0;
        for B = Bv
            j = j + 1; C = 0.1; D = 0.001;
            load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C'  num2str(C*10000) '-D' num2str(D*10000) '.mat']);
            plot(XX(1:end-1), 1-histcounts(RES(:, TT(i+(PS-1)*2)), XX, 'Normalization', 'cdf'), LS{i}, 'LineWidth', 3, ...
                'Color', CG(j, :), 'DisplayName', [LT{i+(PS-1)*2} ', \beta=' num2str(B) ', \gamma=' num2str(C) ...
                ', \delta=' num2str(D) ', \mu=' num2str(round(mean(RES(:, TT(i+(PS-1)*2))), 0))]); hold on;
        end
    elseif i < 5
        j = 0;
        for C = Cv
            j = j + 1; B = 0.4; D = 0.001;
            load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat']);
            plot(XX(1:end-1), 1-histcounts(RES(:, TT(i-2+(PS-1)*2)), XX, 'Normalization', 'cdf'), LS{i-2}, 'LineWidth', 3, ...
                'Color', CP(j, :), 'DisplayName', [LT{i-2+(PS-1)*2} ', \beta=' num2str(B) ', \gamma=' num2str(C) ...
                ', \delta=' num2str(D) ', \mu=' num2str(round(mean(RES(:, TT(i-2+(PS-1)*2))), 0))]); hold on;
        end
    else
        j = 0;
        for D = Dv
            j = j + 1; B = 0.4; C = 0.1;
            load(['_data' num2str(m) '-N' num2str(N) '/T' num2str(T) '-A' num2str(A*1000) '-B' num2str(B*1000) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat']);
            plot(XX(1:end-1), 1-histcounts(RES(:, TT(i-4+(PS-1)*2)), XX, 'Normalization', 'cdf'), LS{i-4}, 'LineWidth', 3, ...
                'Color', CR(j, :), 'DisplayName', [LT{i-4+(PS-1)*2} ', \beta=' num2str(B) ', \gamma=' num2str(C) ...
                ', \delta=' num2str(D) ', \mu=' num2str(round(mean(RES(:, TT(i-4+(PS-1)*2))), 0))]); hold on;
        end
    end
    set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1.5, 'YScale', 'log', 'XTick', 0:20:40, 'YTick', '');
    axis([0 60 0.001 200]);
    legend('Box', 'off', 'FontSize', 16);
end
set(gca, 'XTick', 0:20:60);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTickLabel', '');
set(nexttile(3), 'YTick', [1e-3 1e-2 0.1 1], 'XTickLabel', '');
set(nexttile(5), 'YTick', [1e-3 1e-2 0.1 1]);
saveas(gcf, 'FigS3.png');