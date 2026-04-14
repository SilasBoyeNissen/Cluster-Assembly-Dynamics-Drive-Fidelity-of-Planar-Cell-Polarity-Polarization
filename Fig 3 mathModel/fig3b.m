% GLOBAL PARAMETERS
clear; close all;
N = 1e4;
T = 1e5;
DD = ['{\bf{\color[rgb]{' num2str([27 120 55]/255) '}F}_1}'];
DP = ['{\bf{\color[rgb]{' num2str([27 120 55]/255) '}F}_2}'];
PD = ['{\bf{\color[rgb]{' num2str([153 112 171]/255) '}V}_1}'];
PP = ['{\bf{\color[rgb]{' num2str([153 112 171]/255) '}V}_2}'];
CC = repmat([0.635 0.078 0.184; 0.466 0.674 0.188; 0 0.447 0.741], 2, 1);
PR = {DD, DP, PD, PP};

% Panel B: Time series
clc; tic;
Tv = 7000;
A = 0.99;
B = 0.8;
C = 0.5;
D = 0.001;
rng(2);
S = zeros(Tv, 5);
for i = 2:Tv
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
set(gca, 'Box', 'on', 'FontSize', 18, 'LineWidth', 1, 'XTick', 0:200:1600, 'YTick', -40:20:40, 'YTickLabel', {'40', '20', '0', '20', '40'});
text(410, 40, ['\gamma = ' num2str(C) '; \delta = ' num2str(D)], 'FontSize', 18, 'HorizontalAlignment', 'center')
xlabel('Time ({\itt})', 'FontSize', 20);
ylabel('Cluster size', 'FontSize', 20);
axis([0 800 -45 45]);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig3b.pdf', '-dpdf');
toc;