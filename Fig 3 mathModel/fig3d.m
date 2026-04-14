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

% Panel D: Polarity vs. average cluster size
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
    load(['data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(i)*100) '-B' num2str(Bv(i)*100) '-C' num2str(Cv(i)*10000) '-D' num2str(Dv(i)*10000)]);
    S(i, :) = [mean(RES(:, 4))/mean(RES(:, 3)) mean(RES(:, 6))];
    if mod(i, 5) == 1
        plot(S(i:i+4, 2), S(i:i+4, 1), '.-', 'Color', CC((i+4)/5, :), 'LineWidth', 2.5, 'MarkerSize', 35, ...
            'DisplayName', ['\gamma=' num2str(Cv(i)) ', \delta='  num2str(Dv(i))]); hold on;
    end
end
set(gca, 'Box', 'on', 'FontSize', 18, 'LineWidth', 1.5, 'XTick', 0:20:80, 'YTick', 1:2:9);
legend('Box', 'off', 'Location', 'northwest', 'FontSize', 12);
axis([0 80 0 11], 'square')
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig3d.pdf', '-dpdf');

% GENERATE DATA
function generateData(N, T, Av, Bv, Cv, Dv)
for i = 1:numel(Dv)
    tic;
    A = Av(i); B = Bv(i); C = Cv(i); D = Dv(i);
    if ~isfile(['data/N' num2str(N) '-T' num2str(T) '-A' num2str(A*100) '-B' num2str(B*100) '-C' num2str(C*10000) '-D' num2str(D*10000) '.mat'])
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
        save(['data/N' num2str(N) '-T' num2str(T) '-A' num2str(A*100) '-B' num2str(B*100) '-C' num2str(C*10000) '-D' num2str(D*10000)], 'RES');
    end
    toc;
end
end