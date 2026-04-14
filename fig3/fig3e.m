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

% Panel E: Phase diagram
Av = 0.99*ones(1, 81);
Bv = 0.8*ones(1, 81);
Cv = [zeros(1,9) 1e-4*ones(1,9) 5e-4*ones(1,9) 1e-3*ones(1,9) 5e-3*ones(1,9) 1e-2*ones(1,9) 5e-2*ones(1,9) ones(1,9)/10 ones(1,9)/2];
Dv = repmat([0 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 0.1 0.5], 1, 9);
generateData(N, T, Av, Bv, Cv, Dv);

dx = 9;
figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 480 400]);
S = zeros(dx, dx);
for C = Cv(1:dx:end)
     y = 1;
     for D = Dv(1:9)
         load(['data/N' num2str(N) '-T' num2str(T) '-A' num2str(Av(1)*100) '-B' num2str(Bv(1)*100) '-C' num2str(C*10000) '-D' num2str(D*10000)]);
         S(dx, y) = mean(RES(:, 4))/mean(RES(:, 3));
         y = y + 1;
    end
    dx = dx - 1;
end
imagesc(S);
colormap(flipud(jet));
c = colorbar('Ticks', 1:2:9, 'FontSize', 18, 'Limits', [1 10]);
c.Label.FontSize = 18;
c.Label.String = ['Polarity ' PR{4} '/' PR{3}];
set(gca, 'FontSize', 18, 'XTick', 2:2:8, 'YTick', 2:2:8, 'XTickLabel', Dv(2:2:8), 'YTickLabel', Cv(end-9:-18:1), 'dataAspectRatio', [1 1 1]);
ylabel('\gamma', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('\delta', 'FontSize', 20, 'FontWeight', 'bold');
xtickangle(30)
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig3e.pdf', '-dpdf');

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