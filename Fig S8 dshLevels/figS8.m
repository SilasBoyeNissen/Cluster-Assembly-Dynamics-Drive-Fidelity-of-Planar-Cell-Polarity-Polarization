% Figure S8: Normalized Dsh Levels
figure(1); clf;
x = ["{\bf Dsh^{WT}}" "{\bf Dsh^{N80D}}" "{\bf Dsh^{N80A}}" "{\bf Dsh^{G63D}}" "{\bf W^{1118}}"];
RES = [1 0.201719442 0.513658327; ...
       2 0.266716867 0.460258481; ...
       3 0.135395875 0.085470085; ...
       4 0.386018957 0.200790827; ...
       5 0.018368321 0];
RES(:, 4) = mean(RES(:, 2:3), 2);
bar(x, RES(:, 4)); hold on;
scatter(RES(:, 1), RES(:, 2), 200, 'filled'); hold on;
scatter(RES(:, 1), RES(:, 3), 200, 'filled'); hold on;
ylabel('Normalized Dsh Levels', 'FontWeight', 'bold');
set(gca, 'FontSize', 28);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('figS8', '-dpdf', '-r300');