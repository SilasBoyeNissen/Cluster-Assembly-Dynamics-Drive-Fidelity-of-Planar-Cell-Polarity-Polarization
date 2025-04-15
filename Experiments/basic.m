%% REGISTER channels (0)
clear; tic;
FILE = dir('1-TIF/c*.tif');
for n = 1%:size(FILE, 1)
    F = FILE(n).name(1:end-4);
    CH1 = tiffreadVolume(['1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [1 1]});
    CH2 = tiffreadVolume(['1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [2 2]});
    fixedRefObj = imref2d(size(CH2));
    movingRefObj = imref2d(size(CH1));
    [optimizer, metric] = imregconfig('multimodal');
    tX = mean(fixedRefObj.XWorldLimits) - mean(movingRefObj.XWorldLimits);
    tY = mean(fixedRefObj.YWorldLimits) - mean(movingRefObj.YWorldLimits);
    iT = affinetform2d();
    iT.A(1:2,3) = [tX ; tY];
    tF = imregtform(CH1,movingRefObj,CH2,fixedRefObj,'translation',optimizer,metric,'PyramidLevels',3,'InitialTransformation',iT);
    disp(['File' num2str(n) ' ' F ': ' num2str([tF.Translation sqrt(tF.Translation(1)^2 + tF.Translation(2)^2)])]);
    MOVINGREG.Transformation = tF;
    MOVINGREG.RegisteredImage = imwarp(CH1, movingRefObj, tF, 'OutputView', fixedRefObj, 'SmoothEdges', true);
    MOVINGREG.SpatialRefObj = fixedRefObj;
    CH1all = tiffreadVolume(['1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [1 2 inf]});
    CH2all = tiffreadVolume(['1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [2 2 inf]});
    CH1new = imwarp(CH1all, movingRefObj, tF, 'OutputView', fixedRefObj, 'SmoothEdges', true);
    imwrite(MOVINGREG.RegisteredImage, 'myFile.TIFF', 'Compression', 'none');
    for i = 2:1000%size(CH1new, 3)
        imwrite(CH1new(:, :, i), 'myFile.TIFF', 'writemode', 'append', 'Compression', 'none');
    end
    imwrite(CH1new(:, :, 1001), 'myFile2.TIFF', 'Compression', 'none');
    for i = 1002:2000%size(CH1new, 3)
        imwrite(CH1new(:, :, i), 'myFile2.TIFF', 'writemode', 'append', 'Compression', 'none');
    end
    imwrite(CH1new(:, :, 2001), 'myFile3.TIFF', 'Compression', 'none');
    for i = 2002:size(CH1new, 3)
        imwrite(CH1new(:, :, i), 'myFile3.TIFF', 'writemode', 'append', 'Compression', 'none');
    end
end

%% RENUMBER files (1)
close all; clear;
k = 1;
FO = '1-TIF';
%FO = '2-ROI';
%FO = '7-res';
FILE = dir(['../../' FO '/c*']);
for n = 1:size(FILE, 1)
    F = FILE(n).name;
    try
        movefile(['../../' FO '/' F], ['../../' FO '/' F(1:50) num2str(k, '%03.0f') F(54:end)]);
        k = k + 1;
        disp(F);
    catch
    end
end
disp(['Renumbered ' num2str(k) ' files']);

%% RENUMBER res-folder
FO = '7-res';
k = 1;
FILE = dir(['../../' FO '/c*']);
for n = 1:size(FILE, 1)
    F = FILE(n).name;
    load(['../../' FO '/' F]);
    res(:, 1) = k;
    k = k + 1;
    save(['../../' FO '/' F], 'res');
end

%% DRAW polyline (2)
close all; clear;
ROI = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 486
    F = FILE(n).name(1:end-4);
    ch = 1;
    disp(F);
    figure(1); clf;
    try
        load(['../../3-MAP/' F '.mat']);
        IMG = imadjust(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch ch]}));
        imshow(labeloverlay(IMG, MAP, 'Colormap', [0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1; 1 0 1; 0 0 0], 'Transparency', 0.75));
    catch
        imshow(imadjust(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch ch]})));
    end
    set(gcf, 'Units', 'normalized', 'Outerpos', [0 0.1 1 0.9]); set(gca, 'Position', [0 0 1 1]);
    try
        load(['../../2-ROI/' F '.mat'], 'ROI');
        TEMP = [];
        for i = 1:numel(ROI)
            try
                TEMP = [TEMP; drawpolyline('Color', ROI(i, 1).Color, 'Position', ROI(i, 1).Position)];
            catch
            end
        end
        ROI = TEMP;
        pause();
    catch
    end
    while 1 % ctrl-c
        save(['../../2-ROI/' F '.mat'], 'ROI');
        ROI = [ROI; drawpolyline('color', input(['Press r (red) for proximal, g (green) for anterior or posterior, ' ...
            'b (blue) for distal, y (yellow) for proximal-distal, or c (cyan) for anterior-posterior'], 's'))]; % keyboard input
    end
end

%% CALC map (3)
close all; clear;
W = 3;%9;
%delete('../../3-MAP/*')
FILE = dir('../../1-TIF/c*.tif');
for n = 472:486%1:size(FILE, 1)
    F = FILE(n).name(1:end-4);
    disp(['n=' num2str(n) ': ' F]);
    load(['../../2-ROI/' F '.mat']);
    RAW = tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [1 1]});
    MAP = zeros(size(RAW));
    MAPa = zeros(size(RAW));
    MAPc = zeros(size(RAW));
    MAPp = zeros(size(RAW));
    MAPd = zeros(size(RAW));
    MAPf = zeros(size(RAW));
    for i = 1:numel(ROI)
        if sum(ROI(i, 1).Color) == 2 && ROI(i, 1).Color(1) == 1 % Yellow = Proximal-distal (or null mutant)
            MAPf = MAPf + imdilate(createMask(ROI(i,1), RAW), ones(W, W));
        elseif sum(ROI(i, 1).Color) == 2 && ROI(i, 1).Color(3) == 1 % Cyan = Anterior-posterior
            MAPc = MAPc + imdilate(createMask(ROI(i,1), RAW), ones(W, W));
        elseif ROI(i, 1).Color(1) == 1 % Red = Proximal
            MAPp = MAPp + imdilate(createMask(ROI(i,1), RAW), ones(W, W));
        elseif ROI(i, 1).Color(2) == 1 % Green = Anterior or posterior
            MAPa = MAPa + imdilate(createMask(ROI(i,1), RAW), ones(W, W));
        elseif ROI(i, 1).Color(3) == 1 % Blue = Distal
            MAPd = MAPd + imdilate(createMask(ROI(i,1), RAW), ones(W, W));
        end
    end
    MAP(MAPa > 0) = 2;
    MAP(MAPc > 0) = 5;
    MAP(MAPd > 0) = 1;
    MAP(MAPp > 0) = 3;
    MAP(MAPf > 0) = 4;
    save(['../../3-MAP/' F '.mat'], 'MAP');
end

%% SAVE tif+map (4+5)
close all; clear;
DX = 5; % 3
BACK = 5; % 2
DT = 5;
FR = 1;
TP = 0.75; % 1 % 0.75

%delete('../../4-tif/*')
%delete('../../5-map/*')
FILE = dir('../../1-TIF/c*.tif');
for n = 472:486%1:size(FILE, 1) %[168 191 205 231 235]
    F = FILE(n).name(1:end-4);
    disp(['n=' num2str(n) ': ' F]);
    load(['../../3-MAP/' F '.mat']);
    MAP(bwmorph(MAP > 0, 'remove') == 0) = 0;
    MAP = repelem(MAP, 2, 2, 1);
    for ch = 1:str2double(F(22))
        RAW1 = tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch 2 Inf]});
        RAW1 = RAW1(:, :, 1:10);
        RAW1 = repelem(RAW1, 2, 2, 1);
        RAW = single(RAW1);
        if size(RAW, 3) == 1
            RAW1 = tiffreadVolume(['../../1-TIF/' F '.tif']);
            RAW1 = RAW1(:, :, ch:2:end);
            RAW1 = RAW1(:, :, 1:10);
            RAW1 = repelem(RAW1, 2, 2, 1);
            RAW = single(RAW1);
        end
        if ch == 2
            figure('visible', 'off'); clf;
            imshow(imfuse(h, imadjust(RAW1(:, :, 1)), 'ColorChannels', [1 2 0]));
            print(['../../4-tif/' F '-3.png'], '-dpng','-r0');
        end
        figure('visible', 'off'); clf;
        set(gcf, 'Position', [0 0 size(RAW1, 1:2)]);
        h = imadjust(RAW1(:, :, 1));
        if str2double(F(22)) == 1
            imshow(imfuse(h, zeros(size(RAW1, 1:2)), 'ColorChannels', [2 1 0]));
        else
            imshow(imfuse(h, zeros(size(RAW1, 1:2)), 'ColorChannels', [ch abs(ch-3) 0]));
        end
        print(['../../4-tif/' F '-' num2str(ch) '.png'], '-dpng','-r0');

        TAK = imboxfilt(RAW - FR*imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1);
        SMOO = cumsum(TAK(:, :, 1:DT), 3);
        I = uint16(SMOO(:, :, DT)); % /10
        I(70, 70:101) = max(I, [], 'all');
        if ch == 2
            figure('visible', 'off'); clf;
            HH = imadjust(I, stretchlim(I, [0.0001 0.9999]));
            IMG = imfuse(H, HH, 'ColorChannels', [1 2 0]);
            imshow(labeloverlay(IMG, MAP, 'Colormap', [0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1; 1 0 1; 0 0 0], 'Transparency', TP));
            print(['../../5-map/' F '-3.png'], '-dpng','-r0');
        end
        figure('visible', 'off'); clf;
        set(gcf, 'Position', [0 0 size(RAW1, 1:2)]);
        H = imadjust(I, stretchlim(I, [0.0001 0.9999]));
        if str2double(F(22)) == 1
            IMG = imfuse(H, zeros(size(RAW1, 1:2)), 'ColorChannels', [2 1 0]);
            %IMG = H; % for white colors
        else
            IMG = imfuse(H, zeros(size(RAW1, 1:2)), 'ColorChannels', [ch abs(ch-3) 0]);
        end
        imshow(labeloverlay(IMG, MAP, 'Colormap', [0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1; 1 0 1; 0 0 0], 'Transparency', TP));
        print(['../../5-map/' F  '-' num2str(ch) '.png'], '-dpng','-r0');
    end
    disp(['img' num2str(n) ' took ' num2str(round(toc)) ' sec']);
end

%% MAIN script (6)
close all; clear; rng(1); tic;
DX = 3;
DT = 20;
BACK = 2;
BINW = 0.3;
FMAX = 110;

RES = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 480:486%1:size(FILE, 1)
    F = FILE(n).name(1:end-4); disp(F);
    CH = str2double(F(22));
    load(['../../3-MAP/' F '.mat']);
    DIM = [find(sum(MAP), 1) find(sum(MAP, 2), 1) find(sum(MAP), 1, 'last') find(sum(MAP, 2), 1, 'last')];
    MAP = MAP(DIM(2):DIM(4), DIM(1):DIM(3));
    for ch = 1:CH
        RAW = single(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[DIM(2) 1 DIM(4)], [DIM(1) 1 DIM(3)], [ch CH inf]}));
        if size(RAW, 3) == 1
            RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
            RAW = single(RAW(DIM(2):DIM(4), DIM(1):DIM(3), ch:CH:end));
        end
        SMOO = movmedian(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), [DT 0], 3);
        MASK = (blockproc(SMOO(:, :, 1), [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2]) & MAP);
        [r3, c3] = find(MASK); % puncta centers
        smoo = reshape(permute(SMOO, [3 1 2]), size(RAW, 3), []);
        smor = smoo(:, (MASK)); % puncta centers
        SMOO = SMOO(:, :, 1);
        res = nan(size(r3, 1), 9);
        tank = logical(triu(ones(size(smoo, 1), size(smoo, 1)), 1));
        for i = 1:size(smor, 2) % parfor
            wi = smor(:, i) - smor(:, i)';
            [pwr, f] = pspectrum(histcounts(wi(tank), 'BinWidth', BINW), 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/40]);
            I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
%            if ~isempty(I)
%                res(i, :) = [n str2double(F(17:18)) str2double(F(2:3)) BINW*2*pi/f(I) smor(1, i)/(BINW*2*pi/f(I)) c3(i) r3(i) MAP(r3(i), c3(i)) ch];
%            else %
                [pddf, edges] = histcounts(wi(tank), 'BinWidth', BINW, 'Normalization', 'probability');
                [pwr, f] = pspectrum(pddf, 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/40]);
                I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
                stepSize(F, i, smor, edges, pddf, f, pwr, BINW*2*pi/f(I), ch);
                figure(2);
                figure(1);
                w = waitforbuttonpress;
                if w
                    cs = get(gcf, 'CurrentCharacter');
                    disp([ num2str(n) ': ' num2str(i) ' out of ' num2str(size(smor, 2))])
                end
                res(i, :) = [n str2double(F(17:18)) str2double(F(2:3)) 0 str2double(cs)-1 c3(i) r3(i) MAP(r3(i), c3(i)) ch];
 %           end
        end
        save(['../../7-res/' F '.mat'], 'res'); %
        res(isnan(res(:, 4)), :) = [];
        RES = [RES; res];
    end
    disp(['img' num2str(n) ' took ' num2str(round(toc)) ' sec']);
end
%save('RES.mat', 'RES');

%% Figure 0: LOAD & Prepare
close all; clear;
RES = [];
FILE = dir('../../7-res/c*');
for n = 1:458
    load(['../../7-res/' FILE(n).name]);
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
FILE = dir('../../1-TIF/c*.tif');
for i = 1:471%unique(RES(:, 1))'
    F = FILE(i).name(1:end-4); 
    load(['../../3-MAP/' F '.mat']);
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
    N(N(:, 1) == i, 16:21) = N(N(:, 1) == i, 4:9)./(histcounts(MAP,0.5:6.5)/9*65/1000);
end
N(isnan(N)) = 0;
Nl(Nl == 0) = NaN;
Nle(Nle == 0) = NaN;
XT = 15:5:30;
YT = [0:10:60; 0:5:30; 0:1:6];

%% Figure 1: Single exponentials
pd = 4;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, '1 - Cumulative Distribution Function', 'Fontweight', 'bold', 'FontSize', 24);
for v = [11 12 14 13 15 16]
    TXT = {};
    nt = nexttile;
    for apf = 1:2
        if apf == 1
            nn = N(N(:, 3) == v & N(:, 2) < 24, :);
            txt = ' 15-23 hr APF';
            ft1 = '--';
            ft2 = ':';
            ti = 2;
        else
            nn = N(N(:, 3) == v & N(:, 2) > 23, :);
            txt = ' 24-32 hr APF';
            ft1 = '-';
            ft2 = '-';
            ti = 1;
        end
        stat = zeros(size(nn, 1), 2);
        j = 1;
        for i = nn(:, 1)'
            TXT{end+1} = '';
            y1 = 1-histcounts(RES(RES(:, 1) == i & RES(:, 8) == pd, 5), x, 'Normalization', 'cdf');
            plot(x(1:end-1), y1, ft1, 'Color', [CO6(v-10, :) 1], 'LineWidth', ti); hold on;
            t1 = find(y1==mink(unique(y1), 1), 1);
            G1 = fit(x(FMI:t1)', y1(FMI:t1)', fittype('exp1'));
            stat(j, :) = [G1.a -1/G1.b];
            j = j + 1;
        end
        plot(x, median(stat(:, 1))*exp(-x/median(stat(:, 2))), ft2, 'Color', [CO6(v-10, :) 1], 'LineWidth', ti*4); %
        TXT{end+1} = [txt ', n=' num2str(size(nn, 1)) ', N=' num2str(sum(nn(:, pd+3))) ', ' num2str(round(median(stat(:, 1)), 1)) 'e^{-x/\bf{' num2str(round(median(stat(:, 2)))) '}}']; %
    end
    xlabel('');
    ylabel('');
    if mod(v, 10) == 1 || mod(v, 10) == 3
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 2, 'XTick', 0:20:80, 'YScale', 'log');
    else
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 2, 'XTick', 0:20:80, 'YScale', 'log', 'YTickLabel', '');
    end
    legend(TXT, 'Box', 'off', 'FontSize', 21);
    axis([0 100 YMI 2]);
end
set(gca, 'XTick', 0:20:100);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTick', 20:20:80);
set(nexttile(4), 'YTick', [1e-3 1e-2 0.1 1]);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig1', '-dpdf', '-r300');

%% Figure 2: Clusters grow with age
pd = 4;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 920]);
tl = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 24);
title(tl, ' ');
for v = [11 12 14 13 15 16]
    nt = nexttile;
    errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3); hold on;
    mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
    [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
    plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
    fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    xlabel('');
    ylabel('');
    if mod(v, 10) == 1 || mod(v, 10) == 3
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT);
    else
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
    end
    axis([14 33 0 44.99]);
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
        '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
        '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
        '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
        'FontSize', 20, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
set(gca, 'XTick', XT);
set(nexttile(1), 'YTick', YT(1, :), 'XTick', XT);
set(nexttile(4), 'YTick', YT(1, :));
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig2', '-dpdf', '-r300');

%% Figure 3: Clone boundary
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 24);
for pd = [1 3]
    for v = 23:26
        nt = nexttile;
        errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3); hold on;
        mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
        [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
        plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
        fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xlabel('');
        ylabel('');
        if mod(v, 20) > 3
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
        else
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT);
        end
        axis([14 33 0 64.99]);
        CI = coefCI(mdl);
        [Ypred, Yci] = predict(mdl, 23.5);
        text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
            '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
            '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
            '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
            'FontSize', 20, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end
end
set(gca, 'XTick', XT);
set(nexttile(1), 'YTick', YT(1, :), 'XTick', XT);
set(nexttile(5), 'YTick', YT(1, :));
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig3', '-dpdf', '-r300');

%% Figure 4: Loss-of-function
pd = 4;
figure(1); clf; j = 1;
set(gcf, 'Color', 'w', 'Position', [10 10 1600 820]);
tl = tiledlayout(2, 5, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 24);
title(tl, ' ');
for v = [13 33 15 35 16 36 14 38 14 34]
    cc = 0.4*ones(1, 3);
    nt = nexttile(j);
    f = v;
    if v == 38
        f = 34;
    end
    if v < 20
        j = j + 5;
        cc = CO6(mod(f, 10), :);
    else
        tbl = table(RES(RES(:, 3) == f-20 | RES(:, 3) == v, 3), RES(RES(:, 3) == f-20 | RES(:, 3) == v, 2), RES(RES(:, 3) == f-20 | RES(:, 3) == v, 5),'VariableNames',{'GeneticID','Age','ClusterSize'});
        lme = fitlme(tbl,'ClusterSize~Age+GeneticID+(Age|GeneticID)');
        j = j - 4;
    end
    errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', cc, 'LineWidth', 3); hold on;
    mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
    [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
    plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', cc, 'LineWidth', 3);
    fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], cc, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    xlabel('');
    ylabel('');
    if mod(v, 10) == 3
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT);
    else
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
    end
    axis([14 33 0 39.99]);
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    if v > 20
        text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
            '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
            '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
            '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) newline newline ...
            '{\bf Different from WT?}' newline ' P''-value = ' num2str(round(double(lme.Coefficients(2, 6)), 1, 'significant'))], ...
            'FontSize', 18, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    else
        text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
            '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
            '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
            '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3)))], ...
            'FontSize', 18, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end
end
set(gca, 'XTick', XT);
set(nexttile(6), 'YTick', YT(1, :));
set(nexttile(1), 'YTick', YT(1, :), 'XTick', XT);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig4', '-dpdf', '-r300');

%% Figure 5: 2-color analyze
clear; rng(1); tic;
DX = 5;
BACK = 5;
DT = 5;

DY = 4;
STAT = zeros(0, 5);
FILE = dir('../../1-TIF/c*.tif');
for n = 393:458
    F = FILE(n).name(1:end-4);
    disp(['n=' num2str(n) ': ' F]);
    load(['../../3-MAP/' F '.mat']);
    MAP = repelem(MAP, 2, 2, 1);
    CH1 = img(F, 1, DT, BACK, DX);
    CH2 = img(F, 2, DT, BACK, DX);
    for pd = unique(MAP(MAP>0))'
        ch1 = CH1;
        ch2 = CH2;
        ch1((blockproc(ch1, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP == pd) == 0) = 0;
        ch2((blockproc(ch2, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP == pd) == 0) = 0;
        n1 = imdilate(ch1, ones(DY, DY));
        n2 = imdilate(ch2, ones(DY, DY));
        n12 = unique([n1(n1>0 & n2>0) n2(n1>0 & n2>0)], 'rows');
        STAT = [STAT; n*ones(size(n12, 1), 1) str2double(F(2:3))*ones(size(n12, 1), 1) pd*ones(size(n12, 1), 1) n12];
    end
    disp(['img' num2str(n) ' took ' num2str(round(toc)) ' sec']);
end
save('STAT-fig5.mat', 'STAT');

function I = img(F, ch, DT, BACK, DX)
RAW = single(repelem(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch 2 2*DT]}), 2, 2, 1));
if size(RAW, 3) == 1
    RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
    RAW = single(repelem(RAW(:, :, ch:2:2*DT), 2, 2, 1));
end
I = uint32(sum(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), 3));
end

%% Figure 5: 2-color visualize
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [150 0 200 620]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'none');
title(tl, ' ');
figure5(73, 4, 'Vang', 'Fmi', 1.5, 4.4, 0:0.5:1.5, 0);
figure5(71, 4, 'Fz', 'Vang', 1.5, 7, 0:0.5:1.5, 0);
figure5(72, 4, 'Vang', 'Pk', 1.5, 5, 0:0.5:1.5, 0);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig5a', '-dpdf', '-r300');

figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [400 0 200 600]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'none');
figure5(78, 1, 'Distal Vang', 'Proximal Fmi', 0.6, 4.4, 0:0.3:0.6, 1);
figure5(76, 1, 'Fz', 'Proximal Vang', 0.6, 30, 0:0.3:0.6, 0);
figure5(77, 1, 'Distal Vang', 'Proximal Pk', 0.6, 5, 0:0.3:0.6, 1);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig5b', '-dpdf', '-r300');

figure(3); clf;
set(gcf, 'Color', 'w', 'Position', [650 0 200 600]);
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'none');
figure5(78, 3, 'Proximal Vang', 'Distal Fmi', 2, 4.4, 0:2, 0);
figure5(76, 3, 'Fz', 'Distal Vang', 0.6, 3, 0:0.3:0.6, 1);
figure5(77, 3, 'Proximal Vang', 'Distal Pk', 2, 0.5, 0:2, 1);
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig5c', '-dpdf', '-r300');

function figure5(v, pd, xlab, ylab, xmax, ymax, xt, t)
nexttile;
load('STAT-fig5.mat');
n = numel(unique(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 1)));
n1 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 4))/1e4;
n2 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == pd, 5))/1e4;
scatter(n1, n2, 20, 'filled', 'MarkerFaceAlpha', 0.3);
R = corrcoef(n1, n2);
set(gca, 'Box', 'on', 'FontSize', 17, 'LineWidth', 1, 'XTick', xt);
if t == 0
    text(xmax/2, ymax*1.03, ['n=' num2str(n) ', N=' num2str(numel(n1)) ', r=' num2str(R(2, 1), 1)], 'FontSize', 14, 'HorizontalAlignment', 'center');
else
    text(xmax*0.96, ymax*0.87, ['n=' num2str(n) newline 'N=' num2str(numel(n1)) newline 'r=' num2str(R(2, 1), 1)], 'FontSize', 14, 'HorizontalAlignment', 'right');
end
axis([0 xmax 0 ymax*1.1], 'square');
xlabel(xlab, 'FontWeight', 'bold', 'Color', '#EE220C');
ylabel(ylab, 'FontWeight', 'bold', 'Color', '#1DB100');
end

%% Figure 7: DIX
pd = 4;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [10 10 1600 800]);
tl = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 24);
k = 0;
vv = [61 64 63 62 46 0 48 47];
for v = vv
    k = k + 1;
    nt = nexttile;
    if v == 46
        k = 1;
    end
    if v == 47 || v == 48 || v > 61
        if v > 61
            tbl = table(RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 3), RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 2), RES(RES(:, 3) == vv(k-1) | RES(:, 3) == vv(k), 5),'VariableNames',{'GeneticID','Age','ClusterSize'});
        elseif v == 47
            tbl = table(RES(RES(:, 3) == 47 | RES(:, 3) == 48, 3), RES(RES(:, 3) == 47 | RES(:, 3) == 48, 2), RES(RES(:, 3) == 47 | RES(:, 3) == 48, 5),'VariableNames',{'GeneticID','Age','ClusterSize'});
        elseif v == 48
            tbl = table(RES(RES(:, 3) == 48 | RES(:, 3) == 46, 3), RES(RES(:, 3) == 48 | RES(:, 3) == 46, 2), RES(RES(:, 3) == 48 | RES(:, 3) == 46, 5),'VariableNames',{'GeneticID','Age','ClusterSize'});
        end
        lme = fitlme(tbl,'ClusterSize~Age+GeneticID+(Age|GeneticID)')
    end
    if v > 0
        errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', [CO4(k, :) 1], 'LineWidth', 3); hold on;
        if v == 47 || v == 48 || v > 61
            text(14, 16, ['P'' = ' num2str(round(double(lme.Coefficients(2, 6)), 1, 'significant'))], 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            errorbar(14, 15, 3, 'horizontal', 'Color', 'k', 'LineWidth', 2);
            ax = gca;
            ax.Clipping = 'off';
        end
        mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
        N1 = sum(N(N(:, 3) == v, 7));
        [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
        plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', [CO4(k, :) 1], 'LineWidth', 3);
        fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO4(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xlabel('');
        ylabel('');
        if v ~= 61 && v ~= 46
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
        else
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTick', '');
        end
        axis([14 33 0 29.99]);
        CI = coefCI(mdl);
        [Ypred, Yci] = predict(mdl, 23.5);
        text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, [ ...
            '{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
            '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
            '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
            'FontSize', 20, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    else
        axis([14 33 0 29.99]);
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTickLabel', '');
    end
end
set(gca, 'XTick', XT);
set(nexttile(1), 'YTick', YT(2, :), 'XTick', XT);
set(nexttile(5), 'YTick', YT(2, :));
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
print('fig7', '-dpdf', '-r300');

%% Figure S1: A-P single exponentials
pd = 5;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, '1 - Cumulative Distribution Function', 'Fontweight', 'bold', 'FontSize', 24);
for v = [11 12 14 13 15 16]
    TXT = {};
    nt = nexttile;
    for apf = 1:2
        if apf == 1
            nn = N(N(:, 3) == v & N(:, 2) < 24, :);
            txt = ' 15-23 hr APF';
            ft1 = '--';
            ft2 = ':';
            ti = 2;
        else
            nn = N(N(:, 3) == v & N(:, 2) > 23, :);
            txt = ' 24-32 hr APF';
            ft1 = '-';
            ft2 = '-';
            ti = 1;
        end
        stat = zeros(size(nn, 1), 2);
        j = 1;
        for i = nn(:, 1)'
            TXT{end+1} = '';
            y1 = 1-histcounts(RES(RES(:, 1) == i & RES(:, 8) == pd, 5), x, 'Normalization', 'cdf');
            plot(x(1:end-1), y1, ft1, 'Color', [CO6(v-10, :) 1], 'LineWidth', ti); hold on;
            t1 = find(y1==mink(unique(y1), 1), 1);
            G1 = fit(x(FMI:t1)', y1(FMI:t1)', fittype('exp1'));
            stat(j, :) = [G1.a -1/G1.b];
            j = j + 1;
        end
        plot(x, median(stat(:, 1))*exp(-x/median(stat(:, 2))), ft2, 'Color', [CO6(v-10, :) 1], 'LineWidth', ti*4); %
        TXT{end+1} = [txt ', n=' num2str(size(nn, 1)) ', N=' num2str(sum(nn(:, pd+3))) ', ' num2str(round(median(stat(:, 1)), 1)) 'e^{-x/\bf{' num2str(round(median(stat(:, 2)))) '}}']; %
    end
    xlabel('');
    ylabel('');
    if mod(v, 10) == 1 || mod(v, 10) == 3
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 2, 'XTick', 0:20:40, 'YScale', 'log');
    else
        set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 2, 'XTick', 0:20:40, 'YScale', 'log', 'YTickLabel', '');
    end
    legend(TXT, 'Box', 'off', 'FontSize', 21);
    axis([0 60 YMI 2]);
end
set(gca, 'XTick', 0:20:100);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTick', 20:20:40);
set(nexttile(4), 'YTick', [1e-3 1e-2 0.1 1]);
saveas(gcf, 'figS1.png');

%% Figure S2: A-P grow with age
pd = 5;
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 1100]);
tl = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Average cluster size', 'Fontweight', 'bold', 'FontSize', 24);
for v = [11 12 14 13 15 16]
    nexttile;
    errorbar(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), Nle(Nle(:, 3) == v, pd+3), '.', 'Color', CO6(mod(v, 10), :), 'LineWidth', 3); hold on;
    mdl = fitlm(Nl(Nl(:, 3) == v, 2), Nl(Nl(:, 3) == v, pd+3), 'Weights', Nle(Nle(:, 3) == v, pd+3));
    [Ypred, Yci] = predict(mdl,Nl(Nl(:, 3) == v, 2));
    plot(Nl(Nl(:, 3) == v, 2), Ypred, '-', 'Color', [CO6(mod(v, 10), :) 0.8], 'LineWidth', 3);
    fill([Nl(Nl(:, 3) == v, 2)', fliplr(Nl(Nl(:, 3) == v, 2)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    axis([14 33 0 15]);
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    text(15, 13.5, ['{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
        '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
        '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], 'FontSize', 18);
    set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', XT, 'YTick', 0:2:14, 'YTickLabel', '');
    xlabel('');
    ylabel('');
end
set(gca, 'XTick', XT);
set(nexttile(1), 'YTick', 0:2:14, 'YTickLabel', {'0', '2', '4', '6', '8', '10', '', ''});
set(nexttile(4), 'YTick', 0:2:14, 'YTickLabel', {'0', '2', '4', '6', '8', '10', '', ''});
saveas(gcf, 'figS2.png');

%% Figure S3: Density
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
xlabel(tl, 'Pupal age [hours APF]', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, 'Density [clusters/\mum^2]', 'Fontweight', 'bold', 'FontSize', 24);
for v = [11 12 14 13 15 16]
    nt = nexttile;
    RESs = RES(RES(:, 3) == v & RES(:, 8) > pd, :);
    j = 1;
    sta = zeros(numel(unique(RESs(:, 1))), 3);
    for i = unique(RESs(:, 1))'
        errorbar(unique(RESs(RESs(:, 1) == i, 2))+(rand()/3-1/6), mean(RESs(RESs(:, 1) == i, 10)), std(RESs(RESs(:, 1) == i, 10)), 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3); hold on;
        sta(j, :) = [unique(RESs(RESs(:, 1) == i, 2)) mean(RESs(RESs(:, 1) == i, 10)) numel(RESs(RESs(:, 1) == i, 10))];
        j = j + 1;
    end
    mdl = fitlm(sta(:, 1), sta(:, 2), 'Weights', sta(:, 3));
    [Ypred, Yci] = predict(mdl, sta(:, 1));
    plot(sta(:, 1), Ypred, '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
    fill([sta(:, 1)', fliplr(sta(:, 1)')], [Yci(:, 2)' fliplr(Yci(:, 1)')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    xlabel('');
    ylabel('');
    set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1);
    axis([14 33 0 9]);
    CI = coefCI(mdl);
    [Ypred, Yci] = predict(mdl, 23.5);
    text(min(xlim)+diff(xlim)/30, max(ylim)-diff(ylim)/30, ['{\bf Mean = ' num2str(round(Ypred)) '}, 95%CI [' num2str(floor(Yci(1))) ', ' num2str(ceil(Yci(2))) ']' newline ...
        '{\bf Slope = ' num2str(round(table2array(mdl.Coefficients(2, 1)), 1, 'significant')) '}, 95%CI [' num2str(round(CI(2), 1, 'significant')) ', ' num2str(round(CI(4), 1, 'significant')) ']' newline ...
        '{\bf P-value = ' num2str(round(table2array(mdl.Coefficients(2, 4)), 1, 'significant')) '}, n=' num2str(mdl.NumObservations) ', N=' num2str(sum(N(N(:, 3) == v, pd+3))) ], ...
        'FontSize', 20, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
set(nexttile(1), 'YTick', 0:2:8, 'XTick', '');
set(nexttile(2), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(3), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(4), 'YTick', 0:2:8);
set(nexttile(5), 'YTick', 0:2:8, 'YTickLabel', '');
set(nexttile(6), 'YTick', 0:2:8, 'YTickLabel', '');
saveas(gcf, 'figS3.png');

%% Figure S4: Mobility analyze
rng(1);
DX = 3;
BACK = 2;
DT = 20;

STAT = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 1:456
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
save('STAT-figS4.mat', 'STAT');

%% Figure S4: Mobility visualize
figure(1); clf;
load('STAT-figS4');
TH = 2;
set(gcf, 'Color', 'w', 'Position', [0 0 900 900]);
tl = tiledlayout(4, 2, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 20);
ylabel(tl, 'Log diffusion coefficient [\mum^2/s]', 'Fontweight', 'bold', 'FontSize', 20);
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
        scatter(x-0.01, y, 25, [CO6(mod(v, 10), :)], 'filled', 'MarkerFaceAlpha', 0.3); hold on;
        plot(sta(:, 1), sta(:, 2), '-', 'Color', [CO6(mod(v, 10), :) 1], 'LineWidth', 3);
        fill([sta(:, 1)', fliplr(sta(:, 1)')], [(sta(:, 2)+sta(:, 3))' fliplr((sta(:, 2)-sta(:, 3))')], CO6(mod(v, 10), :), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xlabel('');
        ylabel('');
        set(gca, 'Box', 'on', 'FontSize', 20, 'LineWidth', 1, 'XTick', 0:5:20, 'YTick', -13:2:-5, 'YTickLabel', '');
        axis([0 25 -14 -3]);
        text(max(xlim)-diff(xlim)/30, max(ylim)-diff(ylim)/20, ['n=' num2str(size(nn, 1)) ', N=' num2str(size(x, 1))], 'FontSize', 19, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        if mod(tilenum(nt), 2) == 0
            text(max(xlim)-diff(xlim)/30, max(ylim)-diff(ylim)/6, ['P'' = ' mat2str(round(double(lme.Coefficients(2, 6)), 1, 'significant'))], 'FontSize', 19, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        end
    end
end
set(gca, 'XTick', 0:5:30);
set(nexttile(1), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(3), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(5), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5}, 'XTick', '');
set(nexttile(7), 'YTick', -13:2:-5, 'YTickLabel', {-13, '', -9, '', -5});
saveas(gcf, 'figS4.png');

%% Figure S5: P or D single exponentials
V = [2 900];
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [0 0 1600 900]);
tl = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'tight');
xlabel(tl, 'Number of molecules in clusters', 'Fontweight', 'bold', 'FontSize', 24);
ylabel(tl, '1 - Cumulative Distribution Function', 'Fontweight', 'bold', 'FontSize', 24);
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
                ti = 2;
            else
                nn = N(N(:, 3) == v & N(:, 2) > 23, :);
                txt = ' 24-32 hr APF';
                ft1 = '-';
                ft2 = '-';
                ti = 1;
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
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', 0:20:40, 'YScale', 'log', 'YTickLabel', '');
        else
            set(gca, 'Box', 'on', 'FontSize', 24, 'LineWidth', 1, 'XTick', 0:20:40, 'YScale', 'log');
        end
        legend(TXT, 'Box', 'off', 'FontSize', 19);
        axis([0 60 YMI 5]);
    end
end
set(gca, 'XTick', 0:20:60);
set(nexttile(1), 'YTick', [1e-3 1e-2 0.1 1], 'XTick', 20:20:40);
set(nexttile(5), 'YTick', [1e-3 1e-2 0.1 1]);
saveas(gcf, 'figS5.png');

%% Figure S7: 2-color A or P
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [150 0 200 600]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'none');
figureS7(73, 'Vang', 'Fmi', 2);
figureS7(71, 'Fz', 'Vang', 1.5);
figureS7(72, 'Vang', 'Pk', 1.5);
print('figS7.png', '-dpng','-r0');

function figureS7(v, xlab, ylab, ymax)
nexttile;
load('STAT-fig5.mat');
n = numel(unique(STAT(STAT(:, 2) == v & STAT(:, 3) == 5, 1)));
n1 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == 5, 4))/1e4;
n2 = single(STAT(STAT(:, 2) == v & STAT(:, 3) == 5, 5))/1e4;
scatter(n1, n2, 20, 'filled', 'MarkerFaceAlpha', 0.3);
R = corrcoef(n1, n2);
set(gca, 'Box', 'on', 'FontSize', 17, 'LineWidth', 1, 'XTick', 0:0.2:0.4);
text(0.5/2, ymax*1.03, ['n=' num2str(n) ', N=' num2str(numel(n1)) ', r=' num2str(R(2, 1), 1)], 'FontSize', 14, 'HorizontalAlignment', 'center');
axis([0 0.5 0 ymax*1.1], 'square');
xlabel(xlab, 'FontWeight', 'bold', 'Color', '#EE220C');
ylabel(ylab, 'FontWeight', 'bold', 'Color', '#1DB100');
end

%% Figure S11: Pentamer
clear; rng(1); tic;
DX = 3;
DT = 4;
BACK = 2;
BINW = 1;
FMIN = 80;
FMAX = 160;

RES = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 488:490
    F = FILE(n).name(1:end-4); disp(F);
    CH = str2double(F(22));
    for ch = 1:CH
        RAW = single(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [ch CH inf]}));
        if size(RAW, 3) == 1
            RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
            RAW = single(RAW(:, :, ch:CH:end));
        end
        SMOO = movmedian(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), [DT 0], 3);
        MASK = (blockproc(SMOO(:, :, 1), [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2]));
        [r3, c3] = find(MASK); % puncta centers
        smoo = reshape(permute(SMOO, [3 1 2]), size(RAW, 3), []);
        smor = smoo(:, (MASK)); % puncta centers
        SMOO = SMOO(:, :, 1);
        res = nan(size(r3, 1), 9);
        tank = logical(triu(ones(size(smoo, 1), size(smoo, 1)), 1));
        parfor i = 1:size(smor, 2)
            wi = smor(:, i) - smor(:, i)';
            [pwr, f] = pspectrum(histcounts(wi(tank), 'BinWidth', BINW, 'Normalization', 'probability'), 'FrequencyLimits', [BINW*2*pi/FMAX BINW*2*pi/FMIN]);
            I = find(islocalmax(pwr, 'MaxNumExtrema', 1));
            if ~isempty(I)
                res(i, :) = [n str2double(F(17:18)) str2double(F(2:3)) BINW*2*pi/f(I) smor(1, i)/(BINW*2*pi/f(I)) c3(i) r3(i) 4 ch];
            end
        end
        res(isnan(res(:, 4)), :) = [];
        RES = [RES; res];
    end
end

figure(1); clf;
set(gcf, 'color', 'w')
histogram(RES(:, 5), -0.5:1:15.5)
title(['N = ' num2str(size(RES, 1)) '; \mu = ' num2str(round(mean(RES(:, 5)), 1)) ' +- ' num2str(round(std(RES(:, 5)), 1))]);
xlabel('Cluster size', 'FontWeight', 'bold');
ylabel('Number of clusters', 'FontWeight', 'bold');
set(gca, 'LineWidth', 2, 'FontSize', 20);
saveas(gcf, 'figS11.png')

%% SAVE step size (11)
function stepSize(F, i, smor, edges, pddf, f, pwr, ss, ch)
figure(1); clf;
tiledlayout(1, 3);
set(gcf, 'color', 'w', 'Position', [0 0 1800 600]);

nexttile;
%plot(squeeze(valr(:, i)), 'color', [0 0 1 0.15], 'LineWidth', 0.1); hold on;
plot([0 size(smor, 1)], [0 0], 'color', [0 0 0 0.5], 'LineWidth', 3); hold on;
plot(squeeze(smor(:, i)), 'color', [0, 0.4470, 0.7410 1], 'LineWidth', 4);
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:20:1000, 'YTick', 0:100:500);
xlabel('Frames', 'Fontweight', 'bold');
ylabel('Intensity', 'Fontweight', 'bold');
axis([0 100 -200 500]);

nexttile;
plot(edges(1:end-1), pddf, 'LineWidth', 6);
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:200:800, 'YTick', [1e-6 1e-5 1e-4 1e-3 1e-2], 'YScale', 'log');
xlabel('Pairwise difference', 'Fontweight', 'bold');
ylabel('Probability', 'Fontweight', 'bold');
title(['Channel ' num2str(ch)])
axis([-200 800 1e-6 0.0101]);

nexttile;
INTS = 2*(edges(2)-edges(1))*pi./f;
plot(INTS, pwr, 'LineWidth', 6); hold on;
scatter(round(ss), max(pwr(round(INTS) == round(ss))), 400, 'r', 'filled');
set(gca, 'FontSize', 20, 'LineWidth', 3, 'XTick', 0:40:200);
xlabel('Intensity', 'Fontweight', 'bold');
ylabel('Power', 'Fontweight', 'bold');
axis([0 160 0 1.5e-7], 'auto y');

mkdir(['../../6-StepSize/' F '-' num2str(ch) '/']);
warning('off', 'MATLAB:MKDIR:DirectoryExists');
saveas(gcf, ['../../6-StepSize/' F '-' num2str(ch) '/' num2str(i) '.png']);
end