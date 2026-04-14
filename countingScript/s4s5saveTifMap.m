%% SAVE tif+map. Script #4-5 saves the TIF files' first frame and processed overlaying the map for quick look-up.
close all; clear;
DX = 5; % 3
BACK = 5; % 2
DT = 5;
FR = 1;
TP = 0.75; % 1

delete('../../4-tif/*')
delete('../../5-map/*')
FILE = dir('../../1-TIF/c*.tif');
for n = 1:size(FILE, 1)
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
        I(10, 10:41) = max(I, [], 'all');
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