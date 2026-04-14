%% DRAW polylines. Script #2 used to mark proximal, distal, anterior-posterior boundaries.
close all; clear;
ROI = [];
FILE = dir('../../1-TIF/c*.tif');
for n = 1
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