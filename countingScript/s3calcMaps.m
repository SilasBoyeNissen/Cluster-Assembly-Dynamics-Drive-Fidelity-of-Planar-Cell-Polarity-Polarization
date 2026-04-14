%% CALC maps. Script #3 used to convert ROIs from Script #2 into MATLAB maps for each image file.
close all; clear;
W = 9;
delete('../../3-MAP/*')
FILE = dir('../../1-TIF/c*.tif');
for n = 1:size(FILE, 1)
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