%% REGISTER channels. Script #1 used for registering 2-color images only.
clear; tic;
DX = 5;
BACK = 5;
DT = 5;
FILE = dir('../../1-TIF/c*.tif');
for n = 1
    F = FILE(n).name(1:end-4);
    load(['../../3-MAP/' F '.mat']);
    ch1 = register(F, 1, DT, BACK, DX);
    ch2 = register(F, 2, DT, BACK, DX);
    ch1((blockproc(ch1, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP > 0) == 0) = 0;
    ch2((blockproc(ch2, [1 1], @(b)(max(b.data, [], 'all') == b.data), BorderSize=[2 2], UseParallel=true) & MAP > 0) == 0) = 0;
    [r1, c1] = find(ch1);
    [r2, c2] = find(ch2);
    Idx = knnsearch([r2 c2], [r1 c1]);
    dy = fix(median(r2(Idx) - r1)/2);
    dx = fix(median(c2(Idx) - c1)/2);
    disp([F ': N=' num2str(numel(Idx)) '; dx=' num2str(dx) '; dy=' num2str(dy)]);
    CH1all = circshift(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 1 inf], [1 1 inf], [1 2 inf]}), [dy dx 0]);
    imwrite(CH1all(:, :, 1), 'myFile.TIFF', 'Compression', 'none');
     for i = 2:size(CH1all, 3)
         imwrite(CH1all(:, :, i), 'myFile.TIFF', 'writemode', 'append', 'Compression', 'none');
     end
    imwrite(CH1all(:, :, 2001), 'myFile1.TIFF', 'Compression', 'none');
    for i = 2002:size(CH1all, 3)
        imwrite(CH1all(:, :, i), 'myFile1.TIFF', 'writemode', 'append', 'Compression', 'none');
    end
end

function I = register(F, ch, DT, BACK, DX)
RAW = single(tiffreadVolume(['../../1-TIF/' F '.tif'], 'PixelRegion', {[1 inf], [1 inf], [ch 2 2*DT]}));
if size(RAW, 3) == 1
    RAW = tiffreadVolume(['../../1-TIF/' F '.tif']);
    RAW = single(RAW(:, :, ch:2:2*DT));
end
I = uint32(sum(imboxfilt(RAW - imgaussfilt(RAW, BACK, 'padding', 'symmetric'), DX, 'NormalizationFactor', 1), 3));
end