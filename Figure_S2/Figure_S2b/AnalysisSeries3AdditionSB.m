function [Objects] = AnalysisSeries3AdditionSB(InfoTableThis, PreviewPath, ch1, ch2, ch3, ch4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%    InfoTableThis

    % Segment nuclei
    NucChannelDoG = imfilter(ch1, fspecial('gaussian', 10, 2), 'symmetric')...
        - imfilter(ch1, fspecial('gaussian', 60, 20), 'symmetric'); % it(max(NucChannelCorr, [], 3))
    NucMask = NucChannelDoG > 50; % vol(NucChannelDoG) % vol(ch1)
    %vol(NucMask)


    % Remove small objects which are certainly no nuclei

    NucMask = bwareaopen(NucMask, 200);
    if sum(NucMask(:)) < 15000
        Objects = table();
        return
    end
    ch1AVforNuc = imfilter(ch1, fspecial('average', 5), 'symmetric'); %vol( ch1AVforNuc)
    NucMaskHigh = (ch1AVforNuc > 1000) .* NucMask; % vol(NucMaskHigh)
    NucMaskLow = NucMask .* ~NucMaskHigh; % vol(NucMaskLow)
    NucClassifier = NucMask + 4.* NucMaskLow + 8.* NucMaskHigh; % vol(NucClassifier)

  

    %% 3D strel
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));

    
    % Segment neurons
    NeuroSegmentationIm = ch4; % Tuj1 only % vol(ch4,0,50)
    %vol(NeuroSegmentationIm, 0,50)
    NeuroSegmentationImGlobal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 10, 3), 'symmetric'); % vol(NeuroSegmentationImGlobal,0,50)
    NeuroSegmentationImLocal = imfilter(NeuroSegmentationIm, fspecial('gaussian', 10, 3), 'symmetric')...
        - imfilter(NeuroSegmentationIm, fspecial('gaussian', 20, 6), 'symmetric'); % vol(NeuroSegmentationImLocal,0,20)
    NeuroMaskGlobal = NeuroSegmentationImGlobal > 7; % vol(NeuroSegmentationImGlobal) % vol(NeuroMaskGlobal)
    NeuroMaskLocal = NeuroSegmentationImLocal > 3; % vol(NeuroSegmentationImLocal) % vol(NeuroMaskLocal)
    NeuroMask = NeuroMaskGlobal | NeuroMaskLocal; % Fill gaps % vol(NeuroMask)

    % Remove small objects which are certainly no neurons
    NeuroMask = bwareaopen(NeuroMask, 200);  


    %% TH Mask

    ch2LP = imfilter(ch2, fspecial('gaussian', 10, 1), 'symmetric');
    %vol(ch2LP,0,500)
    TH_Mask = ch2LP > 300;
    %vol(TH_Mask)
    ch2DoG =  imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 21, 7), 'symmetric');
    %vol(ch2DoG, 0, 100, hot)
    TH_LocalMask = ch2DoG > 30;
    %vol(TH_LocalMask, 0, 1)
    TH_Mask = TH_Mask | TH_LocalMask;
    TH_Mask = TH_Mask & ~NucMask; %substracting nuc area & TH overlapping 
    TH_Mask = bwareaopen(TH_Mask, 100); %vol(TH_Mask)
    %it(max(TH_Mask,[],3))
    

    %% Ki67
    %vol(ch3, 0,100)
    ch3LP = imfilter(ch3, fspecial('gaussian', 10, 1), 'symmetric');
    Ki67Mask = ch3LP > 50;
    Ki67Mask = bwareaopen(Ki67Mask, 50);
    %vol(Ki67Mask)
    
    %% Ki67 in TH
    %Ki67inTH = imreconstruct(TH_Mask, Ki67Mask) .*  Ki67Mask;
    
    %% Non Dopaminergic neurons
    nonDANeurons = NeuroMask .* ~TH_Mask;
    %vol(nonDANeurons)
    %% Collect data

    % General Opera type features
    Objects = table();
    InfoTableThisSample = InfoTableThis; 
    Objects.NucVol = sum(NucMask(:));
    Objects.NeuroMask = sum(NeuroMask(:));
    Objects.NeuroMaskbyNucVol = sum(NeuroMask(:)) / sum(NucMask(:));
    Objects.NucHighProportion = sum(NucMaskHigh(:)) / sum(NucMask(:));
    Objects.NucLowProportion = sum(NucMaskLow(:)) / sum(NucMask(:));
    Objects.HoechstInNucMask = sum(sum(sum(ch1 .* uint16(NucMask)))); 
    Objects.THbyNucVol = sum(ch3(:)) / sum(NucMask(:));
    Objects.THbyNeuroVol = sum(ch3(:)) / sum(NeuroMask(:));
    Objects.Ki67byNucVol = sum(sum(sum(uint16(Ki67Mask) .* ch3))) / sum(NucMask(:)); % Remove background and normalize by nuclei
    Objects.Ki67Average = sum(sum(sum(uint16(Ki67Mask) .* ch3))) / sum(Ki67Mask(:)); % Remove background and normalize by Ki67Mask
    Objects.nonDANeuronsbyNucVol = sum(nonDANeurons(:)) / sum(NucMask(:));
    Objects.nonDANeurons = sum(nonDANeurons(:));
    Objects = [InfoTableThisSample, Objects];

  


%     %% save 2D previews
%     
%     PreviewPath = [DataPath, filesep, 'Previews'];
%     Barcode = InfoTableThis.Barcode{:};
%     row = InfoTableThis.Row;
%     column = InfoTableThis.Column;
%     field = InfoTableThis.field;
%     
%     % Scalebar
%     imSize = [size(ch1, 1), size(ch1, 2)];
%     [BarMask, BarCenter] = f_barMask(50, 0.646, imSize, imSize(1)-50, 50, 10);
%     %it(BarMask)
% 
%     % Ki671 (ch3)
%     PreviewKi67 = imoverlay2(imadjust(max(ch3,[],3)), bwperim(max(NucMask,[],3)), [0 0 1]);
%     PreviewKi67 = imoverlay2(PreviewKi67, bwperim(max(Ki67Mask,[],3)), [1 0 0]);
%     PreviewKi67 = imoverlay2(PreviewKi67, BarMask, [1 1 1]);
%     %it(PreviewKi67)
%     
%     % TH (ch2)
%     PreviewTH = imoverlay2(imadjust(max(ch2,[],3)), bwperim(max(TH_Mask,[],3)), [1 0 0]);
%     PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
%     %it(PreviewTH)
%    
%     filename_Ki67 = sprintf('%s\\%s__%03d%03d%03d_Ki67.png', PreviewPath, Barcode, row, column, field); 
%     filename_TH = sprintf('%s\\%s__%03d%03d%03d_TH.png', PreviewPath, Barcode, row, column, field); 
%     
%     imwrite(PreviewKi67, filename_Ki67);
%     imwrite(PreviewTH, filename_TH);
    

end

