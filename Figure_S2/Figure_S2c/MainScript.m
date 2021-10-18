%% Clear Matlab workspace
clear
clc

%% User inputs per run

% New aquisition no overlap one plane
Barcodes = {
% 'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_2_N1_JW403\JW_20170207_60X_DID4_2_N1_JW403'};
% 'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_N1_JW403\JW_20170207_60X_DID4_N1_JW403'};
%   'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_N2_JW403\JW_20170207_60X_DID4_N2_JW403'};
%   'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_N3_JW403\JW_20170207_60X_DID4_N3_JW403'};
%   'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N1_JW403\JW_20170207_60X_DID10_N1_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N2_JW403\JW_20170207_60X_DID10_N2_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N3_JW403\JW_20170207_60X_DID10_N3_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N1_JW403\JW_20170207_60X_DID14_N1_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N2_JW403\JW_20170207_60X_DID14_N2_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N3_JW403\JW_20170207_60X_DID14_N3_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID4_N4_JW403\JW_20170225_60X_DID4_N4_JW403'};
  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID10_N4_JW403\JW_20170225_60X_DID10_N4_JW403'};
%  'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID14_N4_JW403\JW_20170225_60X_DID14_N4_JW403'};

% Thresholds adjusted in the image analysis based on low levels of TUJ1
for b = 1:size(Barcodes, 2) %loop over barcodes
    Barcode = Barcodes{b};
    Objects = {};   
    clear ObjectsAll
    ObjectsCell = {};
    
    %% Document script
    BarcodeName = regexp(Barcode, '.*\\(.*)', 'tokens'); 
    BarcodeName = BarcodeName{:};
    SavePathBegining = {'S:\HCS_Platform\Data\SilviaBolognin\Jonas_Script_CellReportPaper\'};
    SavePathMainDirectory = [SavePathBegining{:}, BarcodeName{:}];
    AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
    SavePath = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
    mkdir(SavePath)
    FileNameShort = mfilename;
    newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
    FileNameAndLocation = [mfilename('fullpath')];
    currentfile = strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup); 
    PreviewPath = [SavePath, filesep, 'Previews'];
    mkdir(PreviewPath)
    Version = version();
    save([SavePath filesep 'MatlabVersion.mat'], 'Version')

    f_LogDependencies(FileNameShort, SavePath); 

    %% Load data
%     [files, info, areaNames, plateLayout, areaMap] = f_Local_Read_Files_Timeseries_Mosaic_BT2('S:\OperaQEHS\OperaArchiveCol', Barcode, 1);
%     InfoTable = struct2table(info);
        InfoTable = f_InfoTable(Barcodes{b});

    for f=1:height(InfoTable)%loop over fields
    % for f=1:height(InfoTable)%loop over fields
        InfoTableThis = InfoTable(f,:);
        cube = readflexcube(InfoTable.files{f}, 'PlaneCount', 4); % Read 4-D image cube
        ch1 = cube.data(:,:,1); %TH % it(ch1(:,:,3)) vol(ch1) for planes:ch1 = cube.data(:,:,:,1);
        ch2 = cube.data(:,:,2); %Hoechst % it(ch2(:,:,3)) vol(ch2)
        ch3 = cube.data(:,:,3); %Map2 % it(ch3(:,:,3)) vol(ch3)
        ch4 = cube.data(:,:,4); %Tuj1 % it(ch4(:,:,3)) vol(ch4)
       
        Objects = AnalysisSeries3AdditionSB(InfoTableThis, PreviewPath, ch1, ch2, ch3, ch4);
        
        if size(Objects,1) == 1
           ObjectsCell{f} = Objects;
        end
            
    end
    
%      if exist('ObjectsCell')
%         ObjectsBarcodes{b} = vertcat(ObjectsCell{:});
%      end
    
      ObjectsAll = vertcat(ObjectsCell{:});
      save([SavePath, filesep, 'data.mat'], 'ObjectsAll');
      writetable(ObjectsAll, [SavePath, '\Objects.csv'], 'WriteRowNames', true)  
      writetable(ObjectsAll, [SavePath, '\Objects.xls'], 'WriteRowNames', true)

end

