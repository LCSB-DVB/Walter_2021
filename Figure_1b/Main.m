%% Clear Matlab workspace
clear
clc

if strcmp(getenv('COMPUTERNAME'), 'LCSB-HCS2')
    delete(gcp('nocreate'))
    parpool(28)
end

%% User inputs per run

% Setup was 20x bin2 for all
Barcodes = {'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170411_20x_DID10_rescue_N3_jw4xx\JW_20170411_20x_DID10_rescue_N3_jw4xx';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID14_N4_JW403\JW_20170225_60X_DID14_N4_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID10_N4_JW403\JW_20170225_60X_DID10_N4_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170225_60X_DID4_N4_JW403\JW_20170225_60X_DID4_N4_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N3_JW403\JW_20170207_60X_DID14_N3_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N3_JW403\JW_20170207_60X_DID10_N3_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_N3_JW403\JW_20170207_60X_DID4_N3_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_N2_JW403\JW_20170207_60X_DID4_N2_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N2_JW403\JW_20170207_60X_DID10_N2_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N2_JW403\JW_20170207_60X_DID14_N2_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID14_N1_JW403\JW_20170207_60X_DID14_N1_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID10_N1_JW403\JW_20170207_60X_DID10_N1_JW403';...
            'S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_2_N1_JW403\JW_20170207_60X_DID4_2_N1_JW403'};
Channels = 4;
SeriesDefinitions = {'Sox1','PARP','Ki67'};

%% Document script
SavePathMainDirectory = 'S:\HCS_Platform\Data\JonasWalter\Differentiation_2D';
AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
SavePath = [SavePathMainDirectory, '\Analysis_', AnalysisTimeStamp];
mkdir(SavePath)
FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
FileNameAndLocation = mfilename('fullpath');
currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

VersionMatlab = version;
save([SavePath, '\', 'MatlabVersion.mat'], 'VersionMatlab');
f_LogDependencies(FileNameShort, SavePath); 


%%%%%%%%%%%%%%%%% Here the barcode loop starts
for n = 1:size(Barcodes,1)

    Barcode = Barcodes{n};
    Objects = {};

    %% Document Matlab state
    Version = ver;
    MegatronPath = pwd;

    %% Prepare output directories
    S1Path = [SavePath, '\SeriesSox1'];
    mkdir(S1Path)
    mkdir([S1Path filesep 'Previews'])
    S2Path = [SavePath, '\SeriesPARP'];
    mkdir(S2Path)
    mkdir([S2Path filesep 'Previews'])
    S3Path = [SavePath, '\SeriesKI67'];
    mkdir(S3Path)
    mkdir([S3Path filesep 'Previews'])

    %% Image analysis

    %InfoTable = f_InfoTable('S:\OperaQEHS\OperaDB\Jonas Walter\JW_20170207_60X_DID4_2_N1_JW403');
    InfoTable = f_InfoTable(Barcodes{n});

    for i=1:height(InfoTable)
        cube = readflexcube(InfoTable.files{i}, 'PlaneCount', 1); % Read 4-D image cube

        ch1 = cube.data(:, :, 1, 1:Channels:size(cube.data, 4)); % 405 Hoechst
        ch1 = squeeze(ch1);
        %vol(ch1)

        ch2 = cube.data(:, :, 1, 2:Channels:size(cube.data, 4)); % 488 TH
        ch2 = squeeze(ch2);
        %vol(ch2)

        ch3 = cube.data(:, :, 1, 3:Channels:size(cube.data, 4)); % 568 Sox1 (stemmnes) or PARP (picnotic/late apoptosis) or Ki67 (cell cycle exit)
        ch3 = squeeze(ch3);
        %vol(ch3)

        ch4 = cube.data(:, :, 1, 4:Channels:size(cube.data, 4)); % Tuj1 only analysed with S1 and S2
        ch4 = squeeze(ch4);
        %vol(ch4,0,100)

        InfoTableThis = InfoTable(i,:);
        AreaName = InfoTableThis.AreaName{:};
        Series = regexp(AreaName, '.*\_S(\d).*', 'tokens');
        try
            Series = Series{:}{:};
        catch
            Series = '1'; % it is S1
        end
        % Series = AreaName(end);
        
        if str2num(Series) == 1
            ObjectsS1 = AnalysisSeries1(InfoTableThis, i, S1Path, ch1, ch2, ch3);
            if size(ObjectsS1,1) == 1
                ObjectsS1Cell{i} = ObjectsS1;
            end
        elseif str2num(Series) == 2
            ObjectsS2 = AnalysisSeries2(InfoTableThis, i, S2Path, ch1, ch2, ch3, ch4);
            if size(ObjectsS2,1) == 1
                ObjectsS2Cell{i} = ObjectsS2;
            end
        elseif str2num(Series) == 3
            ObjectsS3 = AnalysisSeries3(InfoTableThis, i, S3Path, ch1, ch2, ch3, ch4);
            if size(ObjectsS3,1) == 1
                ObjectsS3Cell{i} = ObjectsS3;
            end
        end
        
    end
    
    if exist('ObjectsS1Cell')
        ObjectsS1BarCodes{n} = vertcat(ObjectsS1Cell{:});
    end
    
    if exist('ObjectsS2Cell')
        ObjectsS2BarCodes{n} = vertcat(ObjectsS2Cell{:});
    end
    
    if exist('ObjectsS3Cell')
        ObjectsS3BarCodes{n} = vertcat(ObjectsS3Cell{:});
    end
    
end

if exist('ObjectsS1BarCodes')
    ObjectsAllS1 = vertcat(ObjectsS1BarCodes{:});
    save([S1Path, filesep, 'data.mat'], 'ObjectsAllS1');
    writetable(ObjectsAllS1, [S1Path, '\Objects.csv'], 'WriteRowNames', true)  
    writetable(ObjectsAllS1, [S1Path, '\Objects.xls'], 'WriteRowNames', true)
end

if exist('ObjectsS2BarCodes')
    ObjectsAllS2 = vertcat(ObjectsS2BarCodes{:});
    save([S2Path, filesep, 'data.mat'], 'ObjectsAllS2');
    writetable(ObjectsAllS2, [S2Path, '\Objects.csv'], 'WriteRowNames', true)  
    writetable(ObjectsAllS2, [S2Path, '\Objects.xls'], 'WriteRowNames', true)
end

if exist('ObjectsS3BarCodes')
    ObjectsAllS3 = vertcat(ObjectsS3BarCodes{:});
    save([S3Path, filesep, 'data.mat'], 'ObjectsAllS3');
    writetable(ObjectsAllS3, [S3Path, '\Objects.csv'], 'WriteRowNames', true)  
    writetable(ObjectsAllS3, [S3Path, '\Objects.xls'], 'WriteRowNames', true)
end
