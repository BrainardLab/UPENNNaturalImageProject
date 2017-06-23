% RunAllCalibrationsAndTests.m 
%
% 08/23/2010  gt    Added a flag for removing and regenerating old PGM/rawmat files
% 11/10/2010  dhb   Mostly cosmetic.  On reprocess, delete .ppm and .pgm as well as .PPM and .PGM
%             dhb   More control options.
% 1/5/11      dhb   Set up for full run through.  Will see if it works without crashing ...
% 1/6/11      dhb   Fixed minor bug, MakeCamCal after cd .. rather than before.
% 6/18/17     dhb   Happy Father's Day!

%% Clear and set paths
clear; close all;

%% Select what to do
DODARK = 1;
DOAPERTURE = 1;
DOLINEARITY = 1;
DOISO = 1;
DOCOMPARISON = 1;
DOSPECTRALSENS = 1;
DOMCC = 1;
DOMTF = 1;
DOLUMCHECK = 1;

DELETEONLY = 0;
DELETEOLD = 1;
REPROCESS = 1;

idbPath = getpref('UPENNNaturalImageProject','calibrationImageDir');
dcrawPath = fullfile(fileparts(mfilename('fullpath')),'TESTDCRAW');

%% Start where we live
cd(fullfile(fileparts(mfilename('fullpath'))));

%% If we want to reprocess from scratch, can optionally delete all the old files.
% But sometimes don't want to do that, if we're running only one aspect of this and
% have recently deleted everything.
if (REPROCESS == 1 && DELETEOLD)
    unix(['rm -f ' idbPath '/*/DSC*.mat']);
    unix(['rm -f ' idbPath '/*/DSC*.txt']);
    unix(['rm -f ' idbPath '/*/*.PPM']);
    unix(['rm -f ' idbPath '/*/*.PGM']);
    unix(['rm -f ' idbPath '/*/*.ppm']);
    unix(['rm -f ' idbPath '/*/*.pgm']);

    unix(['rm -f ' idbPath '/*/*/DSC*.mat']);
    unix(['rm -f ' idbPath '/*/*/DSC*.txt']);
    unix(['rm -f ' idbPath '/*/*/*.PPM']);
    unix(['rm -f ' idbPath '/*/*/*.PGM']);
    unix(['rm -f ' idbPath '/*/*/*.ppm']);
    unix(['rm -f ' idbPath '/*/*/*.pgm']);

    unix(['rm -f ' idbPath '/*/*/*/DSC*.mat']);
    unix(['rm -f ' idbPath '/*/*/*/DSC*.txt']);
    unix(['rm -f ' idbPath '/*/*/*/*.PPM']);
    unix(['rm -f ' idbPath '/*/*/*/*.PGM']);
    unix(['rm -f ' idbPath '/*/*/*/*.ppm']);
    unix(['rm -f ' idbPath '/*/*/*/*.pgm']);
end

if (DELETEONLY == 1)
    return;
end

%% DARK CURRENT
if (DODARK)
    cd DARK_IMAGES
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/DARK_IMAGES'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/DARK_IMAGES'], 0);
        ProcessNEFToPGM(1, [idbPath '/DARK_IMAGES_ISOVARY'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/DARK_IMAGES_ISOVARY'], 0);
    end
    EstimateDarkCurrent(idbPath);
    cd ..
end

%% APERTURE
if (DOAPERTURE)
    cd APERTURE
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/APERTURE'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/APERTURE'], 0);
    end
    EstimateApertureLinearity(idbPath)
    cd ..
end

%% EXPOSURE LINEARITY
if (DOLINEARITY)
    cd EXPOSURE_LINEARITY
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/EXPOSURE_LINEARITY'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/EXPOSURE_LINEARITY'], 0);
    end
    EstimateExposureLinearity(idbPath)
    cd ..
end

%% ISO LINEARITY
if (DOISO)
    cd ISO_LINEARITY
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/ISO_LINEARITY_EXPOSURE1'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/ISO_LINEARITY_EXPOSURE1'], 0);
        ProcessNEFToPGM(1, [idbPath '/ISO_LINEARITY_EXPOSURE2'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/ISO_LINEARITY_EXPOSURE2'], 0);
    end
    EstimateISOLinearity(idbPath)
    cd ..
end

%% CAMERA COMPARISON
if (DOCOMPARISON)
    cd CAMERACOMPARISON
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos0'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos0'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos0'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos0'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos+10'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos+10'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos+10'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos+10'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos-10'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/BotsCam/Photos-10'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos-10'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Photos-10'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/BotsCam/Color'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/BotsCam/Color'], 0);
        ProcessNEFToPGM(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Color'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/CAMERACOMPARISON/PhillyCam/Color'], 0);
    end
    CameraComparison(idbPath,'Photos','Photos');
    CameraComparison(idbPath,'Photos0','Photos0');
    CameraComparison(idbPath,'Photos0','Photos-10');
    CameraComparison(idbPath,'Photos0','Photos+10');

    cd ..
end

%% SPECTRAL SENSITIVITY
if (DOSPECTRALSENS)
    cd SPECTRAL_SENSITIVITY
    unix('rm T_camera.mat');
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/SPECTRAL_SENSITIVITY_TRIAL1'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/SPECTRAL_SENSITIVITY_TRIAL1'], 0);
        ProcessNEFToPGM(1, [idbPath '/SPECTRAL_SENSITIVITY_TRIAL2'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/SPECTRAL_SENSITIVITY_TRIAL2'], 0);
        ProcessNEFToPGM(1, [idbPath '/SPECTRAL_SENSITIVITY_2010'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/SPECTRAL_SENSITIVITY_2010'], 0);
    end
    EstimateSpectralSensitivity(idbPath);
    EstimateSpectralSensitivity2010(idbPath);
    cd ..
    MakeCamCal;
end

%% MCC CHECK
if (DOMCC)
    cd MCC
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/MCC'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/MCC'], 0);
        ProcessRawMatToCalibrated(1,[idbPath '/MCC']);
    end
    EstimateMCC(idbPath);
    cd ..
end

%% MTF
if (DOMTF)
    cd MTF
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/MTF/ORIG'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/MTF/ORIG'], 0);
        ProcessRawMatToCalibrated(1,[idbPath '/MTF/ORIG'],0,1);
        ProcessNEFToPGM(1, [idbPath '/MTF/NEW2010'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/MTF/NEW2010'], 0);
        ProcessRawMatToCalibrated(1,[idbPath '/MTF/NEW2010'],0,1);
    end
    EstimateMTF(idbPath,1);
    EstimateMTFNotUsed(idbPath);
    cd ..
end

%% LUMCHECK
if (DOLUMCHECK)
    cd LUMCHECK
    if (REPROCESS == 1)
        ProcessNEFToPGM(1, [idbPath '/LUMCHECK/LUM320'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/LUMCHECK/LUM320'], 0);
        ProcessRawMatToCalibrated(1,[idbPath '/LUMCHECK/LUM320'],0,1);
        ProcessNEFToPGM(1, [idbPath '/LUMCHECK/LUM368'],dcrawPath);
        ProcessPGMToRawMat(1, [idbPath '/LUMCHECK/LUM368'], 0);
        ProcessRawMatToCalibrated(1,[idbPath '/LUMCHECK/LUM368'],0,1);
    end
    EstimateLumCheck(idbPath);
    cd ..
end

return





