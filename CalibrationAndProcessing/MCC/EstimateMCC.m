% EstimateMCC
%
% We took images from MCC and extracted 24 color squares. For each square,
% we have radiometric readings that can be directly transformed into LMS
% isomerizations. We can compare this to the result of our image processing
% pipeline, that takes the photos of MCC, goes from raw camera RGB to LMS
% to isomerizations. We expect consistent results of camera vs
% spectrometer.
%
% Also processes cone info and provides various diagnostic plots.  Outputs
% isomerization rate conversion factors.
%
% Parameters:
%   idbPath -- the path to the root of the calibration image database
%
% 07/08/10  gt  Wrote it.
% 12/17/10  dhb Next revision.
% 12/21/10  dhb Add checks of ProcessRawMatToCalbrated output.
% 06/20/17  dhb Isomerization check misses by 0.1%.  Not sure why.  Not
%               worrying about it, yet.

function EstimateMCC(idbPath)

% Close
close all;

%% Good response range
fitLowResp = 50;
fitHighResp = 16100;

%% Block size
mccBlkSize = 12;

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Clear out old stuff on general principles;
curDir = pwd;

%% If 1, plot where in the image the signal is being extracted from
CHECK = 1;

%% These are the coordinates of the MCC check box being viewed
%% nQ x nQ is the size of the region in pixels
nQ          = 50;
crop_coords = [470 755];

%% Load the MCC spectra
spd_MCC = load('MCCSpectra');

%% Load the photos.
cd([idbPath '/MCC']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
if (CHECK)
    igFig = figure; clf;
end
mccFromCameraRawRGB = zeros(3,length(theFiles));
mccFromCameraDarkRGB = zeros(3,length(theFiles));
mccFromCameraRGB = zeros(3,length(theFiles));
for f = 1:length(theFiles)
    realFile = f;
    [nil,filenameReal] = fileparts(theFiles(realFile).name);
    fprintf('Processing file %s\n',filenameReal);
    
    % Get exposure duration, etc
    imageInfoReal = GetNEFInfo(filenameReal);
    
    % Get exposure duration, etc
    fprintf('\t\tCamera: %s\n',imageInfoReal.whichCamera);
    fprintf('\t\tExposure %g\n',imageInfoReal.exposure);
    fprintf('\t\tfStop %g\n',imageInfoReal.fStop);
    fprintf('\t\tISO %g\n',imageInfoReal.ISO);
    
    % Some checks
    if (imageInfoReal.ISO ~= 1000)
        error('All data should be at ISO 1000');
    end
    if (imageInfoReal.fStop ~= 1.8)
        error('All data should be at f 1.8');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Read in raw image, dark subtract, and truncate negative values
    load([filenameReal '.raw.mat']);
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);
    scaleFactor = GetStandardizingCameraScaleFactor(imageInfoReal);

    % Extract image data and scale to standardized representation
    for k=1:3,
        mccFromCameraRawRGB(k,f) = mean(mean(double(theImage.rawCameraRGB(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,k))));
        mccFromCameraDarkRGB(k,f) = mean(mean(double(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,k))));
        mccFromCameraRGB(k,f)  = scaleFactor*mccFromCameraDarkRGB(k,f);
    end
    
    if (CHECK)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,:) = 1;
        figure(igFig);clf;imagesc(ig);
        if (f == 10)
            cd(curDir);
            imwrite(ig,'PhotoRegions.jpg','jpg');
            cd([idbPath '/MCC']);
        end
    end
    
    % If the calibrated output files have been produced, read them and get values.  This provides a check
    % that the calibrated values are what they should be.
    if (exist([filenameReal '_LMS.mat'],'file'))
        load([filenameReal '_LMS']);
        load([filenameReal '_RGB']);
        load([filenameReal '_LUM']);
        for k=1:3,
            isomFromFile(k,f) = mean(mean(double(LMS_Image(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,k))));
            rgbFromFile(k,f) = mean(mean(double(RGB_Image(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,k))));
        end
        lumFromFile(f) = mean(mean(double(LUM_Image(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ))));
        clear LMS_Image RGB_Image LUM_Image
    end
end
cd(curDir)

%% Check for saturation
if (any(mccFromCameraRawRGB >= 16300))
    fprintf('WARNING: Saturated RGB values somewhere');
end

%% Get predicted RGB from measured spectra
T_MCC = SplineCmf(CamCal.S_camera,CamCal.T_camera,spd_MCC.S);
mccFromRadiometerRGB = T_MCC*spd_MCC.spd;

rgbPlot = figure; clf;
loglog(mccFromRadiometerRGB(1,:), mccFromCameraRGB(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccFromRadiometerRGB(1,index), mccFromCameraRGB(1,index),'ro','MarkerFaceColor','r');
loglog(mccFromRadiometerRGB(2,:), mccFromCameraRGB(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccFromRadiometerRGB(2,index), mccFromCameraRGB(2,index),'go','MarkerFaceColor','g');
loglog(mccFromRadiometerRGB(3,:), mccFromCameraRGB(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccFromRadiometerRGB(3,index), mccFromCameraRGB(3,index),'bo','MarkerFaceColor','b');
loglog([1e4 1e7],[1e4 1e7],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('Standardized RGB from radiometer','fontsize',14);
ylabel('Standardized RGB from camera','fontsize',14);
xlim([1e4 1e7]); ylim([1e4 1e7]);
FigureSave('MCCRGBCheck.pdf',rgbPlot,'pdf');

%% Cone checks
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,spd_MCC.S);
mccFromRadiometerLMS = T_cones*spd_MCC.spd;
mccFromRadiometerRGBLMS = CamCal.M_RGBToLMS*mccFromRadiometerRGB;
mccFromCameraRGBLMS = CamCal.M_RGBToLMS*mccFromCameraRGB;

coneFromCameraRGBPlot = figure; clf;
loglog(mccFromRadiometerLMS(1,:), mccFromCameraRGBLMS(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccFromRadiometerLMS(1,index), mccFromCameraRGBLMS(1,index),'ro','MarkerFaceColor','r');
loglog(mccFromRadiometerLMS(2,:), mccFromCameraRGBLMS(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccFromRadiometerLMS(2,index), mccFromCameraRGBLMS(2,index),'go','MarkerFaceColor','g');
loglog(mccFromRadiometerLMS(3,:), mccFromCameraRGBLMS(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccFromRadiometerLMS(3,index), mccFromCameraRGBLMS(3,index),'bo','MarkerFaceColor','b');
loglog([1e-4 1e-1],[1e-4 1e-1],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('LMS from radiometer','fontsize',14);
ylabel('LMS from camera RGB','fontsize',14);
xlim([1e-4 1e-1]); ylim([1e-4 1e-1]);
FigureSave('MCCLMSFromCameraRGBCheck.pdf',coneFromCameraRGBPlot,'pdf');

coneFromRadiometerRGBPlot  = figure; clf;
loglog(mccFromRadiometerLMS(1,:), mccFromRadiometerRGBLMS(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccFromRadiometerLMS(1,index), mccFromRadiometerRGBLMS(1,index),'ro','MarkerFaceColor','r');
loglog(mccFromRadiometerLMS(2,:), mccFromRadiometerRGBLMS(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccFromRadiometerLMS(2,index), mccFromRadiometerRGBLMS(2,index),'go','MarkerFaceColor','g');
loglog(mccFromRadiometerLMS(3,:), mccFromRadiometerRGBLMS(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccFromRadiometerLMS(3,index), mccFromRadiometerRGBLMS(3,index),'bo','MarkerFaceColor','b');
loglog([1e-4 1e-1],[1e-4 1e-1],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('LMS from radiometer','fontsize',14);
ylabel('LMS from radiometer RGB','fontsize',14);
xlim([1e-4 1e-1]); ylim([1e-4 1e-1]);
FigureSave('MCCLMSFromRadiometerRGBCheck.pdf',coneFromRadiometerRGBPlot,'pdf');

%% Isomerization rates
for f = 1:size(spd_MCC.spd,2);
    [lRate(f),mRate(f),sRate(f)] = IsomerizationsInEyeFunction(spd_MCC.spd(:,f),spd_MCC.S);
end
mccIsomersizationsFromRadiometer =[lRate ; mRate ; sRate];
for k = 1:3
    LMSToIsomerizations(k) = (mccFromRadiometerLMS(k,:)')\(mccIsomersizationsFromRadiometer(k,:)');
end
save LMSToIsomerizations LMSToIsomerizations

isomFromCameraRGBPlot = figure; clf;
loglog(mccIsomersizationsFromRadiometer(1,:),LMSToIsomerizations(1)*mccFromCameraRGBLMS(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(1,index),LMSToIsomerizations(1)*mccFromCameraRGBLMS(1,index),'ro','MarkerFaceColor','r');
loglog(mccFromRadiometerLMS(2,:),LMSToIsomerizations(2)*mccFromCameraRGBLMS(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(2,index),LMSToIsomerizations(2)*mccFromCameraRGBLMS(2,index),'go','MarkerFaceColor','g');
loglog(mccIsomersizationsFromRadiometer(3,:),LMSToIsomerizations(3)*mccFromCameraRGBLMS(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(3,index),LMSToIsomerizations(3)*mccFromCameraRGBLMS(3,index),'bo','MarkerFaceColor','b');
loglog([1e1 1e4],[1e1 1e4],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('LMS isomerizations from radiometer','fontsize',14);
ylabel('LMS isomerizations from camera RGB','fontsize',14);
xlim([1e1 1e4]); ylim([1e1 1e4]);
FigureSave('MCCIsomLMSFromCameraRGBCheck.pdf',isomFromCameraRGBPlot,'pdf');

isomFromRadiometerRGBPlot  = figure; clf;
loglog(mccIsomersizationsFromRadiometer(1,:),LMSToIsomerizations(1)*mccFromRadiometerRGBLMS(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(1,index),LMSToIsomerizations(1)*mccFromRadiometerRGBLMS(1,index),'ro','MarkerFaceColor','r');
loglog(mccIsomersizationsFromRadiometer(2,:),LMSToIsomerizations(2)*mccFromRadiometerRGBLMS(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(2,index),LMSToIsomerizations(2)*mccFromRadiometerRGBLMS(2,index),'go','MarkerFaceColor','g');
loglog(mccIsomersizationsFromRadiometer(3,:),LMSToIsomerizations(3)*mccFromRadiometerRGBLMS(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(3,index),LMSToIsomerizations(3)*mccFromRadiometerRGBLMS(3,index),'bo','MarkerFaceColor','b');
loglog([1e1 1e4],[1e1 1e4],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('LMS isomerizations from radiometer','fontsize',14);
ylabel('LMS isomerizations from radiometer RGB','fontsize',14);
xlim([1e1 1e4]); ylim([1e1 1e4]);
FigureSave('MCCIsomLMSFromRadiometerRGBCheck.pdf',isomFromRadiometerRGBPlot,'pdf');

isomFromRadiometerPlot  = figure; clf;
loglog(mccIsomersizationsFromRadiometer(1,:),LMSToIsomerizations(1)*mccFromRadiometerLMS(1,:),'ro');  hold on;
index = find(mccFromCameraDarkRGB(1,:) < fitHighResp & mccFromCameraDarkRGB(1,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(1,index),LMSToIsomerizations(1)*mccFromRadiometerRGBLMS(1,index),'ro','MarkerFaceColor','r');
loglog(mccIsomersizationsFromRadiometer(2,:),LMSToIsomerizations(2)*mccFromRadiometerLMS(2,:),'go');
index = find(mccFromCameraDarkRGB(2,:) < fitHighResp & mccFromCameraDarkRGB(2,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(2,index),LMSToIsomerizations(2)*mccFromRadiometerLMS(2,index),'go','MarkerFaceColor','g');
loglog(mccIsomersizationsFromRadiometer(3,:),LMSToIsomerizations(3)*mccFromRadiometerLMS(3,:),'bo');
index = find(mccFromCameraDarkRGB(3,:) < fitHighResp & mccFromCameraDarkRGB(3,:) > fitLowResp);
loglog(mccIsomersizationsFromRadiometer(3,index),LMSToIsomerizations(3)*mccFromRadiometerLMS(3,index),'bo','MarkerFaceColor','b');
loglog([1e1 1e4],[1e1 1e4],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('LMS isomerizations from radiometer','fontsize',14);
ylabel('LMS isomerizations from radiometer RGB','fontsize',14);
xlim([1e1 1e4]); ylim([1e1 1e4]);
FigureSave('MCCIsomLMSFromRadiometerCheck.pdf',isomFromRadiometerPlot,'pdf');

%% Check luminance
luminanceFromCamera = CamCal.M_LMSToLum*mccFromCameraRGBLMS;
load T_ss2000_Y2
T_lum = SplineCmf(S_ss2000_Y2,683*T_ss2000_Y2,spd_MCC.S);
luminanceFromRadiometer = T_lum*spd_MCC.spd;

lumPlot = figure; clf;
loglog(luminanceFromRadiometer ,luminanceFromCamera ,'ko','MarkerFaceColor','k'); hold on
loglog([1e0 1e2],[1e0 1e2],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('Luminance from radiometer','fontsize',14);
ylabel('Luminance from camera','fontsize',14);
xlim([1e0 1e2]); ylim([1e0 1e2]);
FigureSave('MCCLuminanceCheck.pdf',lumPlot,'pdf');
save MCCLuminanceCheckData luminanceFromRadiometer luminanceFromCamera

%% Checking that calibrated files contain what we expect
% See ProcessRawMatToCalibrated.m
if (exist('isomFromFile','var'))
    fprintf('Checking calibrated RGB values from file ...');
    if (max(max(abs(mccFromCameraRawRGB(:)-rgbFromFile(:)))) ~= 0)
        fprintf(' BAD!\n');
        fprintf('\tMax absoulte error is %g\n',max(max(abs(mccFromCameraRGB(:)-rgbFromFile(:)))));
        fprintf('\tMean value of data is %g\n',mean(mean(mccFromCameraRGB(:))));
    else
        fprintf(' OK!\n');
    end
    
    fprintf('Checking isomerization values from file ...');
    for k = 1:3
        isomFromHere(k,:) = LMSToIsomerizations(k)*mccFromCameraRGBLMS(k,:);
    end
end
if (max(max(abs((isomFromHere(:)-isomFromFile(:)))./isomFromHere(:))) >= 1e-5)
    fprintf(' BAD!\n');
    fprintf('\tMax absoulte error is %g\n',max(max(abs(isomFromHere(:)-isomFromFile(:)))));
    fprintf('\tMean value of data being compared is %g\n',mean(mean(abs(isomFromHere(:)))));
else
    fprintf(' OK!\n');
end

fprintf('Checking luminance values from file ...');
if (max(max(abs((luminanceFromCamera-lumFromFile)./luminanceFromCamera))) >= 1e-5)
    fprintf(' BAD!\n');
    fprintf('\tMax absoulte error is %g\n',max(max(abs(luminanceFromCamera-lumFromFile))));
else
    fprintf(' OK!\n');
end

end


