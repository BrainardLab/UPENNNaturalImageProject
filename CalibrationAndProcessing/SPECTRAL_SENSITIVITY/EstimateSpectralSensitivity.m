% EstimateSpectralSensitivity
%
% This program estimates the D70 spectral sensitivity
% directly from the photo data.  Assumes that ProcessNEFToPGM
% has been run on both the Trial1 and Trial2 Photos directories
% to produce the raw pgm files.
%
% Parameters:
%   idbPath -- path to the root of the calibration image database
%
% See also
%  ProcessNEFToPGM
%
% 01/25/10  dhb  Wrote this.
% 07/08/10  gt   Updated plots
% 07/18/10  gt, dhb Median rather than mean for dark subtract.
% 12/15/10  dhb  Diagnostics to detect saturated measurements
%           dhb  Get power order convention right for pass 2
%           dhb  Save figures

function EstimateSpectralSensitivity(idbPath)

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Good response range
fitLowResp = 50;
fitHighResp = 16100;

%% Clear out old stuff on general principles;
curDir = pwd;

%% Parameters
extractRows = 457:530;
extractCols = 744:812;
CHECK = 1;

%% Radiometric measurements.  Taken from spreadsheet.  
pass1 = [2.96E-03
6.42E-03
6.60E-03
6.35E-03
8.35E-03
1.06E-02
1.48E-02
1.56E-02
1.92E-02
1.87E-02
2.00E-02
2.63E-02
2.84E-02
3.40E-02
3.71E-02
3.89E-02
3.47E-02
4.28E-02
4.54E-02
4.59E-02
4.92E-02
4.88E-02
5.13E-02
4.51E-02
4.80E-02
3.31E-02
3.87E-02
3.46E-02
2.73E-02
2.23E-02
9.87E-03
];

pass2 = [2.59E-03
5.99E-03
6.59E-03
6.75E-03
8.61E-03
9.83E-03
1.40E-02
1.45E-02
1.85E-02
1.76E-02
1.89E-02
1.97E-02
2.34E-02
2.83E-02
3.56E-02
3.66E-02
3.30E-02
3.48E-02
4.31E-02
4.38E-02
4.56E-02
4.65E-02
3.71E-02
4.25E-02
4.54E-02
3.19E-02
3.63E-02
3.22E-02
2.59E-02
2.03E-02
8.80E-03
];

%% First process the ascending wavelength data
cd([idbPath '/SPECTRAL_SENSITIVITY_TRIAL1']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);


% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
if (CHECK == 1)
    igFig = figure;
end
for f = 1:length(theFiles)/2
    realFile = (f-1)*2+1;
    darkFile = 2*f;
    [nil,filenameReal] = fileparts(theFiles(realFile).name);
    [nil,filenameDark] = fileparts(theFiles(darkFile).name);
    fprintf('Processing file pair index %d, %s and %s\n',f,filenameReal,filenameDark);

    % Get exposure duration, etc
    imageInfoReal = GetNEFInfo(filenameReal);
    imageInfoDark = GetNEFInfo(filenameDark);
    
    % Some checks
    if (imageInfoReal.fStop ~= 1.8 | imageInfoDark.fStop ~= 1.8)
        error('All data should be at fStop 1.8\n');
    end
    if (imageInfoReal.ISO ~= 1000 || imageInfoDark.ISO ~= 1000)
        error('All data should be at ISO 1000');
    end
    if (imageInfoReal.exposure ~= imageInfoDark.exposure)
        error('Real and dark images should have the same exposure duration\n');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard') || ~strcmp(imageInfoDark.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Read in raw image and dark subtract
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    load([filenameDark '.raw.mat']); darkImage = theImage.rawCameraRGB;
    extractRealImage = realImage(extractRows,extractCols,:);
    extractDarkImage = darkImage(extractRows,extractCols,:);
    
    % Diagnostic plot of where we are extracting from.
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(extractRows,extractCols,:) = 1;
        figure(igFig);clf;imagesc(ig);
        if (f == 15)
            cd(curDir);
            imwrite(ig,'FirstPassRegions.jpg','jpg');
            cd([idbPath '/SPECTRAL_SENSITIVITY_TRIAL1']);
        end
    end
    
    % Check for saturation
    for k = 1:3
        maxExtract = max(max((extractRealImage(:,:,k))));
        meanExtract = mean(mean(extractRealImage(:,:,k)));
        medianExtract = median(median(extractRealImage(:,:,k)));
        if (maxExtract > 16300)
            fprintf('\tWARNING: Saturated values channel %d: mean is %g, median is %g, max is %g\n',k,meanExtract,medianExtract,maxExtract);
        end
    end
    
    % Get scale factor
    scaleFactor = GetStandardizingCameraScaleFactor(imageInfoReal);
    
    % Get spectral sensitivity be subtracting dark from real and normalizing by
    % power.
    for k = 1:3
        theRGBReal1(f,k) = median2(extractRealImage(:,:,k));
        theRGBDark1(f,k) = median2(extractDarkImage(:,:,k));
        diffRGB1(f,k) = theRGBReal1(f,k)-theRGBDark1(f,k);
        if (diffRGB1(f,k) < fitLowResp || diffRGB1(f,k) > fitHighResp)
            fprintf('\tWARNING: Median dark subtracted response out of good range, channel %d: %f\n',k,diffRGB1(f,k));
        end
        theSensitivity1(k,f) = scaleFactor*(diffRGB1(f,k))/pass1(f);
    end
end

% Return home
cd(curDir);

%% First process the descending wavelength data
cd([idbPath '/SPECTRAL_SENSITIVITY_TRIAL2']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
for f = 1:length(theFiles)/2
    f1 = length(theFiles)/2 - f + 1;
    realFile = (f1-1)*2+1;
    darkFile = 2*f1;
    [nil,filenameReal] = fileparts(theFiles(realFile).name);
    [nil,filenameDark] = fileparts(theFiles(darkFile).name);
    fprintf('Processing file pair index %d, %s and %s\n',f,filenameReal,filenameDark);

    % Get exposure duration, etc
    imageInfoReal = GetNEFInfo(filenameReal);
    imageInfoDark = GetNEFInfo(filenameDark);
    
    % Some checks
    if (imageInfoReal.fStop ~= 1.8 | imageInfoDark.fStop ~= 1.8)
        error('All data should be at fStop 1.8\n');
    end
    if (imageInfoReal.ISO ~= 1000 || imageInfoDark.ISO ~= 1000)
        error('All data should be at ISO 1000');
    end
    if (imageInfoReal.exposure ~= imageInfoDark.exposure)
        error('Real and dark images should have the same exposure duration\n');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard') || ~strcmp(imageInfoDark.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Read in raw image and dark subtract
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    load([filenameDark '.raw.mat']); darkImage = theImage.rawCameraRGB;
    extractRealImage = realImage(extractRows,extractCols,:);
    extractDarkImage = darkImage(extractRows,extractCols,:);
    
    % Diagnostic plot of where we are extracting from.
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(extractRows,extractCols,:) = 1;
        figure(igFig);clf;imagesc(ig);
        if (f == 15)
            cd(curDir);
            imwrite(ig,'SecondPassRegions.jpg','jpg');
            cd([idbPath '/SPECTRAL_SENSITIVITY_TRIAL2']);
        end
    end
     
    % Check for saturation
    for k = 1:3
        maxExtract = max(max((extractRealImage(:,:,k))));
        meanExtract = mean(mean(extractRealImage(:,:,k)));
        medianExtract = median(median(extractRealImage(:,:,k)));
        if (maxExtract > 16300)
            fprintf('\tWARNING: Saturated values channel %d: mean is %g, median is %g, max is %g\n',k,meanExtract,medianExtract,maxExtract);
        end
    end
    
    % Get scale factor
    scaleFactor = GetStandardizingCameraScaleFactor(imageInfoReal);
    
    % Get spectral sensitivity be subtracting dark from real and normalizing by
    % power.
    for k = 1:3
        theRGBReal2(f,k) = median2(extractRealImage(:,:,k));
        theRGBDark2(f,k) = median2(extractDarkImage(:,:,k));
        diffRGB2(f,k) = theRGBReal2(f,k)-theRGBDark2(f,k);
        if (diffRGB2(f,k) < fitLowResp || diffRGB2(f,k) > fitHighResp)
            fprintf('\tWARNING: Median dark subtracted response out of good range, channel %d: %f\n',k,diffRGB2(f,k));
        end
        theSensitivity2(k,f) = scaleFactor*(diffRGB2(f,k))/pass2(f);
    end
end

% Return home
cd(curDir);

%% Write out average spectral sensitivity in standard form
T_camera = (theSensitivity1 + theSensitivity2)/2;
S_camera = [400 10 31];
save T_camera T_camera S_camera

%% Plot spectral sensitivities
ssPlot = figure; clf; hold on
clx={'r';'g';'b'};
for i=1:3,
    plot(SToWls(S_camera),(theSensitivity1(i,:)'),clx{i},'linewidth',1.3);
    hold on;
    plot(SToWls(S_camera),(theSensitivity2(i,:)'),[clx{i} '--'],'linewidth',1.3);
    plot(SToWls(S_camera),(T_camera(i,:)'),'k','linewidth',1.3);
end
axis square;
set(gca,'fontsize',14);
xlabel('Wavelength (nm)','fontsize',14);
ylabel('Camera sensitivity','fontsize',14);
%title('Camera spectral sensitivity','fontsize',14);
ylim([0 7e7]);
FigureSave(['SpectralSensitivities.pdf'],ssPlot,'pdf');

%% Convert RGB sensitivities to approximate LMS fundamentals  
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S_camera);
load T_ss2000_Y2 
T_Y = SplineCmf(S_ss2000_Y2,683*T_ss2000_Y2,S_camera);
M_RGBToLMS = ((T_camera')\(T_cones'))';
T_cones_check = M_RGBToLMS*(T_camera);
M_LMSToLum = ((T_cones')\(T_Y'))';

%% Plot cone fundamentals
conePlot = figure; clf; hold on
for i=1:3,
    plot(SToWls(S_camera),(T_cones(i,:)'),[clx{i} '-'],'linewidth',1.3);hold on;
end
axis square;
set(gca,'fontsize',14);
xlabel('Wavelength (nm)','fontsize',14);
ylabel('Cone fundamentals','fontsize',14);
ylim([-0.1 1.1])
%title('Norm. cone spectral sensitivity','fontsize',14);
FigureSave('ConeSensitivities.pdf',conePlot,'pdf');

%% Plot approximated cone fundamentals
rgbConePlot = figure; clf; hold on
for i=1:3,
    plot(SToWls(S_camera),(T_cones_check(i,:)'),[clx{i} '-'],'linewidth',1.3); hold on
end
axis square;
set(gca,'fontsize',14);
xlabel('Wavelength (nm)','fontsize',14);
ylabel('Reconstructed cone fundaments','fontsize',14);
ylim([-0.1 1.1])
FigureSave('ApproximatedConeSensitivities.pdf',rgbConePlot,'pdf');

%% Save matrix converting RGB to LMS
save M_RGBToLMS M_RGBToLMS M_LMSToLum

end

function md = median2(inp)
    md = median(inp(:));
end


