% EstimateSpectralSensitivity2010
%
% This program estimates the D70 spectral sensitivity
% directly from the photo data.  Assumes that ProcessNEFToPGM
% has been run on both the 2010 directory.
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
extractRows = 475:515;
extractCols = 800:855;
CHECK = 1;
S_power = [380 4 101];
wls_power = SToWls(S_power);

%% First process the ascending wavelength data
cd([idbPath '/SPECTRAL_SENSITIVITY_2010']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);
if (length(theFiles)/2 ~= 31)
    error('Surprising number of .NEF files');
end

% Loop over files and analyze
if (CHECK == 1)
    igFig = figure;
    powFig = figure; hold on
end
for f = 1:length(theFiles)/2
    wavelength = 400+(f-1)*10;
    filenameReal = sprintf('DSC_%d',wavelength);
    filenameDark = sprintf('DSC_%dd',wavelength);
    matfile = sprintf('m%d',wavelength);
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
    
    % Get power
    load(matfile);
    eval(sprintf('curSpectrum = m%d;',wavelength));
    lightPower(f) = sum(curSpectrum);
    
    % Diagnostic plot of where we are extracting from.
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(extractRows,extractCols,:) = 1;
        figure(igFig);clf;imagesc(ig);
        if (f == 15)
            cd(curDir);
            imwrite(ig,'Regions2010.jpg','jpg');
            cd([idbPath '/SPECTRAL_SENSITIVITY_2010']);
        end
        
        figure(powFig);
        plot(wls_power,curSpectrum,'k');
        drawnow;
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
        theSensitivity(k,f) = scaleFactor*(diffRGB1(f,k))/lightPower(f);
    end
end

% Return home
cd(curDir);

%% Save spectrum figure
FigureSave('MonoSpectra2010.pdf',powFig,'pdf');

%% Write out average spectral sensitivity in standard form
T_camera = theSensitivity;
S_camera = [400 10 31];
save T_camera2010 T_camera S_camera

%% Plot spectral sensitivities
ssPlot = figure; clf; hold on
clx={'r';'g';'b'};
for i=1:3,
    plot(SToWls(S_camera),(T_camera(i,:)'),clx{i},'linewidth',1.3);
end
axis square;
set(gca,'fontsize',14);
xlabel('Wavelength (nm)','fontsize',14);
ylabel('Camera sensitivity','fontsize',14);
%title('Camera spectral sensitivity','fontsize',14);
FigureSave(['SpectralSensitivities2010.pdf'],ssPlot,'pdf');

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
FigureSave('ConeSensitivities2010.pdf',conePlot,'pdf');

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
FigureSave('ApproximatedConeSensitivities2010.pdf',rgbConePlot,'pdf');

%% Save matrix converting RGB to LMS
save M_RGBToLMS2010 M_RGBToLMS M_LMSToLum

end

function md = median2(inp)
    md = median(inp(:));
end


