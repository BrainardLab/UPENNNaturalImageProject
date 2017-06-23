% EstimateApertureLinearity
%
% We took a bunch of images in the dark at various apertures, to check
% whether the image intensity really varies as fStop^2. This will load the
% image set and plot the three color channels against fStop.
%
% Parameters:
%   idbPath -- path to the root of the calibration image database, if
%   empty, use current
%
% 07/08/10  gt   Wrote it.
% 11/10/10  dhb  Add a default path.
% 11/10/10  dhb  Set any pixels less than 0 to 0 after dark subtraction.
%           dhb  Adjust extracted region to be better centered on the ref square.
%           dhb  Print out slope of free fit.
%           dhb  Add exposure check.
% 11/12/10  dhb  Log base 10 for fits (silly, I know)
%           dhb  Find range automatically from response limits
%           dhb  Dark subtract through function.
%           dhb  No legend
% 1/5/10    dhb  Dump text to file, and also mean absolute log fit error.


%% Clear out old stuff on general principles;
function EstimateApertureLinearity(idbPath)

%% Close old figs
close all;

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end
curDir = pwd;

%% If 1, plot where in the image the signal is being extracted from
CHECK = 0;

%% These are the coordinates of the white standard in the exposure
%% linearity photo set; nQ x nQ is the size of the region in pixels
nQ          = 300;
crop_coords = [350 630];

secondary_nQ = 150;
secondary_crop_coords = [450 200];

%% Load images
cd([idbPath '/APERTURE']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
if (CHECK == 1)
    igFig = figure;
end
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
    if (imageInfoReal.exposure ~= 0.004)
        error('All data should be for exposure 0.004');
    end
    if (imageInfoReal.ISO ~= 1000)
        error('All data should be at ISO 1000');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end  
       
    % Load the camera data for the image
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    
    % Read in and dark correct raw image
    load([filenameReal '.raw.mat']);
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);

    % Diagnostic plot of where we are extracting from
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        
        ig(crop_coords(1):2:crop_coords(1)+nQ, crop_coords(2):2:crop_coords(2)+nQ,:) = 1;
        ig(crop_coords(1)+1:2:crop_coords(1)+nQ, crop_coords(2)+1:2:crop_coords(2)+nQ,:) = 0;
        
        ig(secondary_crop_coords(1):2:secondary_crop_coords(1)+secondary_nQ, secondary_crop_coords(2):2:secondary_crop_coords(2)+secondary_nQ,:) = 1;
        ig(secondary_crop_coords(1)+1:2:secondary_crop_coords(1)+secondary_nQ, secondary_crop_coords(2)+1:2:secondary_crop_coords(2)+secondary_nQ,:) = 0;
        
        figure(igFig);clf;imagesc(ig);
        if (f == 5)
            cd(curDir);
            imwrite(ig,'ApertureRegions.jpg','jpg');
            cd([idbPath '/APERTURE']);
        end
    end
    
    % Extract image patches
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
        secondaryCroppedImage(i) = mean(mean(realImage(secondary_crop_coords(1):secondary_crop_coords(1)+secondary_nQ, ...
            secondary_crop_coords(2):secondary_crop_coords(2)+secondary_nQ, i)));
    end
    
    % Build linearity table
    linearityTable(f,1) = imageInfoReal.fStop;
    linearityTable(f,2:4)= croppedImage;
    secondaryLinearityTable(f,1) = imageInfoReal.fStop;
    secondaryLinearityTable(f,2:4) = secondaryCroppedImage; 
    
    % Standardize data, to check how well this routine stabilizes the data
    camFactor = GetStandardizingCameraScaleFactor(imageInfoReal);
    standardTable(f,1:3) = camFactor*croppedImage;
    secondaryStandardTable(f,1:3) = camFactor*secondaryCroppedImage;
end
cd(curDir);


%% Plot main result
% Select the data points to fit on -- exclude the saturated region.
fitLowResp = 50;
fitHighResp = 16100;
warning('off','curvefit:fit:noStartPoint');

fitrange{1} = find((linearityTable(:,2) >= fitLowResp) & (linearityTable(:,2) < fitHighResp));
fitrange{2} = find((linearityTable(:,3) >= fitLowResp) & (linearityTable(:,3) < fitHighResp));
fitrange{3} = find((linearityTable(:,4) >= fitLowResp) & (linearityTable(:,4) < fitHighResp));
% fitrange{1} = 8:size(linearityTable,1);
% fitrange{2} = 6:size(linearityTable,1);
% fitrange{3} = 3:size(linearityTable,1);

f1 = figure; clf;
f3 = figure; clf;
f2 = figure; clf
position = get(f2,'Position');
position(3) = 760; position(4) = 480;
set(f2,'Position',position);
subplot(1,2,1);
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('-2*x+a');
    pf = fit(log10(linearityTable(fitrange{i},1)), log10(linearityTable(fitrange{i},1+i)),ft);
    
    % Plot the basic data
    figure(f1);
    ax1(i)=loglog(linearityTable(:,1), linearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(linearityTable(fitrange{i},1), linearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(linearityTable(fitrange{i},1), 10.^(feval(pf, log10(linearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);

    figure(f2);
    ax(i)=loglog(linearityTable(:,1), linearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(linearityTable(fitrange{i},1), linearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(linearityTable(fitrange{i},1), 10.^(feval(pf, log10(linearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    figure(f3);
    loglog(linearityTable(:,1), standardTable(:,i),'o','color',clx{i});hold on;
    loglog(linearityTable(fitrange{i},1), standardTable(fitrange{i},i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    
    % Compute prediction error as mean of absolute log unit error
    meanLogError(i) = mean(abs(log10(linearityTable(fitrange{i},1+i)) - (feval(pf, log10(linearityTable(fitrange{i},1))))));
    maxLogError(i) = max(abs(log10(linearityTable(fitrange{i},1+i)) - (feval(pf, log10(linearityTable(fitrange{i},1))))));

    % Determine slope if it is allowed to go fre
    ft1 = fittype('-b*x+a');
    pf1 = fit(log10(linearityTable(fitrange{i},1)), log10(linearityTable(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end

figure(f1);
axis square;
set(gca,'fontsize',14);
xlabel('Aperture f-value','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
%title('Aperture linearity (ISO=1000, 1/250s exposure)','fontsize',14);
xlim([1 100]); ylim([10 100000]);

figure(f2);
axis square;
set(gca,'fontsize',14);
xlabel('Aperture f-value','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
xlim([1 100]); ylim([10 100000]);

figure(f3);
axis square;
set(gca,'fontsize',14);
xlabel('Aperture f-value','fontsize',14);
ylabel('Scaled camera RGB (dark subtracted)','fontsize',14);
xlim([1 100]); ylim([1e6 1e8]);

%% Report standard variation in linear units
for i = 1:3
    minStandard = min(standardTable(fitrange{i},i));
    maxStandard = max(standardTable(fitrange{i},i));
    fprintf('Channel %d, min standardized %g, max standardized = %g\n',i,minStandard,maxStandard);
end

%% Plot secondary result
% Select the data points to fit on -- exclude the saturated region. 
fitrange{1} = find((secondaryLinearityTable(:,2) >= fitLowResp) & (secondaryLinearityTable(:,2) < fitHighResp));
fitrange{2} = find((secondaryLinearityTable(:,3) >= fitLowResp) & (secondaryLinearityTable(:,3) < fitHighResp));
fitrange{3} = find((secondaryLinearityTable(:,4) >= fitLowResp) & (secondaryLinearityTable(:,4) < fitHighResp));
% fitrange{1} = 4:size(secondaryLinearityTable,1);
% fitrange{2} = 1:size(secondaryLinearityTable,1)-2;
% fitrange{3} = 1:size(secondaryLinearityTable,1);

figure(f2);
subplot(1,2,2);
clx={'r';'g';'b'};
for i=1:3
    % Plot the basic data
    ax(i)=loglog(secondaryLinearityTable(:,1), secondaryLinearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(linearityTable(fitrange{i},1), secondaryLinearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('-2*x+a');
    pf = fit(log10(secondaryLinearityTable(fitrange{i},1)), log10(secondaryLinearityTable(fitrange{i},1+i)),ft);
    loglog(secondaryLinearityTable(fitrange{i},1), 10.^(feval(pf, log10(secondaryLinearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    % Compute prediction error as mean of absolute log unit error
    meanSecondaryLogError(i) = mean(abs(log10(secondaryLinearityTable(fitrange{i},i+1)) - (feval(pf, log10(secondaryLinearityTable(fitrange{i},1))))));
    maxSecondaryLogError(i) = max(abs(log10(secondaryLinearityTable(fitrange{i},i+1)) - (feval(pf, log10(secondaryLinearityTable(fitrange{i},1))))));
    
    % Determine slope if it is allowed to go fre
    ft1 = fittype('-b*x+a');
    pf1 = fit(log10(secondaryLinearityTable(fitrange{i},1)), log10(secondaryLinearityTable(fitrange{i},1+i)),ft1);
    secondaryFreeslope(i) = pf1.b;
end
axis square;
set(gca,'fontsize',14);
xlabel('Aperture f-value','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
xlim([1 100]); ylim([10 100000]);

%% Save those figures out
FigureSave('ApertureFig.pdf',f1,'pdf');
FigureSave('ApertureFigTwo.pdf',f2,'pdf');
FigureSave('ApertureStandardizedFig.pdf',f3','pdf');
warning('on','curvefit:fit:noStartPoint');



%% Dump info into a text file free slopes
fid = fopen('ApertureLinearityInfo.txt','w');
for i = 1:3
    fprintf(fid,'Free fit slope for channel %d is %0.3f\n',i,freeslope(i));
    fprintf(fid,'Secondary free fit slope for channel %d is %0.2f\n',i,secondaryFreeslope(i));
    fprintf(fid,'Mean absolute log10 error for channel %d is %0.2f\n',i,meanLogError(i));
    fprintf(fid,'Secondary mean absolute log10 error for channel %d is %0.3f\n',i,meanSecondaryLogError(i));
    fprintf(fid,'Max absolute log10 error for channel %d is %0.2f\n',i,maxLogError(i));
    fprintf(fid,'Secondary max absolute log10 error for channel %d is %0.3f\n',i,maxSecondaryLogError(i));
end
fclose(fid);




