% CameraComparison
%
% Compare the responses of the two cameras to the same set of stimuli
%
% 11/16/10  dhb  Started on this, based on the exposure linearity code.

function CameraComparison(idbPath,whichPhilly,whichBots)

% Close
close all;

%% Good response range
fitLowResp = 50;
fitHighResp = 16100;

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Clear out old stuff on general principles;
curDir = pwd;

%% If 1, plot where in the image the signal is being extracted from
CHECK = 1;

%% Choose which series
if (nargin < 2 || isempty(whichPhilly))
    whichPhilly = 'Photos0';
end
if (nargin < 3 || isempty(whichBots))
    whichBots = 'Photos0';
end

%% Analyze the exposure series from the auxiliary (botswana) camera.
cd([idbPath '/CAMERACOMPARISON/BotsCam/' whichBots]);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% These are the coordinates of the white standard in the exposure
% linearity photo set; nQ x nQ is the size of the region in pixels
nQ          = 100;
crop_coords = [450 700];


% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
if (CHECK == 1)
    igFig = figure;
end
for f = 1:length(theFiles)
    realFile = f;
    [~,filenameReal] = fileparts(theFiles(realFile).name);
    fprintf('Processing file %s\n',filenameReal);

    % Get exposure duration, etc
    imageInfoReal = GetNEFInfo(filenameReal);
    
    % Get exposure duration, etc
    fprintf('\t\tCamera: %s\n',imageInfoReal.whichCamera);
    fprintf('\t\tExposure %g\n',imageInfoReal.exposure);
    fprintf('\t\tfStop %g\n',imageInfoReal.fStop);
    fprintf('\t\tISO %g\n',imageInfoReal.ISO);
    
    % Some checks
    if (imageInfoReal.ISO ~= 400)
        error('All data should be at ISO 400');
    end
    if (imageInfoReal.fStop ~= 5.6)
        error('All data should be at fStop 5.6');
    end
    if (~strcmp(imageInfoReal.whichCamera,'auxiliary'))
        error('This should be auxiliary camera\n');
    end
    
    % Load the camera data for the image
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    
    % Read in raw image, dark subtract, and truncate negative values
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);
    
    % Diagnostic plot of where we are extracting from.  There is a surprising
    % amount of image-to-image movement here, but all of the images do seem
    % to be drawing their pixels from the white standard.
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        
        ig(crop_coords(1):2:crop_coords(1)+nQ, crop_coords(2):2:crop_coords(2)+nQ,:) = 1;
        ig(crop_coords(1)+1:2:crop_coords(1)+nQ, crop_coords(2)+1:2:crop_coords(2)+nQ,:) = 0;
        
        figure(igFig);clf;imagesc(ig);
        if (f == 20)
            cd(curDir);
            imwrite(ig,['AuxPhotoRegions_' whichBots '.jpg'],'jpg');
            cd([idbPath '/CAMERACOMPARISON/BotsCam/' whichBots]);
        end
    end
    
    % Select plaquettes from the image
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
    end
    
    % Build the linearity table for different exposure and fstop
    auxLinearityTable(f,1) = imageInfoReal.exposure;
    auxLinearityTable(f,2:4)= croppedImage; 
    
end
cd(curDir);

% Aux camera plot
f1 = figure; clf
position = get(f1,'Position');
position(3) = 680; position(4) = 680;
set(f1,'Position',position);
warning('off','curvefit:fit:noStartPoint');

fitrange{1} = find((auxLinearityTable(:,2) >= fitLowResp) & (auxLinearityTable(:,2) < fitHighResp));
fitrange{2} = find((auxLinearityTable(:,3) >= fitLowResp) & (auxLinearityTable(:,3) < fitHighResp));
fitrange{3} = find((auxLinearityTable(:,4) >= fitLowResp) & (auxLinearityTable(:,4) < fitHighResp));
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope 1 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(auxLinearityTable(fitrange{i},1)), log10(auxLinearityTable(fitrange{i},1+i)),ft);

    figure(f1);
    loglog(auxLinearityTable(:,1), auxLinearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(auxLinearityTable(fitrange{i},1), auxLinearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(auxLinearityTable(fitrange{i},1), 10.^(feval(pf, log10(auxLinearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(auxLinearityTable(fitrange{i},1)), log10(auxLinearityTable(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end
axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
title('Auxiliary Camera');
xlim([10^-4 100]); ylim([1 100000]);

% Dump free slopes and min exposure
for i = 1:3
    fprintf('Free fit slope for auxiliary camera channel %d is %0.2f\n',i,freeslope(i));
end
minExposure = min([auxLinearityTable(fitrange{1},1) ; auxLinearityTable(fitrange{2},1) ; auxLinearityTable(fitrange{3},1)]);
maxExposure = max([auxLinearityTable(fitrange{1},1) ; auxLinearityTable(fitrange{2},1) ; auxLinearityTable(fitrange{3},1)]);
fprintf('Auxiliary camera: min exposure = %f, max = %f\n',minExposure,maxExposure);
FigureSave(['AuxExposureFig_' whichBots '.pdf'],f1,'pdf');

%% Analyze the exposure series from the auxiliary (philly) camera.
cd([idbPath '/CAMERACOMPARISON/PhillyCam/' whichPhilly]);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

nQ          = 100;
crop_coords = [455 690];

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
if (CHECK == 1)
    igFig = figure;
end
if (length(theFiles) ~= size(auxLinearityTable,1))
    error('Mismatch in number of standard and auxiliary camera images');
end
for f = 1:length(theFiles)
    realFile = f;
    [~,filenameReal] = fileparts(theFiles(realFile).name);
    fprintf('Processing file %s\n',filenameReal);

    % Get exposure duration, etc
    imageInfoReal = GetNEFInfo(filenameReal);
    
    % Get exposure duration, etc
    fprintf('\t\tCamera: %s\n',imageInfoReal.whichCamera);
    fprintf('\t\tExposure %g\n',imageInfoReal.exposure);
    fprintf('\t\tfStop %g\n',imageInfoReal.fStop);
    fprintf('\t\tISO %g\n',imageInfoReal.ISO);
    
    % Some checks
    if (imageInfoReal.ISO ~= 400)
        error('All data should be at ISO 400');
    end
    if (imageInfoReal.fStop ~= 5.6)
        error('All data should be at fStop 5.6');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Load the camera data for the image
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    
    % Read in raw image, dark subtract, and truncate negative values
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);
    
    % Diagnostic plot of where we are extracting from.  There is a surprising
    % amount of image-to-image movement here, but all of the images do seem
    % to be drawing their pixels from the white standard.
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        
        ig(crop_coords(1):2:crop_coords(1)+nQ, crop_coords(2):2:crop_coords(2)+nQ,:) = 1;
        ig(crop_coords(1)+1:2:crop_coords(1)+nQ, crop_coords(2)+1:2:crop_coords(2)+nQ,:) = 0;
        
        figure(igFig);clf;imagesc(ig);
        if (f == 20)
            cd(curDir);
            imwrite(ig,['StandardPhotoRegions_' whichPhilly '.jpg'],'jpg');
            cd([idbPath '/CAMERACOMPARISON/PhillyCam/' whichPhilly]);
        end
    end
    
    % Select plaquettes from the image
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
    end
    
    % Build the linearity table for different exposure and fstop
    stanLinearityTable(f,1) = imageInfoReal.exposure;
    if (stanLinearityTable(f,1) ~= auxLinearityTable(f,1))
        error('Mismatch between standard and auxiliary exposure times');
    end
    stanLinearityTable(f,2:4)= croppedImage; 
    
end
cd(curDir);

% Standard camera plot
f2 = figure; clf
position = get(f2,'Position');
position(3) = 680; position(4) = 680;
set(f2,'Position',position);

fitrange{1} = find((stanLinearityTable(:,2) >= fitLowResp) & (stanLinearityTable(:,2) < fitHighResp));
fitrange{2} = find((stanLinearityTable(:,3) >= fitLowResp) & (stanLinearityTable(:,3) < fitHighResp));
fitrange{3} = find((stanLinearityTable(:,4) >= fitLowResp) & (stanLinearityTable(:,4) < fitHighResp));
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope 1 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(stanLinearityTable(fitrange{i},1)), log10(stanLinearityTable(fitrange{i},1+i)),ft);

    figure(f2);
    loglog(stanLinearityTable(:,1), stanLinearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(stanLinearityTable(fitrange{i},1), stanLinearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(stanLinearityTable(fitrange{i},1), 10.^(feval(pf, log10(stanLinearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    % Determine slope if it is allowed to go free
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(stanLinearityTable(fitrange{i},1)), log10(stanLinearityTable(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end
axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
title('Standard Camera');
xlim([10^-4 100]); ylim([1 100000]);

% Dump free slopes and min exposure
for i = 1:3
    fprintf('Free fit slope for standard camera channel %d is %0.2f\n',i,freeslope(i));
end
minExposure = min([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},1) ; stanLinearityTable(fitrange{3},1)]);
maxExposure = max([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},1) ; stanLinearityTable(fitrange{3},1)]);
fprintf('Standard camera: min exposure = %f, max = %f\n',minExposure,maxExposure);
FigureSave(['StandardExposureFig_' whichPhilly '.pdf'],f2,'pdf');

%% Now plot one versus the other and find ratio
f3 = figure; clf
position = get(f3,'Position');
position(3) = 680; position(4) = 680;
set(f3,'Position',position);

fitrange{1} = find((auxLinearityTable(:,2) >= fitLowResp) & (auxLinearityTable(:,2) < fitHighResp) & ...
    (stanLinearityTable(:,2) >= fitLowResp) & (stanLinearityTable(:,2) < fitHighResp));
fitrange{2} = find((auxLinearityTable(:,3) >= fitLowResp) & (auxLinearityTable(:,3) < fitHighResp) & ...
    (stanLinearityTable(:,3) >= fitLowResp) & (stanLinearityTable(:,3) < fitHighResp));
fitrange{3} = find((auxLinearityTable(:,4) >= fitLowResp) & (auxLinearityTable(:,4) < fitHighResp) & ...
    (stanLinearityTable(:,4) >= fitLowResp) & (stanLinearityTable(:,4) < fitHighResp));

ft = '1*x+b';
pfall = fit(log10([stanLinearityTable(fitrange{1},2) ; stanLinearityTable(fitrange{2},3) ; stanLinearityTable(fitrange{3},4)]), ...
	log10([auxLinearityTable(fitrange{1},2) ; auxLinearityTable(fitrange{2},3); auxLinearityTable(fitrange{3},4)]),ft);
overallSlope = 10.^pfall.b;

clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope 1 in log-log through the points of interest to get slope for each channel
    pf = fit(log10(stanLinearityTable(fitrange{i},1+i)),log10(auxLinearityTable(fitrange{i},1+i)),ft);
    freeslope(i) = 10.^pf.b;

    figure(f3);
    loglog(stanLinearityTable(:,1+i), auxLinearityTable(:,1+i),'o','color',clx{i}); hold on
    loglog(stanLinearityTable(fitrange{i},1+i), auxLinearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
end
loglog([stanLinearityTable(fitrange{1},2) ; stanLinearityTable(fitrange{2},3) ; stanLinearityTable(fitrange{3},4)], ...
    10.^feval(pfall,log10([stanLinearityTable(fitrange{1},2) ; stanLinearityTable(fitrange{2},3) ; stanLinearityTable(fitrange{3},4)])), ...
    '-','color','k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('Standard Camera Raw RGB (dark subtracted)','fontsize',14);
ylabel('Auxiliary Camera Raw RGB (dark subtracted)','fontsize',14);
title(sprintf('Slope %0.2f',overallSlope));
xlim([10 100000]); ylim([10 100000]);

% Dump free slopes and min exposure
for i = 1:3
    fprintf('Aux vs standard slope for channel %d is %0.2f\n',i,freeslope(i));
end
fprintf('Overall aux vs standard slope is %0.2f\n',overallSlope);
FigureSave(['CameraComparisonFig1_' whichPhilly '_' whichBots '.pdf'],f3,'pdf');

warning('on','curvefit:fit:noStartPoint');

return

