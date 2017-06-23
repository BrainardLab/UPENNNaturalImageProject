% CameraComparisonMCC
%
% Compare the responses of the two cameras to mcc image data.
%
% 12/08/10  dhb  Started on this, based on the CameraComparison code.

function CameraComparisonMCC(idbPath,whichPhilly,whichBots)

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

%% Choose which series
if (nargin < 2 || isempty(whichPhilly))
    whichPhilly = '0011';
end
if (nargin < 3 || isempty(whichBots))
    whichBots = '0020';
end

%% Get image data for the auxiliary (botswana) camera.
cd([idbPath '/CAMERACOMPARISON/BotsCam/Color']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% Get the image data
filenameReal = ['DSC_' whichBots];
if (CHECK == 1)
end
fprintf('Processing file %s\n',filenameReal);

% Get exposure duration, etc
imageInfoBots = GetNEFInfo(filenameReal);

% Get exposure duration, etc
fprintf('\t\tCamera: %s\n',imageInfoBots.whichCamera);
fprintf('\t\tExposure %g\n',imageInfoBots.exposure);
fprintf('\t\tfStop %g\n',imageInfoBots.fStop);
fprintf('\t\tISO %g\n',imageInfoBots.ISO);
if (~strcmp(imageInfoBots.whichCamera,'auxiliary'))
    error('This should be auxiliary camera\n');
end

% Read in raw image, dark subtract, and truncate negative values
load([filenameReal '.raw.mat']);
CamCal = LoadCamCal(imageInfoBots.whichCamera);
realImage = DarkCorrect(CamCal,imageInfoBots,theImage.rawCameraRGB);
cd(curDir);

%% Read location file, get data, and make a diagnostic image of where we are extracting from.
% Locations were extracted from .JPG produced by camera.  These need to be divided by 2 to deal
% with our demosaicing procedure.  A little tweaking by eye to get stuff centered after initial read.
locations = ReadStructsFromText([filenameReal '_Bots.txt']);
botsXYZ = zeros(3,24);
if (CHECK)
    ig = realImage;
    ig = ig ./ max(ig(:));
end
for j = 1:24
    if (~strcmp(locations(j).Location,sprintf('MCC%d',j)))
        error('Unexpected MCC location name in external standard image text file');
    end
    mccCenterRow = round(locations(j).Y/2);
    mccCenterCol = round(locations(j).X/2);
    mccRowRegion = mccCenterRow-mccBlkSize:mccCenterRow+mccBlkSize;
    mccColRegion = mccCenterCol-mccBlkSize:mccCenterCol+mccBlkSize;
    mccCenterX = mean(mean(realImage(mccRowRegion,mccColRegion,1)));
    mccCenterY = mean(mean(realImage(mccRowRegion,mccColRegion,2)));
    mccCenterZ = mean(mean(realImage(mccRowRegion,mccColRegion,3)));
    botsXYZ(:,j) = [mccCenterX mccCenterY mccCenterZ]';
    if (CHECK)
        ig(mccRowRegion,mccColRegion,:) = 1;
    end
end
if (CHECK)
    igFig = figure;
    figure(igFig);clf;imagesc(ig);
    cd(curDir);
    imwrite(ig,['AuxPhotoRegions_' whichBots '.jpg'],'jpg');
end

%% Repeat for standard camera
cd([idbPath '/CAMERACOMPARISON/PhillyCam/Color']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% Get the image data
filenameReal = ['DSC_' whichPhilly];
if (CHECK == 1)
end
fprintf('Processing file %s\n',filenameReal);

% Get exposure duration, etc
imageInfoPhilly = GetNEFInfo(filenameReal);

% Get exposure duration, etc
fprintf('\t\tCamera: %s\n',imageInfoPhilly.whichCamera);
fprintf('\t\tExposure %g\n',imageInfoPhilly.exposure);
fprintf('\t\tfStop %g\n',imageInfoPhilly.fStop);
fprintf('\t\tISO %g\n',imageInfoPhilly.ISO);
if (~strcmp(imageInfoPhilly.whichCamera,'standard'))
    error('This should be standard camera\n');
end
if (imageInfoPhilly.exposure ~= imageInfoBots.exposure)
    error('Two cameras should use same exposure');
end
if (imageInfoPhilly.fStop ~= imageInfoBots.fStop)
    error('Two cameras should use same fStop');
end
if (imageInfoPhilly.ISO ~= imageInfoBots.ISO)
    error('Two cameras should use same ISO');
end

% Read in raw image, dark subtract, and truncate negative values
load([filenameReal '.raw.mat']);
CamCal = LoadCamCal(imageInfoPhilly.whichCamera);
realImage = DarkCorrect(CamCal,imageInfoPhilly,theImage.rawCameraRGB);
cd(curDir);

%% Read location file, get data, and make a diagnostic image of where we are extracting from.
% Locations were extracted from .JPG produced by camera.  These need to be divided by 2 to deal
% with our demosaicing procedure.  A little tweaking by eye to get stuff centered after initial read.
locations = ReadStructsFromText([filenameReal '_Philly.txt']);
phillyXYZ = zeros(3,24);
if (CHECK)
    ig = realImage;
    ig = ig ./ max(ig(:));
end
for j = 1:24
    if (~strcmp(locations(j).Location,sprintf('MCC%d',j)))
        error('Unexpected MCC location name in external standard image text file');
    end
    mccCenterRow = round(locations(j).Y/2);
    mccCenterCol = round(locations(j).X/2);
    mccRowRegion = mccCenterRow-mccBlkSize:mccCenterRow+mccBlkSize;
    mccColRegion = mccCenterCol-mccBlkSize:mccCenterCol+mccBlkSize;
    mccCenterX = mean(mean(realImage(mccRowRegion,mccColRegion,1)));
    mccCenterY = mean(mean(realImage(mccRowRegion,mccColRegion,2)));
    mccCenterZ = mean(mean(realImage(mccRowRegion,mccColRegion,3)));
    phillyXYZ(:,j) = [mccCenterX mccCenterY mccCenterZ]';
    if (CHECK)
        ig(mccRowRegion,mccColRegion,:) = 1;
    end
end
if (CHECK)
    igFig = figure;
    figure(igFig);clf;imagesc(ig);
    cd(curDir);
    imwrite(ig,['StandardPhotoRegions_' whichPhilly '.jpg'],'jpg');
end

%% Now plot one versus the other and find ratio
f3 = figure; clf
position = get(f3,'Position');
position(3) = 680; position(4) = 680;
set(f3,'Position',position);

warning('off','curvefit:fit:noStartPoint');


fitrange{1} = find((botsXYZ(1,:) >= fitLowResp) & (botsXYZ(1,:) < fitHighResp) & ...
    (phillyXYZ(1,:) >= fitLowResp) & (phillyXYZ(1,:) < fitHighResp));
fitrange{2} = find((botsXYZ(2,:) >= fitLowResp) & (botsXYZ(2,:) < fitHighResp) & ...
    (phillyXYZ(2,:) >= fitLowResp) & (phillyXYZ(2,:) < fitHighResp));
fitrange{3} = find((botsXYZ(3,:) >= fitLowResp) & (botsXYZ(3,:) < fitHighResp) & ...
    (phillyXYZ(3,:) >= fitLowResp) & (phillyXYZ(3,:) < fitHighResp));
stanLinearityTable = phillyXYZ';
auxLinearityTable = botsXYZ';

ft = '1*x+b';
pfall = fit(log10([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},2) ; stanLinearityTable(fitrange{3},3)]), ...
	log10([auxLinearityTable(fitrange{1},1) ; auxLinearityTable(fitrange{2},2); auxLinearityTable(fitrange{3},3)]),ft);
overallSlope = 10.^pfall.b;

clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope 1 through the points of interest and add to plot
    pf = fit(log10(stanLinearityTable(fitrange{i},i)),log10(auxLinearityTable(fitrange{i},i)),ft);
    freeslope(i) = 10.^pf.b;

    figure(f3);
    loglog(stanLinearityTable(:,i), auxLinearityTable(:,i),'o','color',clx{i}); hold on
    loglog(stanLinearityTable(fitrange{i},i), auxLinearityTable(fitrange{i},i),'o','MarkerFaceColor',clx{i},'color',clx{i});
end
loglog([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},2) ; stanLinearityTable(fitrange{3},3)], ...
    10.^feval(pfall,log10([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},2) ; stanLinearityTable(fitrange{3},3)])), ...
    '-','color','k','linewidth',1.3);

% Show a slope of 0.83 too
if (0)
    pfall.b = log10(0.83);
    loglog([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},2) ; stanLinearityTable(fitrange{3},3)], ...
        10.^feval(pfall,log10([stanLinearityTable(fitrange{1},1) ; stanLinearityTable(fitrange{2},2) ; stanLinearityTable(fitrange{3},3)])), ...
        ':','color','r','linewidth',1);
end

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
FigureSave(['CameraComparisonMCCFig_' whichPhilly '_' whichBots '.pdf'],f3,'pdf');

warning('on','curvefit:fit:noStartPoint');

return

