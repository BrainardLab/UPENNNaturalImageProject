% EstimateExposureLinearity
%
% We took a bunch of images in the dark at various exposure times and 
% checked that the response of the camera is proportional to the exposure
% time.
%
% Parameters:
%   idbPath -- the path to the root of the image calibration database
%
% 07/08/10  gt  Wrote it.
% 11/11/10  dhb Default directory added so you don't have to remember it.
% 11/12/10  dhb Add aperature plot.
%           dhb Fit in log10
%           dhb Select fit points based on response.
%           dhb Print minimum aperture

function EstimateExposureLinearity(idbPath)

% Close
close all;

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Clear out old stuff on general principles;
curDir = pwd;

%% If 1, plot where in the image the signal is being extracted from
CHECK = 1;

%% These are the coordinates of the white standard in the exposure
%% linearity photo set; nQ x nQ is the size of the region in pixels
nQ          = 100;
crop_coords = [450 700];

%% Load the data
cd([idbPath '/EXPOSURE_LINEARITY']);
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
    if (imageInfoReal.ISO ~= 200)
        error('All data should be at ISO 200');
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
        if (f == 35)
            cd(curDir);
            imwrite(ig,'ExposureRegions.jpg','jpg');
            cd([idbPath '/EXPOSURE_LINEARITY']);
        end
    end
    
    % Select plaquettes from the image
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
    end
    
    % Build the linearity table for different exposure and fstop
    linearityTable(f,1) = imageInfoReal.exposure;
    linearityTable(f,2) = imageInfoReal.fStop;
    linearityTable(f,3:5)= croppedImage; 
    
end
cd(curDir);

%% Make the basic exposure linearity plot 
% We have a note from elsewhere:
%   Images 56, 112, and 168 correspond to the "Bulb" setting on the camera and are not meaningful.
%   The exif header does give an exposure time for these images, however, and it may be that
%   the camera records and stores the exposure duration for these cases.
fstops = [5 11 16];
for i=1:numel(fstops),
    idx = find(linearityTable(:,2)==fstops(i));
    newLinearityTable{i} = linearityTable(idx,[1 3 4 5]);
end

% Three panel plot, one for each f-number.  This provides nice range
% extension.
f2 = figure; clf
position = get(f2,'Position');
position(3) = 1624; position(4) = 680;
set(f2,'Position',position);
fitLowResp = 50;
fitHighResp = 16100;
warning('off','curvefit:fit:noStartPoint');

fidx = 3;
fitrange{1} = find((newLinearityTable{fidx}(:,2) >= fitLowResp) & (newLinearityTable{fidx}(:,2) < fitHighResp));
fitrange{2} = find((newLinearityTable{fidx}(:,3) >= fitLowResp) & (newLinearityTable{fidx}(:,3) < fitHighResp));
fitrange{3} = find((newLinearityTable{fidx}(:,4) >= fitLowResp) & (newLinearityTable{fidx}(:,4) < fitHighResp));
subplot(1,3,3);
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft);

    figure(f2);
    loglog(newLinearityTable{fidx}(:,1), newLinearityTable{fidx}(:,1+i),'o','color',clx{i});hold on;
    loglog(newLinearityTable{fidx}(fitrange{i},1), newLinearityTable{fidx}(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(newLinearityTable{fidx}(fitrange{i},1), 10.^(feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
        
    % Compute prediction error as mean of absolute log unit error
    meanLogError(fidx,i) = mean(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));
    maxLogError(fidx,i) = max(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));
   
    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end

% Compute minimum exposure analyzed
minExposure(fidx) = min([newLinearityTable{fidx}(fitrange{1},1) ; newLinearityTable{fidx}(fitrange{2},1) ; newLinearityTable{fidx}(fitrange{3},1)]);

axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
title(['f ' num2str(fstops(fidx))],'fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);

fidx = 2;
fitrange{1} = find((newLinearityTable{fidx}(:,2) >= fitLowResp) & (newLinearityTable{fidx}(:,2) < fitHighResp));
fitrange{2} = find((newLinearityTable{fidx}(:,3) >= fitLowResp) & (newLinearityTable{fidx}(:,3) < fitHighResp));
fitrange{3} = find((newLinearityTable{fidx}(:,4) >= fitLowResp) & (newLinearityTable{fidx}(:,4) < fitHighResp));
subplot(1,3,2);
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft);

    figure(f2);
    loglog(newLinearityTable{fidx}(:,1), newLinearityTable{fidx}(:,1+i),'o','color',clx{i});hold on;
    loglog(newLinearityTable{fidx}(fitrange{i},1), newLinearityTable{fidx}(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(newLinearityTable{fidx}(fitrange{i},1), 10.^(feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    % Compute prediction error as mean of absolute log unit error
    meanLogError(fidx,i) = mean(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));
    maxLogError(fidx,i) = max(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));

    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end

% Compute minimum exposure analyzed
minExposure(fidx) = min([newLinearityTable{fidx}(fitrange{1},1) ; newLinearityTable{fidx}(fitrange{2},1) ; newLinearityTable{fidx}(fitrange{3},1)]);

axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
title(['f ' num2str(fstops(fidx))],'fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);

fidx = 1;
fitrange{1} = find((newLinearityTable{fidx}(:,2) >= fitLowResp) & (newLinearityTable{fidx}(:,2) < fitHighResp));
fitrange{2} = find((newLinearityTable{fidx}(:,3) >= fitLowResp) & (newLinearityTable{fidx}(:,3) < fitHighResp));
fitrange{3} = find((newLinearityTable{fidx}(:,4) >= fitLowResp) & (newLinearityTable{fidx}(:,4) < fitHighResp));
subplot(1,3,1);
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft);

    figure(f2);
    loglog(newLinearityTable{fidx}(:,1), newLinearityTable{fidx}(:,1+i),'o','color',clx{i});hold on;
    loglog(newLinearityTable{fidx}(fitrange{i},1), newLinearityTable{fidx}(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(newLinearityTable{fidx}(fitrange{i},1), 10.^(feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
        
    % Compute prediction error as mean of absolute log unit error
    meanLogError(fidx,i) = mean(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));
    maxLogError(fidx,i) = max(abs(log10(newLinearityTable{fidx}(fitrange{i},1+i)) - (feval(pf, log10(newLinearityTable{fidx}(fitrange{i},1))))));

    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(newLinearityTable{fidx}(fitrange{i},1)), log10(newLinearityTable{fidx}(fitrange{i},1+i)),ft1);
    freeslope(i) = pf1.b;
end

% Compute minimum exposure analyzed
minExposure(fidx) = min([newLinearityTable{fidx}(fitrange{1},1) ; newLinearityTable{fidx}(fitrange{2},1) ; newLinearityTable{fidx}(fitrange{3},1)]);

axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
title(['f ' num2str(fstops(fidx))],'fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);

% Save the figure
FigureSave('ExposureFig.pdf',f2,'pdf');
warning('on','curvefit:fit:noStartPoint');

%% Plot of overlay for different apertures.
% This plot isn't perfected for the paper, since we don't show it.  But the
% good overlay for all three color channels confirms our basic aperature analysis.
apertureTable{2} = newLinearityTable{2};
apertureTable{1} = newLinearityTable{1};
apertureTable{1}(:,2:4) = ((fstops(1)/fstops(2))^2)*apertureTable{1}(:,2:4);
apertureTable{3} = newLinearityTable{3};
apertureTable{3}(:,2:4) = ((fstops(3)/fstops(2))^2)*apertureTable{3}(:,2:4);

% Three panel plot, one for each f-number.  This provides nice range
% extension.
f3 = figure; clf
position = get(f3,'Position');
position(3) = 1624; position(4) = 680;
set(f3,'Position',position);

subplot(1,3,1);
clx={'r';'g';'b'};
loglog(apertureTable{1}(:,1), apertureTable{1}(:,2),'o','MarkerFaceColor',clx{1},'color',clx{1});hold on;
loglog(apertureTable{2}(:,1), apertureTable{2}(:,2),'o','MarkerFaceColor',clx{2},'color',clx{2});hold on;
loglog(apertureTable{3}(:,1), apertureTable{3}(:,2),'o','MarkerFaceColor',clx{3},'color',clx{3});hold on;
axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Scaled to f11 RGB (dark subtracted)','fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);
title('Red channel, f 5, 11, and 16');

subplot(1,3,2);
clx={'r';'g';'b'};
loglog(apertureTable{1}(:,1), apertureTable{1}(:,3),'o','MarkerFaceColor',clx{1},'color',clx{1});hold on;
loglog(apertureTable{2}(:,1), apertureTable{2}(:,3),'o','MarkerFaceColor',clx{2},'color',clx{2});hold on;
loglog(apertureTable{3}(:,1), apertureTable{3}(:,3),'o','MarkerFaceColor',clx{3},'color',clx{3});hold on;
axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Scaled to f11 RGB (dark subtracted)','fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);
title('Green channel, f 5, 11, and 16');

subplot(1,3,3);
clx={'r';'g';'b'};
loglog(apertureTable{1}(:,1), apertureTable{1}(:,4),'o','MarkerFaceColor',clx{1},'color',clx{1});hold on;
loglog(apertureTable{2}(:,1), apertureTable{2}(:,4),'o','MarkerFaceColor',clx{2},'color',clx{2});hold on;
loglog(apertureTable{3}(:,1), apertureTable{3}(:,4),'o','MarkerFaceColor',clx{3},'color',clx{3});hold on;
axis square;
set(gca,'fontsize',14);
xlabel('Exposure time (secs)','fontsize',14);
ylabel('Scaled to f11 RGB (dark subtracted)','fontsize',14);
xlim([10^-4 100]); ylim([1 100000]);
title('Blue channel, f 5, 11, and 16');

FigureSave('ExposureApertureCompare.pdf',f3,'pdf');

%% Dump info into a text file free slopes
fid = fopen('ExposureLinearityInfo.txt','w');
for fidx = 3:-1:1
    for i = 1:3
        fprintf(fid,'Free fit slope for channel %d f-%g is %0.2f\n',i,fstops(fidx),freeslope(i));
        fprintf(fid,'Mean absolute log10 error for channel %d f-%g is %0.3f\n',i,fstops(fidx),meanLogError(fidx,i));
        fprintf(fid,'Max absolute log10 error for channel %d f-%g is %0.3f\n',i,fstops(fidx),maxLogError(fidx,i));
    end
    fprintf(fid,'Min exposure analyzed for f-%g was %f\n',fstops(fidx),minExposure(fidx));
    fprintf(fid,'\n');
end
fclose(fid);
