% EstimateISOLinearity
%
% We took a bunch of images of the white standard at different ISO settings
% to check that the camera responds linearly.
%
% Parameters:
%   idbPath -- path to the root of the calibration image database
%
% 07/08/10  gt  Wrote it.
% 11/13/10  dhb Same basic set of changes as to other linearity check programs.

function EstimateISOLinearity(idbPath)

% Close
close all;

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Set current dir
curDir = pwd;

%% If 1, plot where in the image the signal is being extracted from
CHECK = 1;

%% These are the coordinates of the white standard in the exposure
%% linearity photo set; nQ x nQ is the size of the region in pixels
nQ          = 300;
crop_coords = [600 700];


%% First directory
cd([idbPath '/ISO_LINEARITY_EXPOSURE1']);
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
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    if (imageInfoReal.exposure ~= 0.008)
        error('All data should be for exposure 0.008');
    end
    if (imageInfoReal.fStop ~= 2.8)
        error('All data should be for fStop 2.8');
    end
    
    % Load the camera data for the image
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    
    % Read in raw image, dark subtract, and truncate negative values
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);
    
    % Diagnostic plot of where we are extracting from
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(crop_coords(1):2:crop_coords(1)+nQ, crop_coords(2):2:crop_coords(2)+nQ,:) = 1;
        ig(crop_coords(1)+1:2:crop_coords(1)+nQ, crop_coords(2)+1:2:crop_coords(2)+nQ,:) = 0;
        figure(igFig);clf;imagesc(ig);
        if (f == 4)
            cd(curDir);
            imwrite(ig,'ExposureRegions.jpg','jpg');
            cd([idbPath '/ISO_LINEARITY_EXPOSURE1']);
        end
    end
    
    % Select plaquettes from the image
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
    end
    
    % build the linearity table for different exposure and fstop
    linearityTable(f,1) = imageInfoReal.ISO;
    linearityTable(f,2:4)= croppedImage; 
    
end


%% Second directory
cd(curDir);
cd([idbPath '/ISO_LINEARITY_EXPOSURE2']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

% Loop over files and analyze
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
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    if (imageInfoReal.exposure ~= 0.004)
        error('All data should be for exposure 0.004');
    end
    if (imageInfoReal.fStop ~= 2.8)
        error('All data should be for fStop 2.8');
    end
    
    % Load the camera data for the image
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    
    % Read in raw image, dark subtract, and truncate negative values
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    realImage = DarkCorrect(CamCal,imageInfoReal,theImage.rawCameraRGB);
    
    % Diagnostic plot of where we are extracting from
    if (CHECK==1)
        ig = realImage;
        ig = ig ./ max(ig(:));
        ig(crop_coords(1):2:crop_coords(1)+nQ, crop_coords(2):2:crop_coords(2)+nQ,:) = 1;
        ig(crop_coords(1)+1:2:crop_coords(1)+nQ, crop_coords(2)+1:2:crop_coords(2)+nQ,:) = 0;
        figure(igFig);clf;imagesc(ig)
        if (f == 4)
            cd(curDir);
            imwrite(ig,'ExposureRegions2.jpg','jpg');
            cd([idbPath '/ISO_LINEARITY_EXPOSURE2']);
        end
    end
    
    % Select plaquettes from the image
    for i=1:3,
        croppedImage(i)  = mean(mean(realImage(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ, i)));
    end
    
    % build the linearity table for different exposure and fstop
    linearityTable1(f,1) = imageInfoReal.ISO;
    linearityTable1(f,2:4)= croppedImage; 
    
end
cd(curDir);

% Plot
warning('off','curvefit:fit:noStartPoint');
figure;
fitLowResp = 50;
fitHighResp = 16100;
fitrange{1} = find((linearityTable(:,2) >= fitLowResp) & (linearityTable(:,2) < fitHighResp));
fitrange{2} = find((linearityTable(:,3) >= fitLowResp) & (linearityTable(:,3) < fitHighResp));
fitrange{3} = find((linearityTable(:,4) >= fitLowResp) & (linearityTable(:,4) < fitHighResp));
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(linearityTable(fitrange{i},1)), log10(linearityTable(fitrange{i},1+i)),ft);

    loglog(linearityTable(:,1), linearityTable(:,1+i),'o','color',clx{i});hold on;
    loglog(linearityTable(fitrange{i},1), linearityTable(fitrange{i},1+i),'o','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(linearityTable(fitrange{i},1), 10.^(feval(pf, log10(linearityTable(fitrange{i},1)))),'-','color',clx{i},'linewidth',1.3);
    
    % Compute prediction error as mean of absolute log unit error
    meanLogError(1,i) = mean(abs(log10(linearityTable(fitrange{i},1+i)) - (feval(pf, log10(linearityTable(fitrange{i},1))))));
    maxLogError(1,i) = max(abs(log10(linearityTable(fitrange{i},1+i)) - (feval(pf, log10(linearityTable(fitrange{i},1))))));
   
    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(linearityTable(fitrange{i},1)), log10(linearityTable(fitrange{i},1+i)),ft1);
    freeslope(1,i) = pf1.b;
end

fitrange{1} = find((linearityTable1(:,2) >= fitLowResp) & (linearityTable1(:,2) < fitHighResp));
fitrange{2} = find((linearityTable1(:,3) >= fitLowResp) & (linearityTable1(:,3) < fitHighResp));
fitrange{3} = find((linearityTable1(:,4) >= fitLowResp) & (linearityTable1(:,4) < fitHighResp));
clx={'r';'g';'b'};
for i=1:3
    % Fit a line of slope -2 through the points of interest and add to plot
    ft = fittype('1*x+a');
    pf = fit(log10(linearityTable1(fitrange{i},1)), log10(linearityTable1(fitrange{i},1+i)),ft);

    loglog(linearityTable1(:,1), linearityTable1(:,1+i),'s','color',clx{i});hold on;
    loglog(linearityTable1(fitrange{i},1), linearityTable1(fitrange{i},1+i),'s','MarkerFaceColor',clx{i},'color',clx{i});
    loglog(linearityTable1(fitrange{i},1), 10.^(feval(pf, log10(linearityTable1(fitrange{i},1)))),'--','color',clx{i},'linewidth',1.3);
    
    % Compute prediction error as mean of absolute log unit error
    meanLogError(2,i) = mean(abs(log10(linearityTable1(fitrange{i},1+i)) - (feval(pf, log10(linearityTable1(fitrange{i},1))))));
    maxLogError(2,i) = max(abs(log10(linearityTable1(fitrange{i},1+i)) - (feval(pf, log10(linearityTable1(fitrange{i},1))))));

   
    % Determine slope if it is allowed to go fre
    ft1 = fittype('b*x+a');
    pf1 = fit(log10(linearityTable1(fitrange{i},1)), log10(linearityTable1(fitrange{i},1+i)),ft1);
    freeslope(2,i) = pf1.b;
end

axis square;
set(gca,'fontsize',14);
xlabel('ISO','fontsize',14);
ylabel('Raw camera RGB (dark subtracted)','fontsize',14);
xlim([100 10000]); ylim([100 100000]);
FigureSave('ISOLinearity.pdf',gcf,'pdf');
warning('on','curvefit:fit:noStartPoint');


%% Dump info into a text file free slopes
fid = fopen('ISOLinearityInfo.txt','w');
for fidx = 1:2
    for i = 1:3
        fprintf(fid,'Free fit slope for channel %d, set %d, is %0.2f\n',i,fidx,freeslope(fidx,i));
        fprintf(fid,'Mean absolute log10 error for channel %d, set %d, is %0.3f\n',i,fidx,meanLogError(fidx,i));
        fprintf(fid,'Max absolute log10 error for channel %d, set %d, is %0.3f\n',i,fidx,maxLogError(fidx,i));
    end
    fprintf(fid,'\n');
end
fclose(fid);

