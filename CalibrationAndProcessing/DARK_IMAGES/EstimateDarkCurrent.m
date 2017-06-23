% EstimateDarkCurrent
%
% We took a bunch of images in the dark at various exposure times and 
% process these to make a table of RGB dark values for each.
%
% Note that the distribution of dark values across pixels is very skewed:
% some pixels are 'hot' and produce very large values. So it makes sense to
% estimat the dark current either by taking the median (which will neglect
% these outliers) or by taking a mean over all pixels that are BELOW some
% threshold (choosing that as 200 raw values, 99.9% of the pixels are
% below). This method will compute both these measures, but will store as
% dark current table only the median values.
%
% Parameters:
%   idbPath -- path to the root of the calibration image database
%
% 05/18/10  dhb, gt, pg  Wrote it -- we're on a roll now.
% 07/08/10  gt           Updated to compute the means and make the plots.
% 08/20/10  gt           Modified dark table by hand so that below 1s exposure it has fixed darkcurrent values (4,8,0);
% 11/09/10  dhb          Add default path to save typing it over and over.  Fix curDir glitch at end of first section.
% 11/09/10  dhb          Hand modification had R and B reversed from the pattern of the raw data.  Decided to take 
%                        the mean of the medians and use this for exposures <= 1 second.  (It used to be strictly less
%                        than, but to my eye the 1s data is like the less than 1s data.
%           dhb          Save the raw dark table too.
%           dhb          Put in check for all data at ISO 400.
%           dhb          Fix rgb color order in the plot, and save as PDF.

function EstimateDarkCurrent(idbPath)

%% Set default path, corresponds to our setup
if (nargin < 1 || isempty(idbPath))
    idbPath = '../../Images/calibration';
end

%% Clear out old stuff on general principles;
curDir = pwd;

%% First process the ascending wavelength data
cd([idbPath '/DARK_IMAGES']);
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
    if (imageInfoReal.ISO ~= 400)
        error('All data should be at ISO 400');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Store exposure
    darkTable(f,1) = imageInfoReal.exposure;
    
    % Read in raw image and dark subtract
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    
    % Get spectral sensitivity be subtracting dark from real and normalizing by
    % power.

    for k = 1:3
        temp = realImage(:,:,k);
        darkTable(f,k+1) = median(temp(:));
        darkTableMeanBelow200(f,k) = mean(temp(temp<200));
        darkTableFracBelow200(f,k) = sum(temp(:)<200)/numel(temp);
    end
end
cd(curDir);

%% Modify the dark table to be more smooth
disp('Modifying the dark current table by hand for values below 1s exposure');
darkTableRaw = darkTable;
ix = (darkTable(:,1)<=1);
fprintf('Maximum median for exposures less than 1s: %g, %g, %g\n',max(darkTableRaw(ix,2)),max(darkTableRaw(ix,3)),max(darkTableRaw(ix,4)));
fprintf('Dark table fraction below 200: %g, %g, %g\n',darkTableFracBelow200(f,1),darkTableFracBelow200(f,2),darkTableFracBelow200(f,3));
darkTable(ix,2) = round(mean(darkTableRaw(ix,2)));
darkTable(ix,3) = round(mean(darkTableRaw(ix,3)));
darkTable(ix,4) = round(mean(darkTableRaw(ix,4)));

%% Save as text file in case we want a table in the paper
save DarkTable.txt darkTable -ascii

%% Make plots
figure; clf;
subplot(1,2,1); 
semilogx(darkTable(:,1),darkTable(:,2),'r','linewidth',1.3);
hold on
semilogx(darkTable(:,1),darkTable(:,3),'g','linewidth',1.3);
semilogx(darkTable(:,1),darkTable(:,4),'b','linewidth',1.3);
axis square;
set(gca,'fontsize',14); xlabel('Exposure (s)','fontsize',14); ylabel('Median dark current','fontsize',14);
%title('All pixels (100%), f=22, ISO=400','fontsize',14);

subplot(1,2,2); 
semilogx(darkTable(:,1),darkTableMeanBelow200(:,1),'r','linewidth',1.3);
hold on
semilogx(darkTable(:,1),darkTableMeanBelow200(:,2),'g','linewidth',1.3);
semilogx(darkTable(:,1),darkTableMeanBelow200(:,3),'b','linewidth',1.3);
axis square;
set(gca,'fontsize',14); xlabel('Exposure (s)','fontsize',14); ylabel('Mean dark current','fontsize',14);
%title('All pixels < 200 (99.9% pixels), f=22, ISO=400','fontsize',14);
FigureSave('DarkCurrent.pdf',gcf,'pdf');

%% Little check on ISO effect
cd([idbPath '/DARK_IMAGES_ISOVARY']);

% List NEF files of given directory
fileSpec = ['*.NEF'];
theFiles = dir(fileSpec);

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
        %error('All data should be at ISO 1000');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Read in raw image and dark subtract
    load([filenameReal '.raw.mat']); realImage = theImage.rawCameraRGB;
    
    % Get spectral sensitivity be subtracting dark from real and normalizing by
    % power.
    for k = 1:3
        temp = realImage(:,:,k);
        darkValue(k) = median(temp(:));
        meanDarkValue(k) = mean(temp(temp < 100));
    end
    fprintf('\t\tMedian dark triplet %g, %g, %g\n',darkValue(1),darkValue(2),darkValue(3));
    fprintf('\t\tMean dark triplet (pixels < 100) %g, %g, %g\n',meanDarkValue(1),meanDarkValue(2),meanDarkValue(3));

end
cd(curDir)

%% Sort by exposure and save
[nil,index] = sort(darkTable(:,1));
darkTable = darkTable(index,:);
save darkTable darkTable darkTableRaw




