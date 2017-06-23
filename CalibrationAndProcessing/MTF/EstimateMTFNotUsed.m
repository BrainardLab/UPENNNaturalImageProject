% EstimateMTFNotUsed
%
% Process the grating images to estimate LMS/RGB MTF.
% 
% This was set up to take a quick look at the data we 
% collected in December 2010.  In the end, we probably
% don't need this data, although I'm keeping this routine
% around because it contains some of the relevant coordinates
% we identified.
%
% Parameters:
%   idbPath -- the path to the root of the calibration image database
%
% 12/22/10  dhb Merged in reading from raw images with older code.

function EstimateMTFNotUsed(idbPath)

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

%% These are the coordinates of the grating region in each image.
% These were obtained by hand for each image in Photoshop.
TWOPOINT = 1;
if (TWOPOINT)
    coords(1).name = '0002';
    coords(1).rows = 526:526+20;
    coords(1).cols = 606:606+75;
    
    coords(2).name = '0028';
    coords(2).rows = 463:463+9;
    coords(2).cols = 744:744+38;
else
    for increment=1:12
        coords(increment).name='0002';
    end;
    
    coords(1).rows = 424:424+7;
    coords(1).cols = 620:620+19;
    
    coords(2).rows = 434:434+7;
    coords(2).cols = 620:620+19;
    
    coords(3).rows = 449:449+9;
    coords(3).cols = 621:621+17;
    
    coords(4).rows = 460:460+9;
    coords(4).cols = 621:621+16;
    
    coords(5).rows = 475:475+8;
    coords(5).cols = 623:623+13;
    
    coords(6).rows = 488:488+5;
    coords(6).cols = 623:623+13;
    
    coords(7).rows = 499:499+7;
    coords(7).cols = 624:624+12;
    
    coords(8).rows = 513:513+5;
    coords(8).cols = 624:624+11;
    
    coords(9).rows = 418:418+6;
    coords(9).cols = 653:653+10;
    
    coords(10).rows = 428:428+6;
    coords(10).cols = 653:653+10;
    
    coords(11).rows = 438:438+5;
    coords(11).cols = 654:654+8;
    
    coords(12).rows = 449:449+3;
    coords(12).cols = 654:654+8;
end

%% Load the photos.
cd([idbPath '/MTF/NEW2010']);
theDirectory = pwd;
fprintf('Image directory is %s\n',theDirectory);

% Loop over files and analyze
if (CHECK)
    igFig = figure; clf;
    igFig2 = figure; clf;
    igFig3 = figure; clf;
end
for f = 1:length(coords)
    filenameReal = sprintf('DSC_%s',coords(f).name);
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
    if (imageInfoReal.fStop ~= 5.6)
        error('All data should be at f 5.6');
    end
    if (~strcmp(imageInfoReal.whichCamera,'standard'))
        error('This should be standard camera\n');
    end
    
    % Read in image, dark subtract, and truncate negative values
    load([filenameReal '_LMS.mat']);
    load([filenameReal '_RGB.mat']);
    
    % Check extracted location
    if (CHECK)
        ig = LMS_Image;
        ig = ig ./ max(ig(:));
        ig(ig<0) = 0;
        ig(coords(f).rows,coords(f).cols,:) = 1;
        figure(igFig); imagesc(ig);
    end
    
    % Rotate so that we match convention on which rest of code is written.
    for k = 1:3
        LMS_Image1(:,:,k) = rot90(LMS_Image(coords(f).rows,coords(f).cols,k));
        RGB_Image1(:,:,k) = rot90(RGB_Image(coords(f).rows,coords(f).cols,k));
    end
    LMS_Image = LMS_Image1; clear LMS_Image1;
    RGB_Image = RGB_Image1; clear RGB_Image1;
    
    CamCal = LoadCamCal(imageInfoReal.whichCamera);
    RGB_Image = DarkCorrect(CamCal,imageInfoReal,RGB_Image);
    
    % Extract image data and scale to standardized representation
    whichCols = 2:4;
    for k=1:3,
        lmsImage{f,k} = LMS_Image(:,whichCols,k);
        rgbImage{f,k} = RGB_Image(:,whichCols,k);
        lmsGrating{f,k} = sum(lmsImage{f,k},2);
        rgbGrating{f,k} = sum(rgbImage{f,k},2);
    end
    
    % Response range check
    if (any(RGB_Image < fitLowResp) | any(RGB_Image > fitHighResp))
        fprintf('\tWARNING: Some data out of good response range\n');
    end
    
    % More check of location and form of data
    if (CHECK) 
        ig2 = RGB_Image;
        ig2(ig2<0) = 0;
        for k = 1:3
            ig2(:,:,k) = 0.5*ig2(:,:,k)/mean(mean(ig2(:,:,k)));
        end
        ig2(ig2>1) = 1;
        figure(igFig2);
        subplot(1,3,1); imshow(ig2(:,:,1));
        subplot(1,3,2); imshow(ig2(:,:,2));
        subplot(1,3,3); imshow(ig2(:,:,3));
        
        figure(igFig3); clf;
        plot(0.5*rgbGrating{f,1}/mean(rgbGrating{f,1}),'r'); hold on
        plot(0.5*rgbGrating{f,2}/mean(rgbGrating{f,2}),'g');
        plot(0.5*rgbGrating{f,3}/mean(rgbGrating{f,3}),'b');
        hold off
        
        drawnow;
        %f
        %pause;
    end
    
end
cd(curDir)

%% Get power and sf for every image
for i=1:length(coords)
    [lsf(i),lamp(i)] = FindMaxAmp(lmsImage{i,1});
    [msf(i),mamp(i)] = FindMaxAmp(lmsImage{i,2});
    [ssf(i),samp(i)] = FindMaxAmp(lmsImage{i,3});
    
    [rsf(i),ramp(i)] = FindMaxAmp(rgbImage{i,1},1);
    [gsf(i),gamp(i)] = FindMaxAmp(rgbImage{i,2});
    [bsf(i),bamp(i)] = FindMaxAmp(rgbImage{i,3});
end

% LMS Figure
mtfFigLMS = figure; clf; hold on;
plot(lsf,lamp,'ro','markersize',8,'markerfacecolor','r');
plot(msf,mamp,'go','markersize',8,'markerfacecolor','g');
plot(ssf,samp,'bo','markersize',8,'markerfacecolor','b');

set(gca,'xlim',[0.05 0.5]);
set(gca,'ylim',[0 1.2]);
set(gca,'xtick',0:0.10:0.5);
xlabel('Spatial frequency (cycles/pixel)','fontsize',14);
set(gca,'fontsize',14);
ylabel('LMS Amplitude','fontsize',14);
axis square
if (TWOPOINT)
    FigureSave('MTFFigureLMSNotUsedTwo.pdf',mtfFigLMS,'pdf');
else
    FigureSave('MTFFigureLMSNotUsedOne.pdf',mtfFigLMS,'pdf');
end

% RGB Figure
mtfFigRGB = figure; clf; hold on;
plot(rsf,ramp,'ro','markersize',8,'markerfacecolor','r');
plot(gsf,gamp,'go','markersize',8,'markerfacecolor','g');
plot(bsf,bamp,'bo','markersize',8,'markerfacecolor','b');

set(gca,'xlim',[0.05 0.5]);
set(gca,'ylim',[0 1.2]);
set(gca,'xtick',0:0.10:0.5);
xlabel('Spatial frequency (cycles/pixel)','fontsize',14);
set(gca,'fontsize',14);
ylabel('RGB Amplitude','fontsize',14);
axis square
if (TWOPOINT)
    FigureSave('MTFFigureRGBNotUsedTwo.pdf',mtfFigRGB,'pdf');
else
    FigureSave('MTFFigureRGBNotUsedOne.pdf',mtfFigRGB,'pdf');
end

end


%% This is based on Pat's code, that deals with power
% splatter by checking various croppings.  The idea
% is to find the cropping that comes closest to an integer
% number of cycles.
%
% Spatial frequency returned in cycles/pixel
%
% Also assumes that gratings really are horizontal
%
% 12/23/10  dhb  Converted to work on 1D fft, since we were already
%                assuming horizontal input.
%           dhb  Clean up
%           dhb  Fix bug in returned sf

function [sf,amp] = FindMaxAmp(theImage,PLOT)

% Diagnostic figure
if (nargin < 2 || isempty(PLOT))
    PLOT = 0;
end
PLOTALL = 0;
if (PLOT)
    theFig = figure; clf;
end
% Assume horizontal grating, sum over columns to make]
% it one dimensional
theImage = sum(theImage,2);   

% Search over cropping to find maximum power at fundamental
nPhases = min([20 length(theImage)-5]);
for phase=1:nPhases
    
    % Crop image and take fft.  Normalize so that
    % dc component comes out as 1 like it should
    croppedImage=theImage(phase:size(theImage,1),:);    
    croppedImage=croppedImage./mean(croppedImage(:)); 	
    theFFT{phase}=fft(croppedImage);                           
    theFFT{phase}=theFFT{phase}./(length(theFFT{phase}(:)));                
    
    % Get power and amplitude spectra
    powerSpectrum=theFFT{phase}.*conj(theFFT{phase});                 
    ampSpectrum=sqrt(powerSpectrum);                    
    
    % Find non-zero freq with max response, which we hope is fundamental
    powerSpectrumTrunk = powerSpectrum(1:ceil(length(powerSpectrum)/2)); 
    powerSpectrumTrunk(1)=0;                            
    ampSpectrumTrunk = ampSpectrum(1:ceil(length(ampSpectrum)/2));
    whichFreq = find(powerSpectrumTrunk==max(powerSpectrumTrunk));  
    maxAmps(phase) = ampSpectrumTrunk(whichFreq(1));            
    sfs(phase) = (whichFreq(1)-1)./(length(croppedImage));
    theSfs{phase} = (0:(ceil(length(croppedImage)/2)-1))/length(croppedImage);
    ampSpectraTrunk{phase} = ampSpectrumTrunk;
    if (sfs(phase) ~= theSfs{phase}(whichFreq(1)))
        error('Oops on SF calcs');
    end
    
    % Plot
    if (PLOT && PLOTALL)
        figure(theFig); clf;  hold on
        plot(theSfs{phase},ampSpectraTrunk{phase},'r','LineWidth',1.5);
        plot([sfs(phase) sfs(phase)],[0 maxAmps(phase)],'k','LineWidth',2);
        pause(0.1);
    end
end

% Find max over croppings and report
whichPhase = find(maxAmps==max(maxAmps));                  
amp = maxAmps(whichPhase(1));                             
sf = sfs(whichPhase(1)); 

filteredFFT = zeros(size(theFFT{whichPhase}));
filteredFFT(1) = theFFT{whichPhase}(1);
index = find(abs(theFFT{whichPhase}(2:end)) == max(abs(theFFT{whichPhase}(2:end))));
fRange = 0;
f1 = (index(1)+1)-fRange:(index(1)+1)+fRange;
f1 = f1(f1 > 1 & f1 <= length(filteredFFT));
f2 = (index(2)+1)-fRange:(index(2)+1)+fRange;
f2 = f2(f2 > 1 & f2 <= length(filteredFFT));
filteredFFT(f1) = theFFT{whichPhase}(f1);
filteredFFT(f2) = theFFT{whichPhase}(f2);
filteredImage = abs(ifft(filteredFFT.*(length(theFFT{phase}(:)))));
amp = (max(filteredImage)-min(filteredImage))/(max(filteredImage)+min(filteredImage));

% Close plot if neessary
if (PLOT)
    figure(theFig); clf;
    subplot(1,2,1); hold on
    plot(theSfs{whichPhase },ampSpectraTrunk{whichPhase},'r','LineWidth',1.5);
    plot([sf sf],[0 amp],'k','LineWidth',2);
    xlim([0 0.5]); ylim([0 1]);
    subplot(1,2,2);
    plot(croppedImage); hold on
    plot(filteredImage,'r');
    xlim([0 length(croppedImage)]); ylim([0 2]);
    
    %pause;
    close(theFig);

end

end
