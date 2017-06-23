% EstimateLumCheck
%
% Make sure estimated luminance matches directly measured for a bunch of
% images of the same thing, with different shutter speed/aperture combos.
%
% Measured luminance was 368 cd/m2 for LUM368 images, and 320 c/m2 for
% LUM320 images.  368 was under the room fluorescents, while 320 was under
% the slide projector tungsten.
%
% Parameters:
%   idbPath -- the path to the root of the calibration image database
%
% 12/23/10  dhb Wrote.

function EstimateLumCheck(idbPath)

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

%% Recreate the luminance check plot from the MCC data
load ../MCC/MCCLuminanceCheckData
lumPlot = figure; clf;
loglog(luminanceFromRadiometer ,luminanceFromCamera ,'ko','MarkerFaceColor','k'); hold on
loglog([1e0 1e3],[1e0 1e3],'k','linewidth',1.3);
axis square;
set(gca,'fontsize',14);
xlabel('Luminance from radiometer','fontsize',14);
ylabel('Luminance from camera','fontsize',14);
xlim([1e0 1e3]); ylim([1e0 1e3]);

%% Loop over our choices
for which = [320 368]
    
    % These are the coordinates of the MCC check box being viewed
    % nQ x nQ is the size of the region in pixels
    nQ = 50;
    switch (which)
        case 320
            crop_coords = [460 755];
            cd([idbPath '/LUMCHECK/LUM320']);
        case 368
            crop_coords = [455 775];
            cd([idbPath '/LUMCHECK/LUM368']);
        otherwise
            error('Bad choice of which to process');
    end
    
    % Load the photos.
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
        
        if (~strcmp(imageInfoReal.whichCamera,'standard'))
            error('This should be standard camera\n');
        end
        
        % Check response range
        load([filenameReal '_RGB.mat']);
        RGB_Image = RGB_Image(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ,:);
        CamCal = LoadCamCal(imageInfoReal.whichCamera);
        RGB_Image = DarkCorrect(CamCal,imageInfoReal,RGB_Image);
        if (any(RGB_Image < fitLowResp) | any(RGB_Image > fitHighResp))
            fprintf('\tWARNING: Some data out of good response range\n');
        end
        
        % Read in raw image, dark subtract, and truncate negative values
        load([filenameReal '_Lum.mat']);
        if (CHECK)
            ig = LUM_Image;
            ig = ig ./ max(ig(:));
            ig(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ) = 1;
            warning('off','Images:initSize:adjustingMag');
            figure(igFig);clf;imshow(ig);
            warning('on','Images:initSize:adjustingMag');
            if (f == 3)
                cd(curDir);
                switch (which)
                    case 320
                        imwrite(ig,'PhotoRegions320.jpg','jpg');
                        cd([idbPath '/LUMCHECK/LUM320']);
                    case 368
                        imwrite(ig,'PhotoRegions368.jpg','jpg');
                        cd([idbPath '/LUMCHECK/LUM368']);
                    otherwise
                        error('Bad choice of which to process');
                end
            end
            %pause
        end
        
        lumFromFile(f) = mean(mean(double(LUM_Image(crop_coords(1):crop_coords(1)+nQ, crop_coords(2):crop_coords(2)+nQ))));
        figure(lumPlot);
        switch (which)
            case 320
                fprintf('Got luminance of %g cd/m2, expect 320 cd/m2\n\n',lumFromFile(f));
                plot(320,lumFromFile(f),'ro','MarkerFaceColor','r');
            case 368
                fprintf('Got luminance of %g cd/m2, expect 368 cd/m2\n\n',lumFromFile(f));
                plot(368,lumFromFile(f),'go','MarkerFaceColor','g');
            otherwise
                error('Bad choice of which to process');
        end
    end
    cd(curDir)   
end

%% Save the figure
FigureSave('LuminanceCheckPlot.pdf',lumPlot,'pdf');

end


