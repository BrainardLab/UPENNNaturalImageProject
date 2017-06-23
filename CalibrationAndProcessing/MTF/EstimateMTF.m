% EstimateMTF
%
% Process the grating images to estimate LMS MTF
%
% Parameters:
%   idbPath -- the path to the root of the calibration image database
%
% 12/22/10  dhb Merged in reading from raw images with older code.

function EstimateMTF(idbPath,PROCESSRAW)

% Close
close all;

%% Check optional arg
if (nargin < 2 || isempty(PROCESSRAW))
    PROCESSRAW = 0;
end

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

%% Get coords
coords = LoadMTFCoords;

%% Raw processing (a little slow)
if (PROCESSRAW == 1)
    % Loop over files and analyze
    if (CHECK)
        igFig = figure; clf;
        igFig2 = figure; clf;
        igFig3 = figure; clf;
    end
    for f = 1:length(coords)
        filenameReal = sprintf('DSC_%s',coords(f).name);
        if (isempty(coords(f).vert))
            coords(f).vert = 1;
        end
        fprintf('Processing file %d, filename %s, vert %d, plot %d\n',f,filenameReal,coords(f).vert,coords(f).plot);
        
        % Get exposure duration, etc
        cd([idbPath filesep coords(f).dir]);
        imageInfoReal = GetNEFInfo(filenameReal);
        
        % Get exposure duration, etc
        fprintf('\t\tCamera: %s\n',imageInfoReal.whichCamera);
        fprintf('\t\tExposure %g\n',imageInfoReal.exposure);
        fprintf('\t\tfStop %g\n',imageInfoReal.fStop);
        fprintf('\t\tISO %g\n',imageInfoReal.ISO);
        fprintf('\t\tFocus Mode %s\n',imageInfoReal.focusMode);
        
        % Some checks
        if (imageInfoReal.ISO ~= 1600)
            error('All data should be at ISO 1600');
        end
        if (imageInfoReal.fStop ~= 5.6)
            error('All data should be at f 5.6');
        end
        if (~strcmp(imageInfoReal.whichCamera,'standard'))
            error('This should be standard camera\n');
        end
        if (~strcmp(imageInfoReal.focusMode,'AF-S'))
            error('This should be autofocus AF-S\n');
        end
        
        % Read in image, dark subtract, and truncate negative values
        load([filenameReal '_LMS.mat']);
        load([filenameReal '_RGB.mat']);
        cd(curDir)
        
        % Check extracted location
        if (CHECK)
            ig = LMS_Image;
            ig = ig ./ max(ig(:));
            ig(ig<0) = 0;
            ig(coords(f).rows,coords(f).cols,:) = 1;
            figure(igFig); imagesc(ig);
        end
        
        LMS_Image = LMS_Image(coords(f).rows,coords(f).cols,:);
        RGB_Image = RGB_Image(coords(f).rows,coords(f).cols,:);
        CamCal = LoadCamCal(imageInfoReal.whichCamera);
        RGB_Image = DarkCorrect(CamCal,imageInfoReal,RGB_Image);
        
        % Rotate so that we match convention on which rest of code is written.
        if (coords(f).vert == 0)
            for k = 1:3
                LMS_Image1(:,:,k) = rot90(LMS_Image(:,:,k));
                RGB_Image1(:,:,k) = rot90(RGB_Image(:,:,k));
            end
            LMS_Image = LMS_Image1; clear LMS_Image1;
            RGB_Image = RGB_Image1; clear RGB_Image1;
        end
        
        % Extract image data
        nCols = size(LMS_Image,2);
        whichCols = 2:nCols-1;
        for k=1:3,
            lmsImage{f,k} = LMS_Image(:,whichCols,k);
            rgbImage{f,k} = RGB_Image(:,whichCols,k);
            lmsGrating{f,k} = sum(lmsImage{f,k},2);
            rgbGrating{f,k} = sum(rgbImage{f,k},2);
            plotType(f) = coords(f).plot;
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
        end
        
        % Get power and sf for every image
        [lsf(f),lamp(f)] = FindMaxAmp(lmsImage{f,1});
        [msf(f),mamp(f)] = FindMaxAmp(lmsImage{f,2});
        [ssf(f),samp(f)] = FindMaxAmp(lmsImage{f,3});
        [rsf(f),ramp(f)] = FindMaxAmp(rgbImage{f,1},1);
        [gsf(f),gamp(f)] = FindMaxAmp(rgbImage{f,2});
        [bsf(f),bamp(f)] = FindMaxAmp(rgbImage{f,3});
        
        drawnow;
        %pause;
        
    end
    save('ProcessedMTFData','lmsImage','rgbImage','lmsGrating','rgbGrating', 'plotType', ...
        'lsf','msf','ssf','rsf','gsf','bsf', ...
        'lamp','mamp','samp','ramp','gamp','bamp');
    
else
    load ProcessedMTFData
end

%% Get data for fitting
% The variable lets us control both what is plotted, what symbol is used, and
% whether it is included in the fit.
fitMeIndex = 1;
for f = 1:length(coords)   
    switch (coords(f).plot)
        case {1}
            % Accumulate data for fitting
            fitlsf(fitMeIndex) = lsf(f);
            fitmsf(fitMeIndex) = msf(f);
            fitssf(fitMeIndex) = ssf(f);
            fitrsf(fitMeIndex) = rsf(f);
            fitgsf(fitMeIndex) = gsf(f);
            fitbsf(fitMeIndex) = bsf(f);
            
            fitlamp(fitMeIndex) = lamp(f);
            fitmamp(fitMeIndex) = mamp(f);
            fitsamp(fitMeIndex) = samp(f);
            fitramp(fitMeIndex) = ramp(f);
            fitgamp(fitMeIndex) = gamp(f);
            fitbamp(fitMeIndex) = bamp(f);
            fitMeIndex = fitMeIndex + 1;
        otherwise
            
    end
end

%% Fit the MTF data
fitsf = linspace(0.0,0.45,200);
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
vlb = [0 -Inf -Inf 0];
vub = [1 Inf Inf 0.5];
x0 = [1 0.5 0.5 0];
sfs = fitlsf; amps = fitlamp;
truncate = 0;
lx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
lx = fmincon(@InlineMinFunction,lx1,[],[],[],[],vlb,vub,[],options);
[~,lfitMTF] = MinimizeMe(lx,fitsf,ones(size(fitsf)),truncate);

sfs = fitmsf; amps = fitmamp;
truncate = 0;
mx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
mx = fmincon(@InlineMinFunction,mx1,[],[],[],[],vlb,vub,[],options);
[~,mfitMTF] = MinimizeMe(mx,fitsf,ones(size(fitsf)),truncate);

sfs = fitssf; amps = fitsamp;
truncate = 0;
sx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
sx = fmincon(@InlineMinFunction,sx1,[],[],[],[],vlb,vub,[],options);
[~,sfitMTF] = MinimizeMe(sx,fitsf,ones(size(fitsf)),truncate);

sfs = fitrsf; amps = fitramp;
truncate = 0;
rx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
rx = fmincon(@InlineMinFunction,rx1,[],[],[],[],vlb,vub,[],options);
[~,rfitMTF] = MinimizeMe(rx,fitsf,ones(size(fitsf)),truncate);

sfs = fitgsf; amps = fitgamp;
truncate = 0;
gx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
gx = fmincon(@InlineMinFunction,gx1,[],[],[],[],vlb,vub,[],options);
[~,gfitMTF] = MinimizeMe(gx,fitsf,ones(size(fitsf)),truncate);

sfs = fitbsf; amps = fitbamp;
truncate = 0;
bx1 = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);
truncate = 1;
bx = fmincon(@InlineMinFunction,bx1,[],[],[],[],vlb,vub,[],options);
[~,bfitMTF] = MinimizeMe(bx,fitsf,ones(size(fitsf)),truncate);

%% Write out table of parameters
fid = fopen('MTFFitParams.txt','w');
fprintf(fid,'red: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',rx(1),rx(2),rx(3),rx(4));
fprintf(fid,'green: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',gx(1),gx(2),gx(3),gx(4));
fprintf(fid,'blue: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',bx(1),bx(2),bx(3),bx(4));
fprintf(fid,'long: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',lx(1),lx(2),lx(3),lx(4));
fprintf(fid,'medium: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',mx(1),mx(2),mx(3),mx(4));
fprintf(fid,'short: a = %0.2f, b = %0.2f, c = %0.2f, d = %0.2f\n',sx(1),sx(2),sx(3),sx(4));
fclose(fid);

%% Add fits and finalize figure format
mtfFigLMS = figure; clf; hold on;
mtfFigRGB = figure; clf; hold on;
figure(mtfFigLMS);
plot(fitsf,lfitMTF,'r','LineWidth',2);
plot(fitsf,mfitMTF,'g','LineWidth',2);
plot(fitsf,sfitMTF,'b','LineWidth',2);
figure(mtfFigRGB);
plot(fitsf,rfitMTF,'r','LineWidth',2);
plot(fitsf,gfitMTF,'g','LineWidth',2);
plot(fitsf,bfitMTF,'b','LineWidth',2);
for f = 1:length(coords)   
    switch (coords(f).plot)
        case {1}
            figure(mtfFigLMS);
            plot(lsf(f),lamp(f),'ro','markersize',8,'markerfacecolor','r');
            plot(msf(f),mamp(f),'go','markersize',8,'markerfacecolor','g');
            plot(ssf(f),samp(f),'bo','markersize',8,'markerfacecolor','b');
            
            figure(mtfFigRGB);
            plot(rsf(f),ramp(f),'ro','markersize',8,'markerfacecolor','r');
            plot(gsf(f),gamp(f),'go','markersize',8,'markerfacecolor','g');
            plot(bsf(f),bamp(f),'bo','markersize',8,'markerfacecolor','b');
        case {-1}
            figure(mtfFigLMS);
            plot(lsf(f),lamp(f),'ro','markersize',8);
            plot(msf(f),mamp(f),'go','markersize',8);
            plot(ssf(f),samp(f),'bo','markersize',8);
            
            figure(mtfFigRGB);
            plot(rsf(f),ramp(f),'ro','markersize',8);
            plot(gsf(f),gamp(f),'go','markersize',8);
            plot(bsf(f),bamp(f),'bo','markersize',8);
        otherwise
            
    end
end

figure(mtfFigLMS);
set(gca,'xlim',[0.0 0.5]);
set(gca,'ylim',[0 1.2]);
set(gca,'xtick',0:0.10:0.5);
xlabel('Spatial frequency (cycles/pixel)','fontsize',14);
set(gca,'fontsize',14);
ylabel('LMS Amplitude','fontsize',14);
axis square
FigureSave('MTFFigureLMS.pdf',mtfFigLMS,'pdf');

figure(mtfFigRGB);
set(gca,'xlim',[0.0 0.5]);
set(gca,'ylim',[0 1.2]);
set(gca,'xtick',0:0.10:0.5);
xlabel('Spatial frequency (cycles/pixel)','fontsize',14);
set(gca,'fontsize',14);
ylabel('RGB Amplitude','fontsize',14);
axis square
FigureSave('MTFFigureRGB.pdf',mtfFigRGB,'pdf');


%% Pat's grating data
OLDPATDATA = 0;
if (OLDPATDATA)
    load('GratingPatches','RedVals','GreenVals','BlueVals');
    if (CHECK)
        igFig4 = figure; clf;
    end
    for i = 1:length(BlueVals)
        if (CHECK)
            figure(igFig4); clear ig2
            ig2(:,:,1) = RedVals(i).Lum;
            ig2(:,:,2) = GreenVals(i).Lum;
            ig2(:,:,3) = BlueVals(i).Lum;
            ig2(ig2<0) = 0;
            subplot(1,3,1); imshow(ig2(:,:,1)/max(max(ig2(:,:,1))));
            subplot(1,3,2); imshow(ig2(:,:,2)/max(max(ig2(:,:,2))));
            subplot(1,3,3); imshow(ig2(:,:,3)/max(max(ig2(:,:,3))));
            drawnow;
            %pause
        end
        [prsf(i),pramp(i)] = FindMaxAmp(RedVals(i).Lum);
        [pgsf(i),pgamp(i)] = FindMaxAmp(GreenVals(i).Lum);
        [pbsf(i),pbamp(i)] = FindMaxAmp(BlueVals(i).Lum);
    end
    
    % Figure
    pmtfFigRGB = figure; clf; hold on;
    plot(prsf,pramp,'ro','markersize',8,'markerfacecolor','r');
    plot(pgsf,pgamp,'go','markersize',8,'markerfacecolor','g');
    plot(pbsf,pbamp,'bo','markersize',8,'markerfacecolor','b');
    
    set(gca,'xlim',[0.05 0.5]);
    set(gca,'ylim',[0 1.2]);
    set(gca,'xtick',0:0.10:0.5);
    xlabel('Spatial frequency (cycles/pixel)','fontsize',14);
    set(gca,'fontsize',14);
    ylabel('RGB Amplitude (Pat)','fontsize',14);
    axis square
    FigureSave('MTFFigureRGBPat.pdf',pmtfFigRGB,'pdf');
end

%% INLINE FUNCTION TO BE USED FOR MINIMIZATION.
% Inline functions have the feature that any variable they use that is
% not defined in the function has its value inherited
% from the workspace of wherever they were invoked.
    function [f,pred] = InlineMinFunction(x)
        [f,pred] = MinimizeMe(x,sfs,amps,truncate);
    end

end

function [f,pred] = MinimizeMe(x,sfs,amps,truncate)
pred = x(1)*exp(-x(2)*(sfs-x(4)).^2) + (1-x(1))*exp(-x(3)*sfs.^2);
if (truncate)
    index = sfs <= x(4);
    pred(index) = 1;
    pred(pred > 1) = 1;
end
resid = amps-pred;
f = sum(resid.^2);
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
nPhases = min([20 length(theImage)]);
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
