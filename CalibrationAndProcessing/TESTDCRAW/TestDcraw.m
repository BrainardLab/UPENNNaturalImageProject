% TestDcraw
%
% ****** DHB, 6/23/17
% NOTE: This does not currently run.  We seem to be missing an
% executable called dcraw, which the script expects to exist.  Having this
% run right now is not a priority, so I'm not digging into it.  I did move
% the image data out and generate a preference to point to it.  Since the
% script doesn't run, that is not fully tested.
% ******* 
%
% Compare orig .ppm and new .pgm dcraw outputs,
% as well as compare effect of -D versus -d
% options.
%
% Requires SimToolbox (although it would be pretty easy to eliminate this dependency.)
%
% In the table below, "orig" refers to whatever version of dcraw we used in 2004 when we did the
% calibration.  We only have the .ppm files, so we proceed by comparing current output with what
% is in those.
%
% For all comparisons, the values obtained in different ways are related by a simple change of gain,
% although which pixels saturate can vary with the version.
%
% All comparisons were tried with 2 different NEF files, with same conclusion reached for each.  So,
% probably, the differences between versions are independent of image content.
%
% Version 8.99 (rcs 1.432) -- With -d option, max is ~65535 and New->Orig RGB gains 0.26, 0.25, 0.22
%                          -- With -D option, max is ~65535 and New->Orig RGB gains 11.20, 4.00, 4.80
% Version 7.18 (from Pat)  -- Behaves like current (v8.99) for -d except writes 3 plane .ppm rather than 1 plane .pgm
% Version 5.88 (rcs 1.200) -- Like orig version, but max is ~65535 and New->Orig RGB gains 0.25, 0.25, 0.25
% Version 5.71 (rcs 1.180) -- Like orig version
% Version 5.42 (rcs 1.160) -- Max is ~16384 like orig, but New->Orig RGB gains 2.80, 1.00, 1.2
%                          -- Note that the ratio of these gains is the same as for the 8.99 -D option,
%                          -- just scaled down by a factor of 4.

% 1/26/10  dhb  Wrote it.
% 1/27/10  dhb  Much more general, and document comparisons in the header.
% 6/223/17 dhb  Seprate out data and point to it with a pref.

%% Clear and close
clear; close all;

%% Parameters
whichNew = 'OldVer';                   % 'OldVer', 'LittleD', 'BigD'
oldRevision = '571';                   % See options in switch below
imageDir = getpref('UPENNNaturalImageProject','dcrawTestImages');
imageRoot = 'Gasper_DSC_0001';         % Root of any NEF file in the current directory
imagePath = fullfile(imageDir,imageRoot);

%% Create .ppm and .pgm files using new and old
% versions of dcraw
switch oldRevision
    case '718'
        unix(['./dcraw.v718 -d -4 ' imagePath '.NEF']);
    case '588'
        unix(['./dcraw.v588 -d -4 ' imagePath '.NEF']);
    case '571'
        unix(['./dcraw.v571 -d -4 ' imagePath '.NEF']);
    case '542'
        unix(['./dcraw.v542 -d -4 ' imagePath '.NEF']);
    otherwise
        error('Unknown dcraw version requested\n');
end
unix(['mv ' imagePath '.ppm ' imagePath sprintf('.v%s.ppm',oldRevision)]);

unix(['dcraw -d -4 ' imagePath '.NEF']);
unix(['mv ' imagePath '.pgm ' imagePath '.littleD.pgm']);

unix(['dcraw -D -4 ' imagePath '.NEF']);
unix(['mv ' imagePath '.pgm ' imagePath '.bigD.pgm']);

%% Read in output of dcraw we just produced
bigD = double(imread([imagePath '.bigD.pgm']));
littleD = double(imread([imagePath '.littleD.pgm']));
oldverThreePlane = double(imread([imagePath sprintf('.v%s.ppm',oldRevision)]));

%% Check the oldver version for ppm plane identity and then extract one of
% the three identical plane.s
oldver1 = oldverThreePlane(:,:,1);
oldver2 = oldverThreePlane(:,:,2);
oldver3 = oldverThreePlane(:,:,3);
oldver12Diff = max(abs(oldver1(:)-oldver2(:)));
oldver13Diff = max(abs(oldver1(:)-oldver3(:)));
fprintf('Old dcraw version max difference in planes is (1v2) %g, (1v3) %g\n',oldver12Diff,oldver13Diff);
if (oldver12Diff ~= 0 || oldver13Diff ~= 0)
    error('Unexpected difference between planes of old dcraw version ppm file\n');
end

%% Check big versus little D options.  They are in fact different.
bigLittleDiff = max(abs(bigD(:)-littleD(:)));
fprintf('New version -D and -d max difference is %g\n',bigLittleDiff);

%% Set which recently created "new" version to compare with the output from 2004.
switch (whichNew)
    case 'BigD';
        new = bigD;
        figName = sprintf('NewVOrig_%s_%s.png',imagePath,whichNew);
    case ('LittleD');
        new = littleD;  
        figName = sprintf('NewVOrig_%s_%s.png',imagePath,whichNew);
    case ('OldVer')
        new = oldver1(:,1:3039);
        figName = sprintf('NewVOrig_%s_v%s.png',imagePath,oldRevision);
    otherwise
        error('Unknown whichNew option specified\n');
end
clear bigD littleD oldverThreePlane oldver1 oldver2 oldver3

%% Get the original version (written in 2004) and check that
% its three planes are identical.
% Check orig planes 1 versus 2 and 3
origThreePlane = double(imread([imagePath '.orig.ppm']));
orig1 = origThreePlane(:,:,1);
orig2 = origThreePlane(:,:,2);
orig3 = origThreePlane(:,:,3);
orig12Diff = max(abs(orig1(:)-orig2(:)));
orig13Diff = max(abs(orig1(:)-orig3(:)));
fprintf('Orig version max difference in planes is (1v2) %g, (1v3) %g\n',orig12Diff,orig13Diff);
if (orig12Diff ~= 0 || orig13Diff ~= 0)
    error('Unexpected difference between planes of original version ppm file\n');
end

%% Set "orig" version
orig = double(orig1);
clear origThreePlane orig1 orig2 orig3

%% Pull out a region of orig versus new and have a look
% Can make this smaller to speed it up.
% 
% Currently just do the whole image, except need to delete
% last column to make old and new dimensions match
origRegion = orig(:,1:3039);
newRegion = new;

%% This code pops up figures of the image, if you want to stare at the pixels
%figure; imshow(origRegion/max([origRegion(:) ; newRegion(:)]));
%figure; imshow(newRegion/max([origRegion(:) ; newRegion(:)]));

%% Plot each color plane separately.  Requires SimToolbox to generate
% the mask.
simImage = MakeNikonSimImage(origRegion);
mosaicMask = SimCreateMask(simImage.cameraFile,simImage.height,simImage.width);
figure; clf; hold on
index = find(mosaicMask(:,:,1) == 1 & newRegion(:,:) < 0.9*max(max(newRegion(:,:))) &  origRegion(:,:) < 0.9*max(max(origRegion(:,:))));
plot(origRegion(index),newRegion(index),'r+');
redGain = (newRegion(index)\origRegion(index));
index = find(mosaicMask(:,:,2) == 1 & newRegion(:,:) < 0.9*max(max(newRegion(:,:))) & origRegion(:,:) < 0.9*max(max(origRegion(:,:))));
plot(origRegion(index),newRegion(index),'g+');
greenGain = (newRegion(index)\origRegion(index));
index = find(mosaicMask(:,:,3) == 1 & newRegion(:,:) < 0.9*max(max(newRegion(:,:))) &  origRegion(:,:) < 0.9*max(max(origRegion(:,:))));
plot(origRegion(index),newRegion(index),'b+');
blueGain = (newRegion(index)\origRegion(index));
axis([0 max([origRegion(:) ; newRegion(:)]) 0 max([origRegion(:) ; newRegion(:)])]);
axis('square');
xlabel('2004 dcraw pixel value');
ylabel('2009 dcraw pixel value');
title(sprintf('RGB gains %0.2f, %0.2f, %0.2f',redGain,greenGain,blueGain));
saveas(gcf,figName,'png');
fprintf('New->Orig RGB gains %0.2f, %0.2f, %0.2f\n',redGain,greenGain,blueGain);