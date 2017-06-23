% DoTheCheck
%
% Go through a currently processed subset of the Botswana image database
% and compare what's there to what we downloaded from the database server.
% The point of this is to make sure that the way we are now processing the
% Philly images matches up with how the database was processed.
%
% To configure for running this, use:
%   tbUseProject({'UPENNNaturalImageProject'},'reset','full');

% 6/23/17  dhb  Started delving deeper.

%% Clear
clear; close all;

%% Cd into directory containing this routine
cd(fileparts(mfilename('fullpath')));

%% Parameters
percentDiscrepancyThreshold = 1;

%% Do checks on these directories
fromServerPath = getpref('UPENNNaturalImageProject','botswanaCheckFromServer');
checkMePath = getpref('UPENNNaturalImageProject','botswanaCheck');

%% Compare current calibration spectral sensitivity files with those archived from 2011
spectralSensitivityDir = fullfile(fileparts(mfilename('fullpath')),'..','CalibrationAndProcessing','SPECTRAL_SENSITIVITY');
mccDir = fullfile(fileparts(mfilename('fullpath')),'..','CalibrationAndProcessing','MCC');
calibrationDir = getpref('UPENNNaturalImageProject','calibrationDir');
the2011T_camera = load(fullfile(calibrationDir,'CameraSpectralSensFiles2011','T_camera2010'));
theCurrentT_camera = load(fullfile(spectralSensitivityDir,'T_camera2010'));
maxAbsDiff = max(abs(the2011T_camera.T_camera(:) - theCurrentT_camera.T_camera(:)));
if (maxAbsDiff < eps)
    fprintf('2011 and current T_camera files are the same to machine precision\n');
else
    meanAbsValue = mean(abs(the2011T_camera.T_camera(:)));
    maxPercentDiff = 100*maxAbsDiff/meanAbsValue;
    fprintf('2011 and current T_camera files differ at least a little\n');
    fprintf('\tMax abs difference: %0.8e, mean abs value: %0.2g, max percent diff = %0.4f\n',maxAbsDiff,meanAbsValue,maxPercentDiff);
end

%% Compare current LMS fundamentals to those archived from 2012
the2012T_cones = load(fullfile(calibrationDir,'PTBDataFiles2012','T_cones_ss2'));
theCurrentT_cones = load('T_cones_ss2');
maxAbsDiff = max(abs(the2012T_cones.T_cones_ss2(:) - theCurrentT_cones.T_cones_ss2(:)));
if (maxAbsDiff < eps)
    fprintf('2012 and current T_cones_ss2 files are the same to machine precision\n');
else
    meanAbsValue = mean(abs(the2012T_cones.M_RGBToLMS(:)));
    maxPercentDiff = 100*maxAbsDiff/meanAbsValue;
    fprintf('2012 and current T_cones_ss2 files differ at least a little\n');
    fprintf('\tMax abs difference: %0.8e, mean abs value: %0.2g, max percent diff = %0.4f\n',maxAbsDiff,meanAbsValue,maxPercentDiff);
end

%% Compare current LMS transform matrix to that archived from 2011
the2011M_RGBToLMS = load(fullfile(calibrationDir,'CameraSpectralSensFiles2011','M_RGBToLMS2010'));
theCurrentM_RGBToLMS = load(fullfile(spectralSensitivityDir,'M_RGBToLMS2010'));
maxAbsDiff = max(abs(the2011M_RGBToLMS.M_RGBToLMS(:) - theCurrentM_RGBToLMS.M_RGBToLMS(:)));
if (maxAbsDiff < eps)
    fprintf('2011 and current M_RGBToLMS files are the same to machine precision\n');
else
    meanAbsValue = mean(abs(the2011M_RGBToLMS.M_RGBToLMS(:)));
    maxPercentDiff = 100*maxAbsDiff/meanAbsValue;
    fprintf('2011 and current M_RGBToLMS files differ at least a little\n');
    fprintf('\tMax abs difference: %0.8e, mean abs value: %0.2g, max percent diff = %0.4f\n',maxAbsDiff,meanAbsValue,maxPercentDiff);
end

%% Compare current isomerizations to those from 2011
the2011LMSToIsomerizations= load(fullfile(calibrationDir,'IsomerizationsFiles2011','LMSToIsomerizations'));
theCurrentLMSToIsomerizations = load(fullfile(mccDir,'LMSToIsomerizations'));
maxAbsDiff = max(abs(the2011LMSToIsomerizations.LMSToIsomerizations(:) - theCurrentLMSToIsomerizations.LMSToIsomerizations(:)));
if (maxAbsDiff < eps)
    fprintf('2011 and current LMSToIsomerizations files are the same to machine precision\n');
else
    meanAbsValue = mean(abs(the2011LMSToIsomerizations.LMSToIsomerizations(:)));
    maxPercentDiff = 100*maxAbsDiff/meanAbsValue;
    fprintf('2011 and current LMSToIsomerizations files differ at least a little\n');
    fprintf('\tMax abs difference: %0.8e, mean abs value: %0.2g, max percent diff = %0.4f\n',maxAbsDiff,meanAbsValue,maxPercentDiff);
end

%% Derive the standard camera RGBToLMS transformation matrix for a pair of RGB and LMS files
theRGBFileFromServer = load(fullfile(fromServerPath,'out','cd17B_closeup_palm_nut_fresh_shade_sun','DSC_0208_RGB'));
theLMSFileFromServer = load(fullfile(fromServerPath,'out','cd17B_closeup_palm_nut_fresh_shade_sun','DSC_0208_LMS'));
theRGB = ImageToCalFormat(theRGBFileFromServer.RGB_Image);
theLMS = ImageToCalFormat(theLMSFileFromServer.LMS_Image);
M_RGBToLMSFromServer = ((theRGB')\(theLMS'))';
M_RGBToLMSFromServer./the2011M_RGBToLMS.M_RGBToLMS

%Cam_Cal = LoadCamCal(theImage.imageInfo.whichCamera);

% load T_cones_ss2
% T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S_camera);
% load T_ss2000_Y2 
% T_Y = SplineCmf(S_ss2000_Y2,683*T_ss2000_Y2,S_camera);
% M_RGBToLMS = ((T_camera')\(T_cones'))';
% T_cones_check = M_RGBToLMS*(T_camera);
% M_LMSToLum = ((T_cones')\(T_Y'))';

%% Say hello
disp('**************** BEGIN ****************')

%% Get directories to check
theDirectoriesFromServer = dir(fullfile(fromServerPath,'out','*'));
theDirectoriesToCheck = dir(fullfile(checkMePath,'out','*'));

%% Loop over each directory and do the checks for what's inside
indexRGB = 1;
indexLUM = 1;
indexLMS = 1;
for ii = 1:length(theDirectoriesFromServer)
    totalDiscrepancy = 0;
    
    % Check that directories match, and exclude annoying dirs that start
    % with .
    if (theDirectoriesFromServer(ii).name(1) ~= '.')
        if (~strcmp(theDirectoriesFromServer(ii).name,theDirectoriesToCheck(ii).name))
            error('Directories we are supposed to be checking do not match up');
        end
        
        % Announce which directory we are checking.
        fprintf('***** Checking files in directory %s\n',theDirectoriesFromServer(ii).name)
        
        % Get the list of files in the directory
        fromServerDirectoryPath = fullfile(fromServerPath, 'out', theDirectoriesFromServer(ii).name);
        checkMeDirectoryPath = fullfile(checkMePath, 'out', theDirectoriesToCheck(ii).name);
        fromServerDirectory = dir(fromServerDirectoryPath);
        checkMeDirectory = dir(checkMeDirectoryPath);
        
        % Create cells of file names to iterate through. Ignore anything
        % that is not a file.
        fromServerIsDir = [fromServerDirectory.isdir];
        fileNamesFromServer = {fromServerDirectory(~fromServerIsDir).name};
        checkMeIsDir = [checkMeDirectory.isdir];
        fileNamesCheckMe = {checkMeDirectory(~checkMeIsDir).name};
        
        % Compare corresponding LMS, LUM, and RGB files. Note we don't care
        % about the AUX files.
        for jj = 1:size(fileNamesCheckMe, 2)
            
            % Grab the extension.
            [pathstr, name, ext] = ...
                fileparts(fullfile(fromServerDirectoryPath, fileNamesCheckMe{jj}));
            notAUX = isempty(strfind(fileNamesCheckMe{jj}, 'AUX'));
            if (length(name) > 3)
                type = name(end-2:end);
            else
                type = 'IGNORE';
            end
            
            % Check for corresponsing file in theDirectoriesFromServer.  If
            % not, we skip this comparison.
            if (strcmp(ext, '.mat') && ~strcmp(fileNamesCheckMe{jj}(1:2), '._') && notAUX && ...
                    any(strcmp(fileNamesFromServer, fileNamesCheckMe{jj})))
                
                % Read the data from both files
                checkMeData = load(fullfile(checkMeDirectoryPath, fileNamesCheckMe{jj}));
                fromServerData = load(fullfile(fromServerDirectoryPath, fileNamesCheckMe{jj}));
                
                % There are three file types we'd like to compare, and we
                % accumulate each separately.
                switch (type)
                    case 'RGB'
                        checkImage = checkMeData.RGB_Image;
                        fromServerImage = fromServerData.RGB_Image;
                    case 'LUM'
                        checkImage = checkMeData.LUM_Image;
                        fromServerImage = fromServerData.LUM_Image;
                    case 'LMS'
                        checkImage = checkMeData.LMS_Image;
                        fromServerImage = fromServerData.LMS_Image;
                    otherwise
                end
                
                % Find the absolute difference
                maxImageDiscrepancy = max(abs(checkImage(:) - fromServerImage(:)));
                meanImageVal = mean(checkImage(:));
                percentImageDiscrepancy = 100*maxImageDiscrepancy/meanImageVal;
                
                % Print out names of files who do not meet the precision
                % threshold. Note I chose 7 x 10^7 abritrarily based on the
                % histograms.
                if (percentImageDiscrepancy  > percentDiscrepancyThreshold)
                    fprintf('File %s discrepancy exceeds threshold: discrepancy: %0.2g, mean val %0.2g, percent discrepancy %0.1f\n',...
                        fileNamesCheckMe{jj},maxImageDiscrepancy,meanImageVal,percentImageDiscrepancy);
                end
                
                % Accumulate information, keeping image type straight.
                switch (type)
                    case 'RGB'
                        maxImageDiscrepancy_RGB(indexRGB) = maxImageDiscrepancy;
                        meanImageVal_RGB(indexRGB) = meanImageVal;
                        percentImageDiscrepancy_RGB(indexRGB) = percentImageDiscrepancy;
                        indexRGB = indexRGB+1;
                    case 'LUM'
                        maxImageDiscrepancy_LUM(indexLUM) = maxImageDiscrepancy;
                        meanImageVal_LUM(indexLUM) = meanImageVal;
                        percentImageDiscrepancy_LUM(indexLUM) = percentImageDiscrepancy;
                        indexLUM = indexLUM+1;
                    case 'LMS'
                        maxImageDiscrepancy_LMS(indexLMS) = maxImageDiscrepancy;
                        meanImageVal_LMS(indexLMS) = meanImageVal;
                        percentImageDiscrepancy_LMS(indexLMS) = percentImageDiscrepancy;
                        indexLMS = indexLMS+1;
                    otherwise
                end       
            end
        end
    end
end

%% Summarize what hapened for each image type
if (max(abs(percentImageDiscrepancy_RGB) ~= 0))
    RGBfig = figure;
    histogram(percentImageDiscrepancy_RGB);
    title('RGB Percent Discrepancies');
    xlabel('Percent Discrepancy');
    ylabel('Count');
end
fprintf('Maximum absolute RGB percentage discrepancy is %0.2f\n',max(percentImageDiscrepancy_RGB(:)));

if (max(abs(percentImageDiscrepancy_LUM) ~= 0))
    LUMfig = figure;
    histogram(percentImageDiscrepancy_LUM);
    title('LUM Percent Discrepancies');
    xlabel('Percent Discrepancy');
    ylabel('Count');
end
fprintf('Maximum absolute LUM percentage discrepancy is %0.2f\n',max(percentImageDiscrepancy_LUM(:)));

if (max(abs(percentImageDiscrepancy_LMS) ~= 0))
    LMSfig = figure;
    histogram(percentImageDiscrepancy_LMS);
    title('LMS Percent Discrepancies');
    xlabel('Percent Discrepancy');
    ylabel('Count');
end
fprintf('Maximum absolute LMS percentage discrepancy is %0.2f\n',max(percentImageDiscrepancy_LMS(:)));

disp('***************** END *****************')

