% DoTheCheck
%
% Go through a currently processed subset of the Botswana image database
% and compare what's there to what we downloaded from the database server.
% The point of this is to make sure that the way we are now processing the
% Philly images matches up with how the database was processed.
%
% To configure for running this, use:
%   tbUseProject({'UPENNNaturalImageProject'},'reset','full');

%% Clear
clear; close all;

%% Parameters
percentDiscrepancyThreshold = 1;

%% Process the following directories, each directory represents a set of images
fromServerPath = getpref('UPENNNaturalImageProject','botswanaCheckFromServer');
checkMePath = getpref('UPENNNaturalImageProject','botswanaCheck');

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

