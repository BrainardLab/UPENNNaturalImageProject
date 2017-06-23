function UPENNNaturalImageProjectLocalHook
% UPENNNaturalImageProjectLocalHook - Configure for the project
% 
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "Template" from the filename,
% this will get run when you execute tbUseProject('UPENNNaturalImageProject') to set up for
% this project.  You then edit your local copy to match your local machine.
% The call to tbUseProject will do the copy, if there is no local hook
% already present on your machine.
%
% You will need to edit the project prefs to point at the data on your
% computer.

% 6/19/17  dhb, as  Update template for Annie's computer but didn't test.

%% Say hello
fprintf('Running UPENNNaturalImageProject hook\n');

%% Clear prefs
if (ispref('UPENNNaturalImageProject'))
    rmpref('UPENNNaturalImageProject');
end

%% Setup basedir with good guesses
sysInfo = GetComputerInfo();
switch (sysInfo.localHostName)
    case 'eagleray'
        % DHB's desktop
        baseDir = fullfile(filesep,'Volumes','Users1','DropboxLab/UPENNNaturalImageProject');
    case 'Annies-MacBook-Pro-2'
        % Annie's laptop
        baseDir = fullfile(filesep,'/Volumes','Annie','UPENNNaturalImageProject');
    otherwise
        % Some unspecified machine, try user specific customization
        switch(sysInfo.userShortName)
            % Could put user specific things in, but at the moment generic
            % is good enough.
            otherwise
                baseDir = ['/Users/' sysInfo.userShortName 'DropboxLab//UPENNNaturalImageProject'];
        end
end

%% Set preferences

% Images used in the calibration script
setpref('UPENNNaturalImageProject','calibrationImageDir',fullfile(baseDir,'Calibration','Images','calibration'));

% Botswana database
setpref('UPENNNaturalImageProject','botswanaDatabase',fullfile(baseDir,'BotswanaImagesOrig'));

% Philly database
setpref('UPENNNaturalImageProject','phillyDatabase',fullfile(baseDir,'PhillyImages'));

% Botswana check images
setpref('UPENNNaturalImageProject','botswanaCheck',fullfile(baseDir,'BotswanaImagesCheck_061617'));

% Botswana images downloaded from the server
setpref('UPENNNaturalImageProject','botswanaCheckFromServer',fullfile(baseDir,'BotswanaImagesFromServer_061617'));

% Test images for Dcraw tests.
setpref('UPENNNaturalImageProject','dcrawTestImages',fullfile(baseDir,'DcrawTestImages'));


