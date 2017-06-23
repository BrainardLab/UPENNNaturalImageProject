% MakeCamCal
%
% Put together the pieces to make the camera calibration structure
%
% 01/25/10  dhb  Wrote this.
% 08/12/10  gt   Changed into a function.
% 08/20/10  gt   Computes also the auxilliary camera calibration.
% 11/10/10  dhb  Don't return calibration structure -- that could cause confusion because there is no specification of which one.
%                Use LoadCamCal to get the structure.
%           dhb  Move the created files into the D70Toolbox automatically.
% 12/14/10  dhb  Change camera scale factor to 0.80, as per most recent analysis.
% 12/20/10  dhb  Store M_RGBToLMS computed when spectral sensitivities created.  Scale for auxilliary camera.
% 12/21/10  dhb  Use 2010 spectral sensitivities and RGBToLMS.  And, DIVIDE standard by 0.80 to get aux M_RGBToLMS.
%                (Multiplying by 0.80 is still right for spectral sensitivities, I am pretty sure.)

%% Clear out old stuff on general principles;
function MakeCamCal()
    curDir = pwd;

    %% Standard camera
    Cam_Cal.serialNumber = '2000a9a7';
    Cam_Cal.whichCamera = 'standard';

    % Spectral sensitivity and RGB to LMS, etc
    cd SPECTRAL_SENSITIVITY
    load T_camera2010
    Cam_Cal.T_camera = T_camera;
    Cam_Cal.S_camera = S_camera;
    load M_RGBToLMS2010
    Cam_Cal.M_RGBToLMS = M_RGBToLMS;
    Cam_Cal.M_LMSToLum = M_LMSToLum;
    cd(curDir);

    % Dark current
    cd DARK_IMAGES
    load darkTable
    Cam_Cal.darkTable = darkTable;
    cd(curDir);

    % Isomerizations conversion
    cd MCC
    load LMSToIsomerizations
    Cam_Cal.LMSToIsomerizations = LMSToIsomerizations;
    cd(curDir);

    % Save it
    save StandardD70Data Cam_Cal
    disp('Saved standard camera calibration...');
    
    %% Create auxilliarry calibration assuming that the response of the
    % auxilliary camera is 0.80 * (standard camera), given the same
    % spectra.  Other data are the same.
    Cam_Cal.serialNumber = '20004b72';
    Cam_Cal.whichCamera  = 'auxilliary';
    Cam_Cal.T_camera     = Cam_Cal.T_camera * 0.80;
    Cam_Cal.M_RGBToLMS   = Cam_Cal.M_RGBToLMS / 0.80;
    save AuxiliaryD70Data Cam_Cal
    disp('Saved auxilliary camera calibration...');
    
    %% Move files to D70Toolbox
    D70Dir = GetToolboxDirectory('D70Toolbox');
    CamCalDir = [D70Dir filesep 'CameraData'];
    unix(['mv StandardD70Data.mat ' CamCalDir]);
    unix(['mv AuxiliaryD70Data.mat ' CamCalDir]);
end


function toolboxDir = GetToolboxDirectory(toolboxName, suppressWarning)
% toolboxDir = GetToolboxDirectory(toolboxName, [suppressWarning])
%
% Description:
% Returns the root directory of the requested toolbox.
%
% Input:
% toolboxName (string) - Name of the toolbox.
%
% Optional Input:
% suppressWarning (logical) - If true, this function doesn't print out a
%    warning message if the toolbox wasn't found.  Defaults to false.
%
% Output:
% toolboxDir (string) - Root directory of the toolbox or empty if the
%	toolbox isn't found.
%

if nargin < 1 || nargin > 2
	error('toolboxDir = GetToolboxDirectory(toolboxName, [suppressWarning])');
end

if nargin == 1
	suppressWarning = false;
end

p = path;
toolboxDir = [];

% Parse the path.
x = textscan(p, '%s', 'Delimiter', ':');
x = x{1};

% Look at each entry to find the toolbox.
for i = 1:length(x)
	[p, f] = fileparts(x{i});
	
	if ~isempty(strmatch(f, toolboxName, 'exact'))
		toolboxDir = sprintf('%s/%s', p, f);
	end
end

if isempty(toolboxDir) && ~suppressWarning
	fprintf('*** Toolbox "%s" not found on the path.', toolboxName);
end

end



