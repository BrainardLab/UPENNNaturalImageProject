% CheckExposureLinearity
% 
% Read in camera responses as a function of exposure duration and
% check linearity.
%
% 6/9/05    cpr         Wrote it.
% 6/22/05   cpr, dhb    Next generation.
% 8/05/05   cpr, dhb
% 8/08/05   cpr, pbg
% 8/11/05   cpr, pbg

% Clear out old variables
clear all;

% Load the raw data.  This has not been corrected for camera dark response.
% The file contains 168 images.  Numbers 56, 112, and 168 correspond to
% the "Bulb" setting on the camera and are not meaningful.  
load('/Volumes/ColorShare1/CalibrateD70/EXPOSURE_LINEARITY/Data/RGBInfoRaw.txt','-ascii');

% Load it in another form.  This is the f11 data in a more complete form.
%
% In this version, column 1 is exposure time, columns 2:4 are rgb values, and column 5 is the f-number (11 for all)
load('/Volumes/ColorShare1/CalibrateD70/EXPOSURE_LINEARITY/Data/RGBInfo_f11_full.txt','-ascii')

% Check for agreement.  This looks good.
exposureTimes = RGBInfo_f11_full(:,1);
RGBValues_f11_raw = RGBInfoRaw(1:55,:);
RGBValues_f11_full = RGBInfo_f11_full(:,2:4);
figure; clf; hold on
plot(RGBValues_f11_raw(:,1),RGBValues_f11_full(:,1),'r+');
plot(RGBValues_f11_raw(:,2),RGBValues_f11_full(:,2),'g+');
plot(RGBValues_f11_raw(:,3),RGBValues_f11_full(:,3),'b+');
plot([min(RGBValues_f11_raw(:)) max(RGBValues_f11_raw(:))],[min(RGBValues_f11_raw(:)) max(RGBValues_f11_raw(:))],'k');
title('Check Data Consistency');

% Given that RGBValues_f11_raw and RGBValues_f11_full are in agreement, we
% can confidently extract the f5 and f16 data from the raw file.  The
% exposure times are consistent across f numbers, so we use the same
% exposureTimes variable for all three.
RGBValues_f16_raw = RGBInfoRaw(57:111,:);
RGBValues_f5_raw = RGBInfoRaw(113:167,:);

% Read in the dark response matrix. This is 55x4, where the first column is exposure
% time, and the other 3 are rgb dark values.
darkRGBTable = load('DarkCurrentPhillyTable.txt','-ascii');

% Dark correct the RGB values.
for i = 1:length(exposureTimes)
    RGBDark = GetDarkRGB2(exposureTimes(i),darkRGBTable);
    RGBValues_f11(i,:) = RGBValues_f11_raw(i,:)-RGBDark;
    RGBValues_f16(i,:) = RGBValues_f16_raw(i,:)-RGBDark;
    RGBValues_f5(i,:) = RGBValues_f5_raw(i,:)-RGBDark;
end

% Make some plots for each f value
index = find(RGBValues_f11 > 0);
logRGBValues_f11 = NaN*ones(size(RGBValues_f11));
logRGBValues_f11(index) = log10(RGBValues_f11(index));
figure; clf; hold on
plot(log10(exposureTimes),logRGBValues_f11(:,1),'r+')
plot(log10(exposureTimes),logRGBValues_f11(:,2),'g+')
plot(log10(exposureTimes),logRGBValues_f11(:,3),'b+')
title('rgb vs. exp time, f 11');

index = find(RGBValues_f16 > 0);
logRGBValues_f16 = NaN*ones(size(RGBValues_f16));
logRGBValues_f16(index) = log10(RGBValues_f16(index));
figure; clf; hold on
plot(log10(exposureTimes),logRGBValues_f16(:,1),'r+')
plot(log10(exposureTimes),logRGBValues_f16(:,2),'g+')
plot(log10(exposureTimes),logRGBValues_f16(:,3),'b+')
title('rgb vs. exp time, f 16');

index = find(RGBValues_f5 > 0);
logRGBValues_f5 = NaN*ones(size(RGBValues_f5));
logRGBValues_f5(index) = log10(RGBValues_f5(index));
figure; clf; hold on
plot(log10(exposureTimes),logRGBValues_f5(:,1),'r+')
plot(log10(exposureTimes),logRGBValues_f5(:,2),'g+')
plot(log10(exposureTimes),logRGBValues_f5(:,3),'b+')
title('rgb vs. exp time, f 5');

% Now we want to see what part of each graph is linear, how linear it is,
% and throughout what range linearity holds.

% By inspection, choose upper and lower bounds for the linear fit: current
% values chosen by cpr on 6-23-05
logExpMax_f11=-0.4;
logExpMin_f11=-3;

logExpMax_f5=-1.2;
logExpMin_f5=-3.5;

logExpMax_f16=-0.1;
logExpMin_f16=-2.5;

% Find the exact exposure times and indices that lie between the guesses:
UsableExpTimes=zeros(1,6);

fNumbers=[5 11 16];
for i=1:length(fNumbers);
    fNumber=fNumbers(i);
    eval(['uu=log10(exposureTimes)-logExpMin_f' num2str(fNumber) '*ones(size(exposureTimes));'])
    eval(['vv=log10(exposureTimes)-logExpMax_f' num2str(fNumber) '*ones(size(exposureTimes));'])
    a=find(abs(uu)==min(abs(uu)));
    b=find(abs(vv)==min(abs(vv)));
    
    if length(a) ~= 1
        disp('error in a at i=')
        i
    end
 
    if length(b) ~= 1
        disp('error in b at i=')
        i
    end
    
    eval(['logUsableRGBValues_f' num2str(fNumber) '=logRGBValues_f' num2str(fNumber) '(a:b,:);'])
    eval(['logUsableExpTimes_f' num2str(fNumber) '=log10(exposureTimes(a:b));'])
    
    ii=1+2*(i-1);
    UsableExpTimes(ii)=a;
    UsableExpTimes(ii+1)=b;
    
end

% Find the best fitting slope to the lines, report it, and plot it against
% the data.

for i=1:length(fNumbers)
    fNumber=fNumbers(i);
    eval(['X=logUsableExpTimes_f' num2str(fNumber) ';'])
    eval(['Y_r=logUsableRGBValues_f' num2str(fNumber) '(:,1);'])
    eval(['Y_g=logUsableRGBValues_f' num2str(fNumber) '(:,2);'])
    eval(['Y_b=logUsableRGBValues_f' num2str(fNumber) '(:,3);'])
    P_r=polyfit(X,Y_r,1);
    P_g=polyfit(X,Y_g,1);
    P_b=polyfit(X,Y_b,1);
    if i==1
        yInt_f5=[P_r(2) P_g(2) P_b(2)];
    end
    eval(['Slope_f' num2str(fNumber) '_r=P_r(1);'])
    eval(['Slope_f' num2str(fNumber) '_g=P_g(1);'])
    eval(['Slope_f' num2str(fNumber) '_b=P_b(1);'])
    PP_r=polyval(P_r,X);
    PP_g=polyval(P_g,X);
    PP_b=polyval(P_b,X);
    figure; hold on
    eval(['plot(log10(exposureTimes),logRGBValues_f' num2str(fNumber) '(:,1),''r+'')'])
    eval(['plot(log10(exposureTimes),logRGBValues_f' num2str(fNumber) '(:,2),''g+'')'])
    eval(['plot(log10(exposureTimes),logRGBValues_f' num2str(fNumber) '(:,3),''b+'')'])
    plot(X,PP_r,'r')
    plot(X,PP_g,'g')
    plot(X,PP_b,'b')
    hold off
    eval(['ww=[''rgb vs. exp time, f'' num2str(fNumber) ] ;'])
    title(ww)
end

for i=1:length(fNumbers)
    fNumber=fNumbers(i);
    eval(['Slope_f' num2str(fNumber) '_r'])
    eval(['Slope_f' num2str(fNumber) '_g'])
    eval(['Slope_f' num2str(fNumber) '_b'])
end

% from the slopes of log pixel value vs log exp time, get pixel value as
% function of exposure time for the f5 aperture, and then set other
% aperture values to this one.
% the pixel value y = exp(slope*log10(exp time))

% set bounds for linearity
% exposureTimes(7)=1/2000 sec.
% exposureTimes(40)=1 sec.

MinExpTime=exposureTimes(7);
MaxExpTime=exposureTimes(40);

% calculate the average slope for each channel:
% 
% Slope_r=(Slope_f5_r+Slope_f11_r+Slope_f16_r)/3;
% Slope_g=(Slope_f5_g+Slope_f11_g+Slope_f16_g)/3;
% Slope_b=(Slope_f5_b+Slope_f11_b+Slope_f16_b)/3;

% read in exposure time and aperture:
% for each exp time, to adjust to f5 aperture, multiply by


% Dutch:  I'm commenting this out, and redoing these calculations. -Pat
% fNum=readmein;
% ExpTime=readmein;
% img=readmein;
% 
% ApFac=(5/fNum)^2;
% img=ApFac*img;
% Using F5


% img(:,:,1) = img(:,:,1)/exp(Slope_f5_r*log10(ExpTime)+yInt(1));
% img(:,:,2) = img(:,:,2)/exp(Slope_f5_g*log10(ExpTime)+yInt(2));
% img(:,:,3) = img(:,:,3)/exp(Slope_f5_b*log10(ExpTime)+yInt(3));




fNums = [1.8 2.0 2.2 2.5 2.8 3.2 3.5 4.0 4.5 5.0 5.6 6.3 7.1 8.0 9.0 10.0 11.0 13.0 14.0 16.0 18.0 20.0 22.0];

% % want function that converts fNums and exposureTimes to a pixel value:
% % first, want bounds on exposure times for each f-number:
% % MaxExpTime(10) is the maximum exposure time corresponding to fnum=5,
% % because fNums(10)=5
% 



%%%%%%%%  Pat's Code Starts Here %%%%%%%%%


% ApFac=(5/fNum)^2; Dutch:  I think this should be ApFac=(fNum/5)^2
% img=ApFac*img;

figure; clf; hold on
plot(log10(exposureTimes),logRGBValues_f11(:,1),'r+')
plot(log10(exposureTimes),logRGBValues_f11(:,2),'g+')
plot(log10(exposureTimes),logRGBValues_f11(:,3),'b+')
plot(log10(exposureTimes),logRGBValues_f16(:,1),'ro')
plot(log10(exposureTimes),logRGBValues_f16(:,2),'go')
plot(log10(exposureTimes),logRGBValues_f16(:,3),'bo')
plot(log10(exposureTimes),logRGBValues_f5(:,1),'r*')
plot(log10(exposureTimes),logRGBValues_f5(:,2),'g*')
plot(log10(exposureTimes),logRGBValues_f5(:,3),'b*')

logRGBValues_f11=log10(((11/5)^2)*(10.^(logRGBValues_f11)));
logRGBValues_f16=log10(((16/5)^2)*(10.^(logRGBValues_f16)));

figure; clf; hold on
plot(log10(exposureTimes),logRGBValues_f11(:,1),'r+')
plot(log10(exposureTimes),logRGBValues_f11(:,2),'g+')
plot(log10(exposureTimes),logRGBValues_f11(:,3),'b+')
plot(log10(exposureTimes),logRGBValues_f16(:,1),'ro')
plot(log10(exposureTimes),logRGBValues_f16(:,2),'go')
plot(log10(exposureTimes),logRGBValues_f16(:,3),'bo')
plot(log10(exposureTimes),logRGBValues_f5(:,1),'r*')
plot(log10(exposureTimes),logRGBValues_f5(:,2),'g*')
plot(log10(exposureTimes),logRGBValues_f5(:,3),'b*')

% Linearize Corrected RGB Values
LC_RGB_Values_f5=10.^(logRGBValues_f5);
LC_RGB_Values_f11=10.^(logRGBValues_f11);
LC_RGB_Values_f16=10.^(logRGBValues_f16);

% Calculate RGB Values Per exposure time
RGB_f5_per_exp=LC_RGB_Values_f5./[exposureTimes exposureTimes exposureTimes];
RGB_f11_per_exp=LC_RGB_Values_f11./[exposureTimes exposureTimes exposureTimes];
RGB_f16_per_exp=LC_RGB_Values_f16./[exposureTimes exposureTimes exposureTimes];

figure; clf; hold on
plot(log10(exposureTimes),log10(RGB_f5_per_exp(:,1)),'r+')
plot(log10(exposureTimes),log10(RGB_f5_per_exp(:,2)),'g+')
plot(log10(exposureTimes),log10(RGB_f5_per_exp(:,3)),'b+')
plot(log10(exposureTimes),log10(RGB_f11_per_exp(:,1)),'ro')
plot(log10(exposureTimes),log10(RGB_f11_per_exp(:,2)),'go')
plot(log10(exposureTimes),log10(RGB_f11_per_exp(:,3)),'bo')
plot(log10(exposureTimes),log10(RGB_f16_per_exp(:,1)),'r*')
plot(log10(exposureTimes),log10(RGB_f16_per_exp(:,2)),'g*')
plot(log10(exposureTimes),log10(RGB_f16_per_exp(:,3)),'b*')

% Load Cone Spectral Sensitivity Function
load T_cones_ss2;
wls = SToWls(S_cones_ss2);

% Load Cone Luminance Sensitivity Function
load T_ss2000_Y2

% Load Camera Spectral Sensitivity Function
CameraRGB=load('/Volumes/ColorShare1/CalibrateD70/EXPOSURE_LINEARITY/Data/CameraRGB.txt');
Camera_wls=CameraRGB(:,1);
Camera_SS=CameraRGB(:,2:4)';

% Use Spline Fitting for Camera Wls
Camera_SS=SplineCmf(Camera_wls,Camera_SS,wls);  %% !!! Changed SplineCmf - Ask David about this !!!

% Find Linear Transform from Camera Values to Cone Values
i=find(wls<=400 | wls>=700);
Camera_SS(:,i)=[];
T_cones_ss2(:,i)=[];
wls(i,:)=[];

LMS_Transform = ((Camera_SS')\T_cones_ss2')';
Ypred = LMS_Transform*Camera_SS;
figure;
plot(wls,Ypred');
hold on;
plot(wls,T_cones_ss2,'.');

% Let's see what happens if we compute some cone responses to some surfaces
load sur_nickerson
surfaces = SplineSrf(S_nickerson,sur_nickerson,wls);
load spd_D65
light = SplineSpd(S_D65,spd_D65,wls);
actualLMS = T_cones_ss2*diag(light)*surfaces;
actualRGB = Camera_SS*diag(light)*surfaces;
cameraLMS = LMS_Transform*actualRGB;
figure; clf; hold on
plot(actualLMS(1,:),cameraLMS(1,:),'r+');
plot(actualLMS(2,:),cameraLMS(2,:),'g+');
plot(actualLMS(3,:),cameraLMS(3,:),'b+');
axis('square');
axis([0 max(actualLMS(:)) 0 max(actualLMS(:))]);
plot([0 max(actualLMS(:))], [0 max(actualLMS(:))],'k');


