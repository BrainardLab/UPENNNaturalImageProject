# UPENNNaturalImageProject
This contains code used to process the UPENN Natural Image Database.

We brought this into gitHub based on the archives (DHB disk ColorShare3-B) from the work we did when we published the paper and UPENN Natural Image Database.

The database is here: tofu.psych.upenn.edu/~upennidb/. The paper is: Tkacik G et al, "Natural images from the birthplace of the human eye", PLoS ONE 6: e20409 (2011).

We then began working on new things.  Here is information about how things were brought over and where things are, followed by the original README text from 2011.

6/16/17 - DHB

  a) Removed D70Toolbox since it is now part of BrainardLabToolbox.
  b) Needed to install ExifTool to get code to run.  Got this from the web and installed version ExifTool-10.56.  Put the installer into Calibration/CalibrateD70/EXIFTOOLINSTALLER.
  c) Updated a few things to run.  Needed to add an end statement to a one line conditional in D70Toolbox/ProcessPGMToRawMat.  I think this is because Matlab used to treat the “end” as implicitly being there for one line conditionals, but doesn’t any more.
  d) Edited directory strings in process_multipics, and change exifinfo to imfinfo. (Function exifinfo no longer exists in Matlab.)  Also edit so it doesn’t look for picture descriptor files, since we don’t have those and the code that looked for them was a bit brittle about directory names.  This means that that AUX files won’t match what we got when we originally processed the database for the server.
  e) Rename Images -> BotswanaImagesOrig.  These are, I believe, the processed images that we uploaded to the database.
  f) Downloaded cd01A, cd59A and cd17B_closeup_palm_nut_fresh_shade_sun from the database server and them it into BotswanaImagesFromServer_061617/out
  g) Copy over all the Philly raw images to PhillyImages/in and add an input.txt which has the name of each directory. Run process_multipics over the Philly raw images.
  h) Get raw input Botwswana image data and put it into BotswanaImagesCheck_061617/in.  Got two directories from RawImageData/BOTSWANA (DISC1 and DISC59) and called these cd01A and cd59A in the check directory.  Got cd17B_closeup_palm_nut_fresh_shade_sun from BotswanaImagesOrig/in because we don’t seem to have them anywhere else, and deleted all but the original NEF and JPG.  Run process_multipics over these.
  i) Rename Calibration/CalibrateD70 -> CalibrationAndProcessing.
  j) Delete from data archive any folders called xOld.... These remain on the archive disk.
  k) Move the data folders to dropbox UPENNNaturalImage project.
  
2/4/11 - DHB

This directory contains an archive of our natural image database.  Folders as follows

D70Toolbox -- Toolbox as of the time I ran the image processing over the whole database.

Images -- The database directory.  The subfolder "in" has the images from all the CDs, plus the intermediate products of the processing.  The subfolder "out" is what the database server uses.  There are two .csv files in the in subfolder that contain annotations.  The "outppm" folder contains the output .ppm (raw demosaiced) files.  These may be useful for us.  Here is how to read them:
  % In ProcessPGMToRawMat we read these with code:
  %      mosaicPGMName = [filename '.ppm'];
  %      mosaicMonoImage = double(imread(mosaicPGMName));
  %      mosaicMonoImage = mosaicMonoImage(1:2014,1:3038);
  % You don't actually need to chop off the image, but it is what we did.
  % Starting at the upper left, the mosaic is B G B G ...
  %                                           G R G R ...
  %                                           B G B G ...
  %                                           G R G R ...
  %                                              ...

Note that in the database on the web, Gasper removed some images that had people in them, for privacy reasons: cd02A/DSC_0002, cd02A/DSC_0051, cd03A/DSC_0021, cd06A/DSC_0001.

Calibration -- This is the calibration software as well as the calibration image data.  The routine "process_multipics" is what one sets going to process the whole database.  "RunAllCalibrationsAndTests" runs the whole calibration suite and generates the camera calibration files for the D70 toolbox.  This is also under SVN, except I think for the calibration images themselves.  There is a copy of the two subdirectories here (CalibrateD70, Images) in ColorShare3/xArchives/NaturalImageProject.

RawImageData -- The actual images taken in Botswana and Philadelphia, pretty much as they came off the camera.  I also have a copy of this on ColorShare3/xArchives/NaturalImageProject/RawImageData.  I deleted all the old .ppm files from here.

xOldCalibrationCode -- Calibration code from tofu as of January 2011.  Old.  Don't use.  Not sure exactly why I'm keeping it around.  I guess you never know.

xOldMTFImages -- I found these on the tofu server.  It looks like the "pat1" directory contains the set of images we used to calibrate the MTF (NaturalImageDatabase/Images/calibration/MTF/Orig).  The "pat" directory is dated about a week earlier and at least some of the images are not well aligned with the camera in terms of their rotations.  There are a bunch of JPG images, and then NEF images taken right after.  There are only 4 NEF images.  I did not go back and process these, I think they are rightly viewed as a first pass.  I deleted the .ppm files here because heaven knows how they were produced.

xOldOut_122711 -- The out directory from tofu as of 012711.  At least, I'm pretty sure that is what this is.  Unlikely that we would ever want this.

