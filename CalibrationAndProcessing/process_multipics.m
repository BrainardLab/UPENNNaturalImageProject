% process_multipics
%
% A script that processes the image database and calibrates each image.
%
% You need to have dcraw and exiftool installed.
%
% This depends on BrainardLabToolbox, Psychtoolbox, and SimToolbox.  If you
% use the ToolboxToolbox (and we know you want to), just type:
%   tbUseProject('UPENNNaturalImageProject','reset','full');
%
% June 2017.  Not currently producing annotations becasue that code was
% brittle to the structure of the Botswana database input.  Could fix, but
% not important currently.
%
% 02/2008 -- gt created
% 06/2010 -- gt updated
% 06/17/17 -- dhb  update to run on Philly images.
%             dhb  had to replace exifread with imfinfo
% 06/20/17    dhb  use preferences to get paths, add switch for various
%                  things to do.

%% Clear and close
clear; close all;

%% What images to process?
whichProcess = 'botswanaCheck';
switch (whichProcess)
    case 'botswanaDatabase'
        % This directory is our archive of what's on the server.  Don't
        % reprocess unless you're really sure you want to do this.
        error('Do this with great caution');
        dirprefixpath = getpref('UPENNNaturalImageProject','botswanaDatabase');
    case 'phillyDatabase'
        dirprefixpath = getpref('UPENNNaturalImageProject','phillyDatabase');
    case 'botswanaCheck'
         dirprefixpath = getpref('UPENNNaturalImageProject','botswanaCheck');
end

%% Path to dcraw binary.  Be sure to use the right version.
dcrawpath = fullfile(fileparts(mfilename('fullpath')),'TESTDCRAW');

%% Process the following directories, each directory represents a CD (cds)
disp('**************** BEGIN   process_multipics.m ****************')
cds=textread(fullfile(dirprefixpath,'in','input.txt'),'%s');

%% Eliminate those entries that contain a comment character (#) as the first character 
delete_cds = [];
for i=1:length(cds)
    if (strcmp(cds{i}(1),'#'))
        delete_cds = [delete_cds i];
    end
end
cds(delete_cds)=[];

%% Loop over directories
for ncds = 1:length(cds),
    % Say hello
    disp(sprintf('\n\n---------------- BEGIN        CD %d of %d      ----------------\n\n',ncds, length(cds)));

    % Get the input and put directories, create output if needed
    indir=sprintf('%s/in/%s',dirprefixpath, cds{ncds});
    outdir=sprintf('%s/out/%s',dirprefixpath, cds{ncds});
    if (exist(outdir) ~= 7),  disp(sprintf('Creating output directory %s.\n',outdir)); mkdir(outdir); end;

    % Find all the NEF files
    drr = dir(sprintf('%s/*.NEF',indir));
    disp(sprintf('Processing input directory %s', indir));
    
    % Preprocessing -- normalize the file names, build a list of images
    NumPics = length(drr);
    PicNums = zeros(1,NumPics);
    for i=1:NumPics,
        nm = drr(i).name;
        disp(sprintf('Queueing file %s.', nm));
        ix1 = find(nm == '_');
        if (isempty(ix1))
            % the file does not have a DSC prefix, rename it
            ix1=0; 
            movefile(sprintf('%s/%s', indir, drr(i).name),sprintf('%s/DSC_%s.NEF', indir, nm(1:4)));
            movefile(sprintf('%s/%s', indir, strrep(drr(i).name, 'NEF','JPG')),sprintf('%s/DSC_%s.JPG', indir, nm(1:4)));
        end;
        ix2 = find(nm == '.');
        if (isempty(ix1) || isempty(ix2))
            error('Wrong filename format.');
        end
        nmm = nm(ix1+1:ix1+1+3);
        PicNums(i) = int32(str2num(nmm));
    end

    % Process the chosen path here
    ProcessNEFToPGM(1, indir, dcrawpath);
    ProcessPGMToRawMat(1, indir, 0);
    ProcessRawMatToCalibrated(1, indir, 0);
  
    % Postprocessing -- extract metadata, move files to their respective
    % locations
    %
    % June 2017. This was causing me some headaches as the code is brittle with
    % respect to directory names.  We don't need this right now, so it
    % currently doesn't do anything.  The setting of fi = -1 is the step
    % that causes the skip.
    for i_Num=1:NumPics

        % June 2017, had to replace exifread with imfinfo.
        %metad=exifread(sprintf('%s/DSC_%s.JPG', indir, sprintf('%04d', PicNums(i_Num))));
        metad=imfinfo(sprintf('%s/DSC_%s.JPG', indir, sprintf('%04d', PicNums(i_Num))));
 
        clear Image;
        AUX_StringIn=sprintf('%s/%s_AUX.mat', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        load(AUX_StringIn);
        
        Image.metadata=metad;
        Image.cd = cds{ncds};
        Image.name = sprintf('DSC_%s', sprintf('%04d', PicNums(i_Num)));
        Image.timestamp = date();
        
        % Try to find picture descriptors
        %
        % June 2017.  This is where we avoid doing anything.  Would need to
        % handle directory names better in commented out line to make this work.
        %fi=fopen(sprintf('%s/in/cd%c.csv',dirprefixpath,cds{ncds}(5)));
        fi = -1;
        if (fi == -1)
            Image.ann_distance = [];
            Image.ann_object = [];
            Image.ann_tripod = [];
            Image.ann_scene = [];
            Image.ann_time = [];
            Image.ann_light = [];
            Image.ann_date = [];
            Image.ann_notes = [];
        else
            ttx=textscan(fi, '%s%s%s%s%s%s%s%s%s%s','delimiter',',');
            for i=1:length(ttx{1}),
                cdix = sprintf('cd%02d',str2num(char(ttx{1}(i))));
                if (strcmp(cdix, cds{ncds}(1:4)) == 0) continue; end;
                ffq = strfind(char(ttx{2}(i)),'to');
             %   ffq
                if (isempty(ffq))
                    urange = str2num(char(ttx{2}(i)));
                    lrange = str2num(char(ttx{2}(i)));
                else
                    sx = char(ttx{2}(i));
                    lrange = str2num(sx(1:ffq-1));
                    urange = str2num(sx(ffq+2:end));
                end
              %  urange
              %  lrange
                if (PicNums(i_Num) >= lrange & PicNums(i_Num) <= urange)
                    % we have found the descriptor
                    Image.ann_distance = char(ttx{3}(i));
                    Image.ann_object = char(ttx{4}(i));
                    Image.ann_tripod = char(ttx{5}(i));
                    Image.ann_scene = char(ttx{6}(i));
                    Image.ann_time = char(ttx{7}(i));
                    Image.ann_light = char(ttx{8}(i));
                    Image.ann_date = char(ttx{9}(i));
                    Image.ann_notes = char(ttx{10}(i));
                    break;
                else
                    % leave the fields empty
                    Image.ann_distance = [];
                    Image.ann_object = [];
                    Image.ann_tripod = [];
                    Image.ann_scene = [];
                    Image.ann_time = [];
                    Image.ann_light = [];
                    Image.ann_date = [];
                    Image.ann_notes = [];
                end
            end
            clear ttx;
            fclose(fi);
        end
        
        LMS_String=sprintf('%s/%s_LMS.mat', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        LUM_String=sprintf('%s/%s_LUM.mat', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        NEF_String=sprintf('%s/%s.NEF', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        RGB_String=sprintf('%s/%s_RGB.mat', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        AUX_String=sprintf('%s/%s_AUX.mat', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        JPG_String=sprintf('%s/%s.JPG', outdir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
   
        LMS_StringIn=sprintf('%s/%s_LMS.mat', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        LUM_StringIn=sprintf('%s/%s_LUM.mat', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        NEF_StringIn=sprintf('%s/%s.NEF', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        RGB_StringIn=sprintf('%s/%s_RGB.mat', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        JPG_StringIn=sprintf('%s/%s.JPG', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        
        save(AUX_String,'Image');
        copyfile(JPG_StringIn, JPG_String);
        copyfile(NEF_StringIn, NEF_String);
        movefile(LMS_StringIn, LMS_String);
        movefile(LUM_StringIn, LUM_String);
        movefile(RGB_StringIn, RGB_String);
        
        % Delete the .raw.mat and .ppm to conserve space, since they are
        % uniquely regenerated from NEF
        RAW_StringIn=sprintf('%s/%s.raw.mat', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        PPM_StringIn=sprintf('%s/%s.ppm', indir, strcat( 'DSC_', num2str4(PicNums(i_Num))));
        delete(RAW_StringIn);
        delete(PPM_StringIn);
    end
    
    disp(sprintf('\n\n---------------- DONE         CD %d of %d      ----------------\n\n',ncds, length(cds)));

end
