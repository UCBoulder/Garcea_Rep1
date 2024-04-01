% The main goal of 1A is to find the nuclear mask based on the green channel
% (usually TAG-AF488) signal

clear all; clc;
% pc_server_path = fullfile('Z:','Alex Holland','Matlab Analysis','Garcea_AH230305'); % from a PC accessing the server
% mac_server_path = fullfile('/Volumes/garcea-group','Alex Holland'); % example
% mac_local_path = fullfile('/Users','alex','Desktop','AH230226 iPOND'); % example
pc_local_path = 'C:\Users\alxho\OneDrive\Desktop\iPOND Github';
current_path = pc_local_path; % to make it easy to switch between PC and Mac, local vs server
cd(current_path); 

% Designating input and output folders

% Input_folder with tiff images (cropped using ImageJ)
Input_folder = fullfile(current_path,'Subset AH230226 github INPUT'); 

% Output folder with output images (only for aim1C)
Output_folder = fullfile(current_path,'Subset AH230226 github OUTPUT'); 

% NUCLEAR SEGMENTATION PARAMETERS TO ADJUST
LowAreaBound = 300; % area debris in px^2 to exclude small objects
HighAreaBound = 50000; % area debris in px^2 to exclude large objects
% large objects exclusion not actually needed since the cells were already chosen
% manually
fspecial_radius=6; %ADJUST default=6 for C57 MEFs cells, =2 for WOP cells
is_WOP=0; % if you have WOP cells use =1; for C57 use =0 (default)
% ADJUST manually thresh_bs - see section below
% "Image Overlay TEST" visualize bask_bs1
thresh_bs=600 ; % use about 2xmax value of empty region (image J histogram)


% Choose a few cells by entering the Tiff_name below:
cd(Input_folder); 
Tiff_name = 'r02c06f03p02-ch2sk1fk1fl1-1.tif'; %TAG or green channel to segment since DAPI omitted
% make sure to use .tif or .tiff as appropriate in the file name


Im_double=double(imread(Tiff_name)); % this is our raw image converted to a double

% normlog enhances the image display
normlog_bs=mat2gray(log(Im_double));
% imshow(normlog_bs, []);
% "imtool(normlog_bs,'DisplayRange',[])" in the command window

% binarize using the thresh_bs parameter
mask_bs=imbinarize(Im_double,thresh_bs); % filled
mask_bs1=mask_bs; %test the threshold and vary thresh_bs 

% Alex filtering based on shape
if is_WOP==0
    filter_round=fspecial('disk',fspecial_radius);
    mask_bs=imfilter(mask_bs,filter_round,'same');
    mask_bs2=mask_bs; %test
else
end

if is_WOP==1
    % Good for WOP cells %Active contour (AlexNov2021)
    mask_bs=activecontour(Im_double,mask_bs,'Chan-vese','SmoothFactor',0,'ContractionBias',0);
    mask_bs3=mask_bs; %test
else
end

% Doug Peters size exclusion
% binary, excludes small/large objects, mask_bs here is filled
mask_bs=xor(bwareaopen(mask_bs,LowAreaBound),bwareaopen(mask_bs,HighAreaBound));
% mask_bs=bwareaopen(mask_bs,LowAreaBound); % if cells are already handpicked
mask_bs4=mask_bs; %test

% Fill-in nucleoli
mask_filled=imfill(mask_bs4, 'holes');

% Open the final mask for better visual
mask_open=bwperim(mask_filled);

% Image Overlay TEST (sequential mask_bs1,2,3,4) FILLED
Overlay_filled = imoverlay(normlog_bs, mask_filled, 'red');

%Image Overlay FINAL- OPEN MASK
Overlay_open = imoverlay(normlog_bs, mask_open, 'red');

C=[Overlay_filled Overlay_open];
imshow(C);

% % If you need to troubleshoot each step, uncomment this section
% mask_test = mask_bs4; % input sequentially mask_bs1,2,3,4, mask_filled
% Overlay_open_test = imoverlay(normlog_bs, mask_test, 'green');
% imshow(Overlay_open_test);
