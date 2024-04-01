% The main goal of 1C is to estimate PCC in nucleus region, outside of
% nucleoli, only if the VRC fraction is above area threshold (1% default)
% and output figures to quality check

clear all; clc;

% PATHS
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

% IMAGE DISPLAY PARAMETERS TO ADJUST
background_ch2 = 300; % visualize output for infected vs health cells
background_ch3 = 600; %  visualize output for infected vs health cells
MinAdj_ch2=1*background_ch2/65536;
MinAdj_ch3=1*background_ch3/65536;
cutoff = 0.3; % 0.5 works if signal is dim in the channel; 0.4 for iPOND
gamma = 0.6; % default is 0.1 for immunoplaque; 0.4 for iPOND

% VRC and NUCLEOLI PARAMETERS TO ADJUST
VRC_factor = 2; % VRC > VRC_factor * median_TAG - default =2
Nucleoli_factor = 0.6; % Nucleoli < Nucleoli_factor * median_TAG - default =0.6
VRC_fraction_thresh = 0.01; % fraction of nucleus area above which PCC - default =0.01
% is performed over the 'no hole' region [region excluding nucleoli holes]
% By default 0.02, namely we perform analysis only if the cell has 2% or
% more of the nucleus area as a VRC.

% INPUT FILENAME
% Choose a few cells by entering the Tiff_name ch2 below:
cd(Input_folder); 
ch2_name = 'r03c04f24p06-ch2sk1fk1fl1-1.tif'; % SPECIFY
% ch2_name = 'r02c08f02p04-ch2sk1fk1fl1-1.tif'; % SPECIFY
% ch2_name = 'r02c06f03p02-ch2sk1fk1fl1-1.tif'; % SPECIFY
% make sure to use .tif or .tiff as appropriate in the file name

% finds the corresponding ch3 filename once ch2 is given
ch2_string = string(ch2_name);

ch2_root = extractBetween(ch2_string,1,9);
ch2_root_char = char(ch2_root);
ch2_im_number = extractBetween(ch2_string,26,27);
ch2_im_number_char = char(ch2_im_number);

ch3_pat_char = [ch2_root_char, '*','-ch3sk1fk1fl1', ch2_im_number_char, '.tif'];

Ach3 = dir(ch3_pat_char);
ch3_name = Ach3.name;  % this is the ch3 name, even if z is not the same

% Nucleus segmentation based on ch2 signal
Tiff_name = ch2_name;

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

C1=[Overlay_filled Overlay_open]; % Image display option
%imshow(C1);

% Image overlay nucleus with mask
Nucleus_log = imoverlay(normlog_bs, mask_open, 'blue');

Im_ch2 = imread(ch2_name);
ch2_adjG = imadjust(Im_ch2, [MinAdj_ch2 cutoff], [0 1], gamma);
Overlay_adjG_ch2 = imoverlay(ch2_adjG, mask_open, 'green');

Im_ch3 = imread(ch3_name);
ch3_adjG = imadjust(Im_ch3, [MinAdj_ch3 cutoff], [0 1], gamma);
Overlay_adjG_ch3 = imoverlay(ch3_adjG, mask_open, 'red');

C2 = [Nucleus_log Overlay_adjG_ch2 Overlay_adjG_ch3]; % For image display

% we get TAG stats within nuclear mask
DAPI.cc = bwconncomp(mask_filled,8); % you need this to use regionprops

I.Ch2 = Im_ch2;
I.Ch3 = Im_ch3;

f = 1; % f will span all sites, all wells

RP(f,1).Ch2 = regionprops(DAPI.cc,I.Ch2,'Area','PixelIdxList','PixelValues','MeanIntensity');
RP(f,1).Ch3 = regionprops(DAPI.cc,I.Ch3,'Area','PixelIdxList','PixelValues','MeanIntensity');
%  you need to look at the structure to find out that PixelValues is the 3rd entry


VRC_mask = false(size(Im_ch2, 1), size(Im_ch2, 2)); % initialize the VRC mask as black
PCC_mask = false(size(Im_ch2, 1), size(Im_ch2, 2)); % initialize the PCC mask as black ('no holes')

t = 0; % need to make sure value goes back to zero in-between loops

if ~isempty(DAPI.cc.PixelIdxList) % make sure that there are nuclei in FOV
    for t = 1:length(DAPI.cc.PixelIdxList) % t spans the number of nuclei in site
        if t>1, clearvars is_px_VRC is_px_no_Hole; else end
        ch2.Average(t,1) = RP(f,1).Ch2(t,1).MeanIntensity; % average TAG in nucleus t
        ch2.Median(t,1) = median(RP(f,1).Ch2(t,1).PixelValues); % median TAG in nucleus t
        ch2.Mode(t,1) = mode(RP(f,1).Ch2(t,1).PixelValues); % median TAG in nucleus t
        ch2.VRC_thresh(t,1) = VRC_factor*ch2.Median(t,1); % VRC threshold
        ch2.Hole_thresh(t,1) = Nucleoli_factor*ch2.Median(t,1); % nucleoli threshold
        ch2.area(t,1) = RP(f,1).Ch2(t,1).Area; % area (namely number of pixels) in nucleus t

        for k = 1:length(RP(f,1).Ch2(t,1).PixelValues) % spans all pixels in nucleus t
            if (RP(f,1).Ch2(t,1).PixelValues(k) > ch2.VRC_thresh(t,1))
                is_px_VRC(k) = 1; % is the pixel in the VRC?
                VRC_mask(RP(f,1).Ch2(t,1).PixelIdxList(k)) = true;
            else
                is_px_VRC(k) = 0;
            end

            ch2.VRC_area(t,1) = sum(is_px_VRC(:));
            ch2.VRC_fraction(t,1) = ch2.VRC_area(t,1)/ch2.area(t,1);

            if (RP(f,1).Ch2(t,1).PixelValues(k) > ch2.Hole_thresh(t,1) )
                is_px_no_Hole(k) = 1; % the pixel outside the nucleoli hole
                PCC_mask(RP(f,1).Ch2(t,1).PixelIdxList(k)) = true;
            else
                PCC_mask(RP(f,1).Ch2(t,1).PixelIdxList(k)) = false;
                is_px_no_Hole(k) = 0; % the pixel is within nucleoli hole
            end

        end

        ch2.VRC_area(t,1) = sum(is_px_VRC(:));
        ch2.VRC_fraction(t,1) = ch2.VRC_area(t,1)/ch2.area(t,1);

        PCC_indices = find(is_px_no_Hole == 1);
        PCC_ch2 = RP(f,1).Ch2(t,1).PixelValues(PCC_indices);
        PCC_ch3 = RP(f,1).Ch3(t,1).PixelValues(PCC_indices);

        if (ch2.VRC_fraction(t,1) > VRC_fraction_thresh)
            ch2.PCC(t,1) = corr2(PCC_ch2, PCC_ch3);
            ch2.isPCC(t,1) = 1;
        else
            ch2.PCC(t,1) = NaN;
            ch2.isPCC(t,1) = 0;
        end

    end
else
    ch2.isPCC(:,1) = 0;
    ch2.VRC_fraction = 'no nucleus in FOV';
    ch2.PCC = 'no nucleus in FOV';
end

% Necessary to avoid dim error in case of empty VRC or PCC region
% displays blue colored overlay if no VRC present
% displays red colored overlay if PCC not calculated (VRC area too small)
if VRC_mask == false(size(Im_ch2, 1), size(Im_ch2, 2))
    VRC_mask =  true(size(Im_ch2, 1), size(Im_ch2, 2)); %otherwise dim errors
else  end

if (sum(ch2.isPCC(:,1)) == 0)|(PCC_mask == false(size(Im_ch2, 1), size(Im_ch2, 2)))
    PCC_mask =  true(size(Im_ch2, 1), size(Im_ch2, 2)); %otherwise dim errors
else  end  % displays colored image if no PCC area is present

Holes_mask = mask_filled - PCC_mask; % 'nucleoli holes' are defined as within nucleus AND outside PCC_mask

if Holes_mask == false(size(Im_ch2, 1), size(Im_ch2, 2))
    Holes_mask =  true(size(Im_ch2, 1), size(Im_ch2, 2)); %otherwise dim errors
else  end


Overlay_VRC_ch2 = labeloverlay(ch2_adjG, VRC_mask,'Transparency',0.2,'Colormap','jet'); % 0.1 opaque; 0.7default
Overlay_nucleoli_ch2 = labeloverlay(ch2_adjG, Holes_mask,'Transparency',0.4,'Colormap','jet'); % 0.1 opaque; 0.7default
Overlay_PCC_ch2 = labeloverlay(ch2_adjG, PCC_mask,'Transparency',0,'Colormap','summer'); % 0 opaque; 0.8default ; 1 transparent
Overlay_PCC_ch3 = labeloverlay(ch3_adjG, PCC_mask,'Transparency',0,'Colormap','autumn'); % 0 opaque; 0.9default; 1 transparent
C3 = [Overlay_VRC_ch2 Overlay_PCC_ch2 Overlay_PCC_ch3]; % Image display option

% Final plot to be saved for Quality control for each datapoint

tiledlayout(2,2);
nexttile
imshow(Overlay_VRC_ch2);
title('TAg channel - VRCs');
nexttile
imshow(Overlay_nucleoli_ch2);
title('TAg channel - nucleoli');
nexttile
imshow(Overlay_PCC_ch2);
title('TAg channel - PCC area');
nexttile
imshow(Overlay_PCC_ch3);
title('iPOND channel - PCC area');


% Select image in imshow argument to adjust parameters
% C1 if tweaking 'NUCLEAR SEGMENTATION PARAMETERS TO ADJUST'
% C2 if tweaking 'IMAGE DISPLAY PARAMETERS TO ADJUST'
% C3 if tweaking 'VRC NUCLEOLI PARAMETERS TO ADJUST'

% imshow(C3);



disp('VRC area fractions = ');
ch2.VRC_fraction
disp('PCC no holes = ');
ch2.PCC

cd(Output_folder);
savefig([ch2_root_char,'-snapshot',ch2_im_number_char,'.fig']);
