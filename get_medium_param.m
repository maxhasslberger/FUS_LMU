function [medium, focus_coords, input_ct] = get_medium_param(ct_filename, focus_coords)
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022; Adjusted by
% Max Hasslberger, Technical University of Munich, 21 May 2023
arguments
    ct_filename char
    focus_coords (1,3) {mustBeNumeric}
%     bowl_coords (1,3) {mustBeNumeric}
end
% arguments (Repeating)
%     varargin
% end

%% Simulation parameters
% Please change these if not using 500 kHz transducer.

% medium parameters
c_min               = 1500;     % sound speed [m/s]
c_max               = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
rho_min             = 1000;     % density [kg/m^3]
rho_max             = 1900;     % max. skull density [kg/m3]
alpha_power         = 1.43;     % Robertson et al., PMB 2017 usually between 1 and 3? from Treeby paper
alpha_coeff_water   = 0;        % [dB/(MHz^y cm)] close to 0 (Mueller et al., 2017), see also 0.05 Fomenko et al., 2020?
alpha_coeff_min     = 4;        %
alpha_coeff_max     = 8.7;      % [dB/(MHz cm)] Fry 1978 at 0.5MHz: 1 Np/cm (8.7 dB/cm) for both diploe and outer tables

%% Prepare skull & simulation, check model
%%% Skull preparation
hu_min 	= 300;	% minimum Hounsfield Unit in CT image
hu_max 	= 2000;	% maximum Hounsfield Unit in CT image

% Load CT image (nifti format)
% voxel size = 1 x 1 x 1 mm3, matrix size: varies, mostly 176x256x256
input_ct = niftiread(ct_filename);
% t1_img = niftiread(t1_filename);
% header = niftiinfo(ct_filename);
input_ct = double(input_ct);
input_ct = flip(input_ct, 1); % mirror x axis for voxel space

% voxelPlot(double(input_ct > 1500));

focus_space = zeros(size(input_ct));
focus_space(focus_coords(1), focus_coords(2), focus_coords(3)) = 1;
% voxelPlot(double(input_ct > 1500 | focus_space));

figure;
imagesc(imrotate(squeeze(input_ct(:, focus_coords(2), :) > 500 | focus_space(:, focus_coords(2), :)), 90));
figure;
imagesc(imrotate(squeeze(input_ct(focus_coords(1), :, :) > 500 | focus_space(focus_coords(1), :, :)), 90));
figure;
imagesc(imrotate(permute(squeeze(input_ct(:, :, focus_coords(3)) > 500 | focus_space(:, :, focus_coords(3))), [1 2]), 90));

% update hu_max
ct_max = max(input_ct(:));
if ct_max < hu_max
    hu_max = ct_max;
end
clear ct_max;

% truncate CT HU (see Marsac et al., 2017)
skull_model = input_ct;
skull_model(skull_model < hu_min) = 0; % only use HU for skull acoustic properties
skull_model(skull_model > hu_max) = hu_max;

% pad images by 100 on each side
tmp_model = zeros(size(skull_model,1)+200, size(skull_model,2)+200, size(skull_model,3)+200);
tmp_model(101:size(skull_model,1)+100, ...
    101:size(skull_model,2)+100, ...
    101:size(skull_model,3)+100) = skull_model;
tmp_focus = focus_coords+100;

% centre on focus
% grid size = 256x256x256, new focus coords = [128,128,128]
shift_idx = [tmp_focus(1)-128+1,tmp_focus(1)+128;
    tmp_focus(2)-128+1,tmp_focus(2)+128; ...
    tmp_focus(3)-128+1,tmp_focus(3)+128];
model = tmp_model(shift_idx(1,1):shift_idx(1,2), ...
    shift_idx(2,1):shift_idx(2,2), ...
    shift_idx(3,1):shift_idx(3,2));

shift_x = 128-focus_coords(1);
shift_y = 128-focus_coords(2);
shift_z = 128-focus_coords(3);

% bowl_coords     = bowl_coords + [shift_x, shift_y, shift_z];	% centre of rear surface of transducer [grid points]
focus_coords    = focus_coords + [shift_x, shift_y, shift_z];  % point on the beam axis of the transducer [grid points]

%%% Medium properties
% assign medium properties for skull
% derived from CT HU based on Marsac et al., 2017 & Bancel et al., 2021
medium.density = rho_min + (rho_max - rho_min) * ...
                (model - 0) / (hu_max - 0);
medium.sound_speed = c_min + (c_max - c_min) * ...
                    (medium.density - rho_min) / (rho_max - rho_min);
medium.alpha_coeff = alpha_coeff_min + (alpha_coeff_max - alpha_coeff_min) * ...
                    (1 - (model - hu_min) / (hu_max - hu_min)).^0.5;

% assign medium properties for non-skull (brain, soft tissue, modelled as water)
medium.density(model == 0) = rho_min;
medium.sound_speed(model == 0) = c_min;
medium.alpha_coeff(model == 0) = alpha_coeff_water;

medium.alpha_power = alpha_power;

