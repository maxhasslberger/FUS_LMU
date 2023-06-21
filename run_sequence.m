%% Init params
close all;
clear;

% Select conditions
% % subj_id = 'theresa';
% % subj_id = 'boris';
% subj_id = 'FUN0001';
% subj_id = 'FUN0003';
subj_id = 'FUN0004';
% subj_id = 'FUN0005';
% subj_id = 'FUN0006';
% subj_id = 'FUN0007';
% subj_id = 'FUN0009';
% subj_id = 'FUN0010';
% subj_id = 'FUN0011';
% subj_id = 'FUN0012';

sham_cond = true;

% Simulation options
acoustic_sim = true;
thermal_sim = false;
% write_localite_file = false;

% General param
transducer = 'CTX500';

filepath = '../Scans';
output_dir = fullfile('../Results');

if strcmp(subj_id, 'FUN0001')

    t1_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9.nii');
    ct_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-41, -16, 59]; % 001
    else
    focus_coords_mm_orig = [-8, -15, 21]; % 001 sham
    pressure = 45000;
    end
    
    offset = [96, 126, 126]; % 001

elseif strcmp(subj_id, 'FUN0003')
    
    t1_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9.nii');
    ct_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-43, -10, 64]; % 003
    else
    focus_coords_mm_orig = [-17, -13, 25]; % 003 sham
    pressure = 45000;
    end

    offset = [96, 127, 126]; % 003

elseif strcmp(subj_id, 'FUN0004')

    t1_filename = fullfile(filepath, 'FUN0004t1_T1w_MPRAGE_t1_20230520152450_9.nii');
    ct_filename = fullfile(filepath, 'FUN0004t1_T1w_MPRAGE_t1_20230520152450_9_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-37, -5, 17]; % real
    focus_coords_mm_orig = [-34, -5, 17]; % real
    else
    focus_coords_mm_orig = [-19, -8, -15];
    pressure = 42272;
    end

    offset = [96, 127, 126];
    bowl_coord_axis_origin = [-50, -10, 67];

elseif strcmp(subj_id, 'FUN0005')

    t1_filename = fullfile(filepath, 'FUN0005t1_T1w_MPRAGE_t1.nii');
    ct_filename = fullfile(filepath, 'FUN0005t1_T1w_pseudoCT.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-25, -22, 6]; % 005
    else
    focus_coords_mm_orig = [-13, -20, -40]; % 005 sham
    pressure = 45000;
    end

    offset = [96, 113, 168]; % 005

elseif strcmp(subj_id, 'FUN0006')

    t1_filename = fullfile(filepath, 'FUN0006t1_T1w_MPRAGE_t1.nii');
    ct_filename = fullfile(filepath, 'FUN0006t1_T1w_MPRAGE_t1_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-24, -22, 6]; % real
    focus_coords_mm_orig = [-21, -22, 6]; % real
    else
    focus_coords_mm_orig = [-30, -10, -16];
    pressure = 42272;
    end

    offset = [96, 114, 168];
    bowl_coord_axis_origin = [-36, -37, 53];

elseif strcmp(subj_id, 'FUN0007')

    t1_filename = fullfile(filepath, 'FUN0007T1_T1w_MPRAGE_t1_20230525161748_11.nii');
    ct_filename = fullfile(filepath, 'FUN0007T1_T1w_MPRAGE_t1_20230525161748_11pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-29, -5, 56]; % 007
    else
    focus_coords_mm_orig = [-29, -5, 00]; % 007 sham
    pressure = 42272;
    end

    offset = [96, 114, 127]; % 007
    bowl_coord_axis_origin = [-59, -13, 93]; % 007

elseif strcmp(subj_id, 'FUN0009')

    t1_filename = fullfile(filepath, 'FUN0009t1_mprage_T1w_MPRAGE_t1_20230601151318_10.nii');
    ct_filename = fullfile(filepath, 'FUN0009t1_mprage_T1w_MPRAGE_t1_20230601151318_10pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-33, -4, 28]; % real
    focus_coords_mm_orig = [-29, -4, 28]; % real
    else
    focus_coords_mm_orig = [-9, -12, 24];
    pressure = 42272;
    end

    offset = [99, 116, 142];
    bowl_coord_axis_origin = [-61, -12, 76];

elseif strcmp(subj_id, 'FUN0010')

    t1_filename = fullfile(filepath, 'FUN0010_T1W_MPRAGE_T1_0010_T1w_MPRAGE_t1_20230603074612_10.nii');
    ct_filename = fullfile(filepath, 'FUN0010T1_T1w.nii_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-34, -15, 60]; % real
    focus_coords_mm_orig = [-34, -15, 64]; % real
    else
    focus_coords_mm_orig = [-9, -12, 24];
    pressure = 42272;
    end

    offset = [96, 127, 126];
    bowl_coord_axis_origin = [-56, -22, 111];

elseif strcmp(subj_id, 'FUN0011')

    t1_filename = fullfile(filepath, 'FUN0011_T1W_MPRAGE_T1_0010_T1w_MPRAGE_t1_20230603091615_10.nii');
    ct_filename = fullfile(filepath, 'FUN0011T1_T1w.nii_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-39, -22, 52];
    else
    focus_coords_mm_orig = [-7, -18, 9];
    pressure = 42272;
    end

    offset = [99, 121, 147];
    bowl_coord_axis_origin = [-60, -27, 100];

elseif strcmp(subj_id, 'FUN0012')

    t1_filename = fullfile(filepath, 'FUN0012_T1W_MPRAGE_T1_0010_T1w_MPRAGE_t1_20230603160918_10.nii');
    ct_filename = fullfile(filepath, 'FUN0012T1_T1w.nii_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-34, -22, 38]; % real
    focus_coords_mm_orig = [-32, -22, 38]; % real
    else
    focus_coords_mm_orig = [-1, -16, 9];
    pressure = 42272;
    end

    offset = [96, 126, 125];
    bowl_coord_axis_origin = [-52, -26, 86];

% elseif strcmp(subj_id, 'boris')
% 
%     t1_filename = fullfile(filepath, 'boris_t1w.nii');
%     ct_filename = fullfile(filepath, 'boris_pct.nii');
%     
%     focus_coords_mm_orig = [-41, -16, 59]; % boris
% 
% elseif strcmp(subj_id, 'theresa')
%     
% 
%     t1_filename = fullfile(filepath, 'FUNtest01t1_T1w_MPRAGE_t1.nii');
%     ct_filename = fullfile(filepath, 'FUNtest01t1_T1w_pseudoCT.nii');
% 
%     focus_coords_mm_orig = [-27, -21, 40]; % theresa
end

if ~sham_cond %%%%%%%%%%%% per subject!
    min_pad_offset = 2; % same for all
    add_offset = 5.5; % FUN12, 11, 4(, 10?)  % additional offset as heterogeneous medium deforms focal spot
    add_offset = 9.5; % FUN9
    add_offset = 11.5; % FUN9
%     add_offset = 5.0; % FUN6
else
    min_pad_offset = 12; % FUN12, 11, 10, 9
    add_offset = 10.5; % FUN12, 11
    add_offset = 4.5; % FUN12, 11
%     add_offset = 14.5; % FUN10, 9
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinate transformation into matlab space
focus_coords_mm = focus_coords_mm_orig .* [-1, 1, 1] + offset;

% Convert into kgrid space
dxyz = [1.0, 1.0, 1.0] * 1e-3; % m
% bowl_coords = bowl_coords_mm * 1e-3 ./ dxyz;
% focus_coords_mm = focus_coords_mm_orig .* [-1, 1, 1] + [192, 256, 256] / 2;
focus_coords = focus_coords_mm * 1e-3 ./ dxyz;

pulse_length = 20e-3; % s
pulse_rep_freq = 5; % Hz
stim_dur = 80; % s

% %% Write transformation matrix into position file for Localite
% if write_localite_file
%     position_dir = fullfile('Localite_pos');
%     t_matrix = get_transducer_transform(bowl_coords*dxyz, focus_coords*dxyz);
%     
%     position_filename = fullfile(position_dir, 't_position.txt');
%     if exist(position_filename, 'file') ~= 2
%         fopen(position_filename, 'w');
%     end
%     dlmwrite(position_filename, t_matrix, 'delimiter', ' ');
% end


% Get grid and medium
[medium, focus_coords_rel, input_ct] = get_medium_param(ct_filename, focus_coords);

%% Tranducer positioning and param
close all;

if ~sham_cond
%     bowl_coord_axis = [132, 127, 133]; % 001, 003, 005, 007 #############
%     bowl_coord_axis = [132, 127, 136]; % 010
    bowl_coord_axis = [132, 127, 134]; % 009
%     bowl_coord_axis = [132, 127, 137]; % 011
%     bowl_coord_axis = [132, 127, 138]; % 012
%     bowl_coord_axis = [131, 125, 137]; % 006
%     bowl_coord_axis = [131, 127, 137]; % 004
%     bowl_coord_axis = [-1, 128, 128]; % Find angle
else
    bowl_coord_axis = bowl_coord_axis_origin .* [-1, 1, 1] + offset;
    bowl_coord_axis = bowl_coord_axis - (focus_coords - focus_coords_rel);
end

[bowl_coords_rel, transducer_angle, pad_offset_mm, focus_depth] ...
    = get_transducer_position(medium, focus_coords_rel, bowl_coord_axis, min_pad_offset, add_offset);

bowl_coords = bowl_coords_rel + (focus_coords - focus_coords_rel);
bowl_coords_mm = ((bowl_coords / 1e-3 .* dxyz) - offset) .* [-1, 1, 1] % mm

focus_depth_mm = focus_depth / 1e-3 * dxyz(1); % mm
focus_depth_mm_NeuroFUS = focus_depth_mm - 13 % Distance to device face

if ~sham_cond
    [pressure, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
    pressure
else
    [~, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
end

% [pressure, isppa_lut] = get_source_pressure(transducer, phase, focus_depth_mm, isppa_device)

transducer_angle
pad_offset_mm

%% Run simulation and store results in output_dir
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth_mm, transducer, ...
    'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, 'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, ...
    'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);

