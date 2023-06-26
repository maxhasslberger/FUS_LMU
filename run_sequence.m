%% Init params
close all;
clear;

% Select conditions
subj_id = 'theresa';
% % subj_id = 'boris';
% subj_id = 'FUN0001';
% subj_id = 'FUN0002';
% subj_id = 'FUN0003';
% subj_id = 'FUN0004';
% subj_id = 'FUN0005';
% subj_id = 'FUN0006';
% subj_id = 'FUN0007';
% subj_id = 'FUN0009';
% subj_id = 'FUN0008';
% subj_id = 'FUN0010';
% subj_id = 'FUN0011';
% subj_id = 'FUN0012';
% subj_id = 'FUN0013';
% subj_id = 'FUN0014';

sham_cond = false;

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
    pressure = 42272;
    end
    
    offset = [96, 126, 126]; % 001

elseif strcmp(subj_id, 'FUN0002')
    
    t1_filename = fullfile(filepath, 'FUN0002t1_T1w_MPRAGE_t1_20230520110230_9.nii');
    ct_filename = fullfile(filepath, 'FUN0002t1_T1w_MPRAGE_t1_20230520110230_9_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-43, -10, 38];
    focus_coords_mm_orig = [-42, -10, 38];
    else
    focus_coords_mm_orig = [-15, -17, 1];
    pressure = 42272;
    end

    offset = [96, 112, 155];
    bowl_coord_axis_origin = [-68, -24, 78];


elseif strcmp(subj_id, 'FUN0003')
    
    t1_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9.nii');
    ct_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-43, -10, 64]; % 003
    else
    focus_coords_mm_orig = [-13, -10, 26]; % 003 sham
    pressure = 42272;
    end

    offset = [96, 127, 126]; % 003
    bowl_coord_axis_origin = [-73, -18, 101];

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
    focus_coords_mm_orig = [-27, -10, -17];
    focus_coords_mm_orig = [-11, -13, -31];
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

elseif strcmp(subj_id, 'FUN0008')

    t1_filename = fullfile(filepath, 'FUN0008t1_T1w_MPRAGE_t1.nii');
    ct_filename = fullfile(filepath, 'FUN0008t1_T1w_MPRAGE_t1_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-31, -17, 31]; % real
    focus_coords_mm_orig = [-30, -17, 31]; % real
    else
    focus_coords_mm_orig = [-19, -19, -8];
    pressure = 42272;
    end

    offset = [94, 119, 134];
    bowl_coord_axis_origin = [-59, -31, 73];

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
    focus_coords_mm_orig = [-41, -17, 33];
    focus_coords_mm_orig = [-40, -17, 46]; % two different t1 scans (change in z)!
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
    focus_coords_mm_orig = [-33, -22, 38]; % real
    else
    focus_coords_mm_orig = [-1, -16, 9];
    pressure = 42272;
    end

    offset = [96, 126, 125];
    bowl_coord_axis_origin = [-52, -26, 86];

elseif strcmp(subj_id, 'FUN0013')

    t1_filename = fullfile(filepath, 'FUN0013T1_T1w_MPRAGE_.nii');
    ct_filename = fullfile(filepath, 'FUN0013T1_T1w_MPRAGE_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-40, 14, 54]; % real
    focus_coords_mm_orig = [-39, 14, 54]; % real
    else
    focus_coords_mm_orig = [-16, 5, 26];
    pressure = 42272;
    end

    offset = [96, 127, 127];
    bowl_coord_axis_origin = [-67, 0, 96];

elseif strcmp(subj_id, 'FUN0014')

    t1_filename = fullfile(filepath, 'FUN0014T1_FUN0014T1_T1w_MPRAGE_t1_20230620172460_10.nii');
    ct_filename = fullfile(filepath, 'FUN0014T1_T1w_MPRAGE_t1_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-36, -4, 37]; % real
    focus_coords_mm_orig = [-34, -4, 37]; % real
    else
    focus_coords_mm_orig = [-10, -9, 9];
    pressure = 42272;
    end

    offset = [96, 127, 127];
    bowl_coord_axis_origin = [-63, -18, 79];

elseif strcmp(subj_id, 'boris')

    t1_filename = fullfile(filepath, 'boris_t1w.nii');
    ct_filename = fullfile(filepath, 'boris_pct.nii');
    
    if ~sham_cond
    focus_coords_mm_orig = [-27, -21, 40]; % real
    focus_coords_mm_orig = [-25, -21, 40]; % real
    else
    focus_coords_mm_orig = [-8, -19, -1];
    pressure = 42272;
    end

    offset = [96, 127, 127];
    bowl_coord_axis_origin = [-0, -0, 0];

elseif strcmp(subj_id, 'theresa')
    

    t1_filename = fullfile(filepath, 'FUNtest01t1_T1w_MPRAGE.nii');
    ct_filename = fullfile(filepath, 'FUNtest01t1_T1w_MPRAGE_pct.nii');

    if ~sham_cond
    focus_coords_mm_orig = [-27, -21, 40]; % real
    focus_coords_mm_orig = [-27, -21, 40]; % real
    else
    focus_coords_mm_orig = [-8, -19, -1];
    pressure = 42272;
    end

    offset = [96, 109, 131];
    bowl_coord_axis_origin = [-0, -0, 0];
end

if ~sham_cond %%%%%%%%%%%% per subject!
    min_pad_offset = 2; % same for all
    add_offset = 5.5; % FUN12, 4, 10, 2, 8, Theresa  % additional offset as heterogeneous medium deforms focal spot
%     add_offset = 9.5; % FUN11
%     add_offset = 11.5; % FUN9
%     add_offset = 5.0; % FUN6
%     add_offset = 3.; % FUN2
else
    min_pad_offset = 12; % FUN12, 11, 10, 9
%     add_offset = 0; % FUN6
%     add_offset = 5.5; % FUN6
%     add_offset = 10.5; % FUN12, 11, 2, 8
    add_offset = 7.5; % FUN13
    add_offset = 4.5; % FUN13
%     add_offset = 14.5; % FUN10, 9
end

isppa_device = 20; % W/cm^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adj_vec = [1, 1, 1];

% Coordinate transformation into matlab space
focus_coords_mm = focus_coords_mm_orig .* adj_vec + offset;

% Convert into kgrid space
dxyz = [1.0, 1.0, 1.0] * 1e-3; % m
% bowl_coords = bowl_coords_mm * 1e-3 ./ dxyz;
% focus_coords_mm = focus_coords_mm_orig .* adj_vec + [192, 256, 256] / 2;
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
[medium, focus_coords_rel, input_ct] = get_medium_param(t1_filename, ct_filename, focus_coords);

%% Tranducer positioning and param
close all;

if ~sham_cond
%     bowl_coord_axis = [124, 127, 133]; % 001, 003, 005, 007 #############
%     bowl_coord_axis = [124, 127, 136]; % 010
%     bowl_coord_axis = [124, 127, 134]; % 009
    bowl_coord_axis = [124, 126, 134]; % 011, 008, 002, 13, 14
%     bowl_coord_axis = [124, 127, 138]; % 012
%     bowl_coord_axis = [125, 125, 137]; % 006
%     bowl_coord_axis = [125, 127, 137]; % 004
%     bowl_coord_axis = [-1, 128, 128]; % Find angle
else
    bowl_coord_axis = bowl_coord_axis_origin .* adj_vec + offset;
    bowl_coord_axis = bowl_coord_axis - (focus_coords - focus_coords_rel);
%     bowl_coord_axis = [131, 126, 134]; % FUN0006
end

[bowl_coords_rel, transducer_angle, pad_offset_mm, focus_depth] ...
    = get_transducer_position(medium, focus_coords_rel, bowl_coord_axis, min_pad_offset, add_offset);

bowl_coords = bowl_coords_rel + (focus_coords - focus_coords_rel);
bowl_coords_mm = ((bowl_coords / 1e-3 .* dxyz) - offset) .* adj_vec % mm

focus_depth_mm = focus_depth / 1e-3 * dxyz(1); % mm
focus_depth_mm_NeuroFUS = focus_depth_mm - 13 % Distance to device face

if ~sham_cond
%     [pressure, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
    [pressure, phase] = generate_driving_params(focus_depth_mm, transducer, isppa_device); % [Pa, deg]
    pressure
else
%     [~, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
    [~, phase] = generate_driving_params(focus_depth_mm, transducer, isppa_device); % [Pa, deg]
end

transducer_angle
pad_offset_mm

%% Run simulation and store results in output_dir
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth_mm, transducer, ...
    'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, 'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, ...
    'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);

