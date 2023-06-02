%% Init params
close all;
clear;
% Data

% subj_id = 'theresa';
% subj_id = 'boris';
% subj_id = 'FUN0001';
% subj_id = 'FUN0003';
% subj_id = 'FUN0005';
subj_id = 'FUN0007';

filepath = '../Scans';

% t1_filename = fullfile(filepath, 'FUNtest01t1_T1w_MPRAGE_t1.nii');
% t1_filename = fullfile(filepath, 'boris_t1w.nii');
% t1_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9.nii');
% t1_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9.nii');
% t1_filename = fullfile(filepath, 'FUN0005t1_T1w_MPRAGE_t1.nii');
t1_filename = fullfile(filepath, 'FUN0007T1_T1w_MPRAGE_t1_20230525161748_11.nii');

% ct_filename = fullfile(filepath, 'FUNtest01t1_T1w_pseudoCT.nii');
% ct_filename = fullfile(filepath, 'boris_pct.nii');
% ct_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9_pct.nii');
% ct_filename = fullfile(filepath, 'FUN0003t1_T1w_MPRAGE_t1_20230520133911_9_pct.nii');
% ct_filename = fullfile(filepath, 'FUN0005t1_T1w_pseudoCT.nii');
ct_filename = fullfile(filepath, 'FUN0007T1_T1w_MPRAGE_t1_20230525161748_11_pct.nii');

output_dir = fullfile('../Results');

% Simulation options
acoustic_sim = true;
thermal_sim = false;
% write_localite_file = false;

% Transducer param
transducer = 'CTX500';

% focus_coords_mm_orig = [-27, -21, 40]; % theresa
% focus_coords_mm_orig = [-27, -21, 40]; % theresa new
% focus_coords_mm_orig = [-41, -16, 59]; % boris
% focus_coords_mm_orig = [-41, -6, 39]; % 001 fmri
% focus_coords_mm_orig = [-41, -16, 59]; % 001
% focus_coords_mm_orig = [-8, -15, 21]; % 001 sham
% focus_coords_mm_orig = [-43, -10, 64]; % 003
% focus_coords_mm_orig = [-17, -13, 25]; % 003 sham
% focus_coords_mm_orig = [-25, -22, 6]; % 005
% focus_coords_mm_orig = [-13, -20, -40]; % 005 sham
focus_coords_mm_orig = [-25, -22, 6]; % 007%%%%%
% focus_coords_mm_orig = [-13, -20, -40]; % 007 sham

% offset = [96, 126, 126]; % 001
% offset = [96, 127, 126]; % 003
% offset = [96, 113, 168]; % 005
offset = [96, 114, 127]; % 007

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


%% Tranducer positioning
bowl_coord_axis = [132, 127, 133]; % 001, 003, 005
% bowl_coord_axis = [-1, 128, 128];
[bowl_coords, transducer_angle, skull_offset_mm, add_offset] = get_transducer_position(medium, focus_coords_rel, bowl_coord_axis);
bowl_coords = bowl_coords + (focus_coords - focus_coords_rel);
bowl_coords_mm = ((bowl_coords / 1e-3 .* dxyz) - offset) .* [-1, 1, 1] % mm

transducer_angle
skull_offset_mm

% Get electrical params
focus_depth_mm = ceil(norm(focus_coords / 1e-3 .* dxyz - bowl_coords / 1e-3 .* dxyz)-add_offset); % mm
[pressure, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
focus_depth_mm_NeuroFUS = focus_depth_mm - 13 % Distance to device face
% [pressure, isppa_lut] = get_source_pressure(transducer, phase, focus_depth_mm, isppa_device)

%% Run simulation and store results in output_dir
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth_mm, transducer, ...
    'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, 'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, ...
    'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);

