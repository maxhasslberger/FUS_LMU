%% Input params
close all;
clear;
% Data

% subj_id = 'theresa';
% subj_id = 'boris';
% subj_id = 'FUN0001';
subj_id = 'FUN0005';

filepath = '../Scans';

% t1_filename = fullfile(filepath, 'FUNtest01t1_T1w_MPRAGE_t1.nii');
% t1_filename = fullfile(filepath, 'boris_t1w.nii');
% t1_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9.nii');
t1_filename = fullfile(filepath, 'FUN0005t1_T1w_MPRAGE_t1.nii');

% ct_filename = fullfile(filepath, 'FUNtest01t1_T1w_pseudoCT.nii');
% ct_filename = fullfile(filepath, 'boris_pct.nii');
% ct_filename = fullfile(filepath, 'FUN0001t1_T1w_MPRAGE_20230520084537_9_pct.nii');
ct_filename = fullfile(filepath, 'FUN0005t1_T1w_pseudoCT.nii');
output_dir = fullfile('../Results');

% Simulation options
acoustic_sim = true;
thermal_sim = false;
% write_localite_file = false;

% Transducer param
transducer = 'CTX500';

% focus_coords_mm_orig = [-27, -21, 40]; % mm
% focus_coords_mm_orig = [-41, -16, 59]; % mm
% focus_coords_mm_orig = [-41, -6, 39]; % mm
focus_coords_mm_orig = [-25, -21, 7]; % mm

% Std to structural MRI
% focus_coords_mm_orig = 

isppa_device = 10; % W/cm^2

pulse_length = 20e-3; % s
pulse_rep_freq = 5; % Hz
stim_dur = 80; % s


% Convert into kgrid space
dxyz = [1.0, 1.0, 1.0] * 1e-3; % m
% bowl_coords = bowl_coords_mm * 1e-3 ./ dxyz;
focus_coords_mm = focus_coords_mm_orig .* [-1, 1, 1] + [192, 256, 256] / 2;
focus_coords = focus_coords_mm * 1e-3 ./ dxyz;

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


%% Get grid and medium
[medium, focus_coords_rel, input_ct] = get_medium_param(ct_filename, focus_coords);


%% Tranducer positioning
[bowl_coords, transducer_angle] = get_transducer_position(medium, focus_coords_rel);
bowl_coords = bowl_coords + (focus_coords - focus_coords_rel);
bowl_coords_mm =  ((bowl_coords / 1e-3 .* dxyz) - [192, 256, 256] / 2) .* [-1, 1, 1]; % mm

transducer_angle

% Get electrical params
focus_depth_mm = floor(norm(focus_coords / 1e-3 .* dxyz - bowl_coords / 1e-3 .* dxyz)); % mm
[pressure, phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
focus_depth_mm_NeuroFUS = focus_depth_mm - 13 % Distance to device face
% [pressure, isppa_lut] = get_source_pressure(transducer, phase, focus_depth_mm, isppa_device)

%% Run simulation and store results in output_dir
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth_mm, transducer, ...
    'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, 'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, ...
    'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);

