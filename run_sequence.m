%% Input params

% Data
subj_id = 'boris';
filepath = '../Scans';
t1_filename = fullfile(filepath, 'boris_t1w.nii');
ct_filename = fullfile(filepath, 'boris_pct.nii');
output_dir = fullfile('../Results');

% Simulation options
acoustic_sim = true;
thermal_sim = false;
write_localite_file = false;

% Transducer param
transducer = 'CTX500';

% focus_coords = [128, 170, 182]; % mm
% bowl_coords = [128, 170, 250]; % mm
% focus_coords = [99, 161, 202]; % mm
% bowl_coords = [90, 193, 262]; % mm

focus_coords = [130, 130, 150]; % mm
bowl_coords = focus_coords + [70, 0, 0];

isppa_device = 10; % W/cm^2

pulse_length = 20e-3; % s
pulse_rep_freq = 5; % Hz
stim_dur = 80; % s


% Convert into kgrid space
dxyz = [1.0, 1.0, 1.0] * 1e-3; % m
% bowl_coords = bowl_coords * 1e-3 ./ dxyz;
focus_coords = focus_coords * 1e-3 ./ dxyz;

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


%% Run simulation to get grid and medium
[kgrid, medium] = ...
    tussim_skull_3D('', t1_filename, ct_filename, '', focus_coords, bowl_coords, 50, transducer);
% dxyz = [kgrid.dx, kgrid.dy, kgrid.dz];


%% Tranducer positioning
stim_angle = 30; % deg - reference z-axis

bowl_coords = get_transducer_position(medium, focus_coords, stim_angle);

% Get remaining control params
focus_depth = floor(norm(focus_coords / 1e-3 .* dxyz - bowl_coords / 1e-3 .* dxyz)); % mm
[~, phase] = get_driving_params(focus_depth, transducer); % [Pa, deg]
[pressure, isppa_lut] = get_source_pressure(transducer, phase, focus_depth, isppa_device);

%% Run simulation and store results in output_dir
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer, 'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, ...
    'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, 'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);


