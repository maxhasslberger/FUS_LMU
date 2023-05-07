%% Sim params
% Data
subj_id = 'sub-test01';
filepath = '..\Scans';
t1_filename = fullfile(filepath, 'sub-test01_t1w.nii');
ct_filename = fullfile(filepath, 'sub-test01_pct.nii');
output_dir = fullfile('..\Results');
position_dir = fullfile('Localite_pos');

acoustic_sim = true;
thermal_sim = false;

% Transducer param
transducer = 'CTX500';

bowl_coords = [90, 193, 262]; % mm
focus_coords = [99, 161, 202]; % mm
% bowl_coords = [128, 128, 256];
% focus_coords = [128, 128, 128];
focus_depth = 60; % mm

pulse_length = 20e-3; % s
pulse_rep_freq = 5; % Hz
stim_dur = 80; % s
[pressure, phase] = get_driving_params(focus_depth, transducer); % [Pa, deg]

get_space = false;
dxyz = [1.0, 1.0, 1.0] * 1e-3; % m

%% Run simulation to get grid and medium
if get_space
[kgrid, medium] = ...
    tussim_skull_3D('', t1_filename, ct_filename, '', [99, 161, 202], [90, 193, 262], 60, 'CTX500');
dxyz = [kgrid.dx, kgrid.dy, kgrid.dz];
end


% Position / Phase optimization?


%% Write transformation matrix into position file for Localite
t_matrix = get_transducer_transform(bowl_coords*1e-3, focus_coords*1e-3);

position_filename = fullfile(position_dir, 't_position.txt');
if exist(position_filename, 'file') ~= 2
    fopen(position_filename, 'w');
end
dlmwrite(position_filename, t_matrix, 'delimiter', ' ');


%% Run simulation and store results in output_dir
bowl_coords = bowl_coords * 1e-3 ./ dxyz;
focus_coords = focus_coords * 1e-3 ./ dxyz;

[kgrid, medium] = ...
    tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer, 'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, ...
    'PulseLength', pulse_length, 'PulseRepFreq', pulse_rep_freq, 'StimDuration', stim_dur, 'SourcePressure', pressure, 'SourcePhase', phase);


