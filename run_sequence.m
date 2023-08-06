%% Init params
close all;
clear;

% Select conditions
% subj.id = 'FUNtest01';
% subj.id = 'FUN0000';
% subj.id = 'FUN0001';
subj.id = 'FUN0002';
% subj.id = 'FUN0003';
% subj.id = 'FUN0004';
% subj.id = 'FUN0005';
% subj.id = 'FUN0006';
% subj.id = 'FUN0007';
% subj.id = 'FUN0009';
% subj.id = 'FUN0008';
% subj.id = 'FUN0010';
% subj.id = 'FUN0011';
% subj.id = 'FUN0012';
% subj.id = 'FUN0013';
% subj.id = 'FUN0014';

explorative_sim = true;
% condition = "real";
condition = "sham";

% Simulation options
acoustic_sim = true;
thermal_sim = false;
% write_localite_file = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General param
filepath = '../Scans';
output_dir = fullfile('../Results');

sham_cond = condition == "sham";
disp("Processing subject " + subj.id + "...")
subj_filename = fullfile('subject_init_params/', strcat(subj.id, "_", condition, '.mat'));

if ~exist(subj_filename, 'file') || explorative_sim
%% Modify params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subj.t1_filename = fullfile(filepath, 'FUN0002t1_T1w_MPRAGE_t1_20230520110230_9.nii');
    subj.ct_filename = fullfile(filepath, 'FUN0002t1_T1w_MPRAGE_t1_20230520110230_9_pct.nii');

    subj.offset = [96, 112, 155];

    if sham_cond
        subj.bowl_coord_axis_origin = [-68, -24, 78];
        subj.pressure = 42272;
    end

    % Explorative params
    expl.focus_coords_mm_orig = [-8, -19, -1]; % Desired FT or sham focus position - adjust for large diffractions introduced by the skull
    expl.min_pad_offset = 12; % min. pad offset (added to hair offset (~3 mm))
    expl.add_offset = 6; % additional pad offset after focus_depth computation
    expl.bowl_coord_axis = [124, 126, 134]; % Adjust relative bowl position (focus at (128, 128, 128))
    expl.isppa_device = 20; % W/cm^2 - FDA ISPTA limit = 720 mW/cm^2


    % Fix params
    subj.dxyz = [1.0, 1.0, 1.0] * 1e-3; % m

    subj.transducer = 'CTX500';
    subj.pulse_length = 20e-3; % s
    subj.pulse_rep_freq = 5; % Hz
    subj.stim_dur = 80; % s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if explorative_sim && exist(subj_filename, 'file')
        load(subj_filename, 'subj'); % Discard prev. entered params
    end
    subj.expl = expl; % Overwrite explorative params

    % Coordinate transformation into matlab space
    focus_coords_mm = subj.expl.focus_coords_mm_orig + subj.offset;
    
    % Convert into kgrid space
    subj.focus_coords = focus_coords_mm * 1e-3 ./ subj.dxyz;

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
    [medium, focus_coords_rel, ~] = get_medium_param(subj.t1_filename, subj.ct_filename, subj.focus_coords);

    %% Tranducer positioning and param
%     close all;

    if sham_cond
        subj.expl.bowl_coord_axis = subj.bowl_coord_axis_origin + subj.offset;
        subj.expl.bowl_coord_axis = subj.expl.bowl_coord_axis - (subj.focus_coords - focus_coords_rel);
    end

    [bowl_coords_rel, subj.transducer_angle, subj.pad_offset_mm, focus_depth] ...
        = get_transducer_position(medium, focus_coords_rel, subj.expl.bowl_coord_axis, subj.expl.min_pad_offset, subj.expl.add_offset);
    
    subj.bowl_coords = bowl_coords_rel + (subj.focus_coords - focus_coords_rel);
    
    subj.focus_depth_mm = focus_depth / 1e-3 * subj.dxyz(1); % mm
    
    if ~sham_cond
    %     [subj.pressure, subj.phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
        [subj.pressure, subj.phase] = generate_driving_params(subj.focus_depth_mm, subj.transducer, subj.expl.isppa_device); % [Pa, deg]
    else
        [~, subj.phase] = get_driving_params(subj.focus_depth_mm, subj.transducer); % [Pa, deg]
    %     [~, subj.phase] = generate_driving_params(focus_depth_mm, transducer, isppa_device); % [Pa, deg]
    end
    

    save(subj_filename, 'subj');
else
    load(subj_filename, 'subj');
end

% Display important params
bowl_coords_mm = ((subj.bowl_coords / 1e-3 .* subj.dxyz) - subj.offset); % mm
focus_depth_mm_NeuroFUS = subj.focus_depth_mm - 13; % Distance to device face

disp("Bowl_coords_mm: " + num2str(bowl_coords_mm))
disp("focus_depth_mm_NeuroFUS: " + num2str(focus_depth_mm_NeuroFUS))
disp("k-wave Pressure: " + num2str(subj.pressure))
disp("transducer_angle: " + num2str(subj.transducer_angle))
disp("pad_offset_mm: " + num2str(subj.pad_offset_mm))
fprintf("\n\n")

%% Run simulation and store results in output_dir
tussim_skull_3D(subj.id, subj.t1_filename, subj.ct_filename, output_dir, subj.focus_coords, subj.bowl_coords, subj.focus_depth_mm, subj.transducer, ...
    'RunAcousticSim', acoustic_sim, 'RunThermalSim', thermal_sim, 'PulseLength', subj.pulse_length, 'PulseRepFreq', subj.pulse_rep_freq, ...
    'StimDuration', subj.stim_dur, 'SourcePressure', subj.pressure, 'SourcePhase', subj.phase);

fprintf("\n\n")
