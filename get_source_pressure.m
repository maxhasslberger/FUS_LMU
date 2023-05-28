function [pressure, isppa_lut] = get_source_pressure(transducer_id, phases, focus_depth, isppa_device)

if ~strcmp(transducer_id,'CTX500')
    disp("No data available for this transducer");
    return;
end

amp = 66283;%45000;%30000:15000:60000; % TBD: Dynamic stepsize optimization
filename = fullfile('driving_params/', strcat('params_dis_', num2str(focus_depth), 'mm.mat'));

if ~exist(filename, 'file')
% simulation settings
DATA_CAST = 'single';       % set to 'single' or 'gpuArray-single' to speed up computations
MASK_PLANE = 'xy';          % set to 'xy' or 'xz' to generate the beam pattern in different planes

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set total number of grid points not including the PML
Nx = 256;    % [grid points]
Ny = 128;     % [grid points]
Nz = 128;     % [grid points]
% Nx = 128;    % [grid points]
% Ny = 64;     % [grid points]
% Nz = 64;     % [grid points]

% set desired grid size in the x-direction
x = 180e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny, Nz);      % [m/s]
medium.density = 1000*ones(Nx, Ny, Nz);          % [kg/m^3]
% medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
% medium.alpha_power = 1.5;
% medium.BonA = 6;

% create the time array
t_end = 100e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);
% kgrid.makeTime(medium.sound_speed);
% kgrid.dt = checkStability(kgrid, medium);

% =========================================================================
% DEFINE THE INPUT SOURCE
% =========================================================================

% CTX-500 transducer params
transducer.n_elements = 4; % number of elements in the transducer
transducer.Elements_ID_mm = [0, 32.9184, 46.1264, 56.0324];
transducer.Elements_OD_mm = [32.3596, 45.5676, 55.5244, 64.008];
transducer.curv_radius_mm = 63.20; % radius of curvature of the bowl 
transducer.dist_to_plane_mm = 52.38; % distance to the transducer plane from the geometric focus
% [Pa] (72850 calibrated values at 30 W/cm^2 free-water Isppa, 84250 at 40 W/cm^2, 94100 for 50 W/cm^2)
transducer.source_phase_deg = phases; % source phase [deg] (calibrated values at 20 W/cm^2 free-water Isppa)
transducer.source_freq_hz = 500e3; % [Hz] the central frequency



%% Position + Phase Optimization
t_position = round([Nx/10, Ny/2, Nz/2]);
t_focuspos = round([Nx, Ny/2, Nz/2]);
  
[source.p_mask, t_label] = transducer_setup(transducer, t_position, t_focuspos, [Nx, Ny, Nz], dx*1000);

% Setup source
% define the input signal for each element

p_mask_source_p = t_label;
p_mask_source_p(p_mask_source_p(:,:,:) == 0) = [];
p_mask_source_p = reshape(p_mask_source_p,[],1);


% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
switch MASK_PLANE
    case 'xy'
        
        % define mask
        sensor.mask(:, :, Nz/2) = 1;
        
        % store y axis properties        
        Nj = Ny;
        j_vec = kgrid.y_vec;
        j_label = 'y';
        
    case 'xz'
        
        % define mask
        sensor.mask(:, Ny/2, :) = 1;
        
        % store z axis properties
        Nj = Nz;
        j_vec = kgrid.z_vec;
        j_label = 'z';
        
end 
% sensor.mask = ones(Nx, Ny, Nz);
% set the record mode such that only the rms and peak values are stored
% sensor.record = {'p_rms', 'p_max'};
sensor.record = {'p_max'};


% =========================================================================
% RUN THE SIMULATIONS for different pressure amplitudes
% =========================================================================
lut.isppa = zeros(length(amp), 1);
lut.pressure = zeros(length(amp), 1);
lut_idx = 1;

for amp_x = amp

transducer.source_amp = amp_x * ones(1,4);
input_signal = createCWSignals(kgrid.t_array, transducer.source_freq_hz, transducer.source_amp, transducer.source_phase_deg/180*pi);

source.p = zeros(length(p_mask_source_p),length(input_signal));

for i = 1 : length(p_mask_source_p)
    source.p(i, :) = input_signal(p_mask_source_p(i), :);
end


% set the input settings
% input_args = {'PMLInside', true, 'PlotPML', false, 'DisplayMask', source.p_mask, ...
%     'DataCast', DATA_CAST, 'DataRecast', true, 'PlotScale', [-1/2, 1/2] * transducer.source_amp(1)};
input_args = {'PMLInside', true, 'PlotPML', false, 'DisplayMask', source.p_mask, ...
    'PlotScale', [-1/2, 1/2] * transducer.source_amp(1)};

% run the simulation
% sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
% sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor);

% compute Isppa
x_start = round((focus_depth-5)*1e-3/dx);
sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);
[max_pressure, max_idx] = max(sensor_data.p_max(x_start:end, :), [], 'all');
Isppa = max_pressure^2 / (2 * max(medium.density(:)) * max(medium.sound_speed(:))) * 1e-4; % W/cm^2

lut.isppa(lut_idx) = Isppa;
lut.pressure(lut_idx) = amp_x;
lut_idx = lut_idx + 1;

end

% Store in LUT
% save(filename, 'lut');

% =========================================================================
% COMPUTE THE BEAM PATTERN
% =========================================================================

% reshape the returned rms and max fields to their original position

% plot the beam pattern using the pressure maximum
% voxelPlot(double(source.p_mask));

figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_max * 1e-6);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Total Beam Pattern Using Maximum Of Recorded Pressure');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [MPa]');
axis image;

else
    load(filename, 'lut');
end

[~,idx]=min(abs(lut.isppa-isppa_device));
pressure = lut.pressure(idx);
isppa_lut = lut.isppa(idx);

end
    
