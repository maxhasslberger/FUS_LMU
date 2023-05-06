clearvars;
close all;

% simulation settings
DATA_CAST = 'single';       % set to 'single' or 'gpuArray-single' to speed up computations
MASK_PLANE = 'xy';          % set to 'xy' or 'xz' to generate the beam pattern in different planes

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set total number of grid points not including the PML
Nx = 128;    % [grid points]
Ny = 64;     % [grid points]
Nz = 64;     % [grid points]

% set desired grid size in the x-direction
x = 150e-3;                  % [m]

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
medium.sound_speed = 1540*ones(Nx, Ny, Nz);      % [m/s]
medium.density = 1000*ones(Nx, Ny, Nz);          % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
% medium.BonA = 6;

% create the time array
% t_end = 100e-6;                  % [s]
% kgrid.makeTime(medium.sound_speed, [], t_end);
kgrid.makeTime(medium.sound_speed);
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
% transducer.source_amp = 84250; % [Pa] (72850 calibrated values at 30 W/cm^2 free-water Isppa, 84250 at 40 W/cm^2, 94100 for 50 W/cm^2)
transducer.source_amp = 84250 * ones(1,4); % [Pa] (72850 calibrated values at 30 W/cm^2 free-water Isppa, 84250 at 40 W/cm^2, 94100 for 50 W/cm^2)
transducer.source_phase_deg = [0.0, 307.0, 262.0, 209.0]; % source phase [deg] (calibrated values at 20 W/cm^2 free-water Isppa)
transducer.source_freq_hz = 500e3; % [Hz] the central frequency



%% Position + Phase Optimization
% transducer.source_phase_deg = [0.0, 0.0, 0.0, 0.0];
transducer.source_phase_deg = [0.0, 307.0, 262.0, 209.0];
t_position = round([1, Ny/2, Nz/2]);
t_focuspos = round([Nx, Ny/2, Nz/2]);



  
[source.p_mask, t_label] = transducer_setup(transducer, t_position, t_focuspos, [Nx, Ny, Nz], dx*1000);

% % Setup source
% % create the input signal using toneBurst
% tone_burst_cycles = 5;
% tone_burst_offset = 1.5*1e-6 / kgrid.dt / transducer.n_elements * (1:transducer.n_elements); % Test max delay = 1.5 us
% 
% input_signal = transducer.source_amp(1) * toneBurst(1/kgrid.dt, transducer.source_freq_hz, tone_burst_cycles, ...
%     'SignalOffset', tone_burst_offset);

% Setup source
% define the input signal for each element
input_signal = createCWSignals(kgrid.t_array, transducer.source_freq_hz, transducer.source_amp, transducer.source_phase_deg/180*pi);

p_mask_source_p = t_label;
p_mask_source_p(p_mask_source_p(:,:,:) == 0) = [];
p_mask_source_p = reshape(p_mask_source_p,[],1);

source.p = zeros(length(p_mask_source_p),length(input_signal));

for ii = 1 : length(p_mask_source_p)
    source.p(ii, :) = input_signal(p_mask_source_p(ii), :);
end

% source.p = zeros(sum(source.p_mask(:)), size(input_signal, 2));
% sum_elements = 0;
% for i = 1:transducer.n_elements
%     condition = t_label == i;
%     condition = sum(condition(:));
%     sum_elements = sum_elements + condition;
%     source.p(1+sum_elements:condition+sum_elements, :) = repmat(input_signal(i, :), condition, 1);
% end

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
% source.p = (transducer.source_amp ./ (medium.sound_speed * medium.density)) .* source.p;

% Introduce sonicated object
grid_focus = round(transducer.curv_radius_mm / 1e3 / dx);
position1 = {round([grid_focus, Ny/2, Nz/2])};
% r1 = Ny/2;
r1 = 50;
delta1 = 3;
pseuod_pos = {round([Nx/4, Ny/3, Nz/2])};
r2 = 5;
delta2 = r2;
% n_sub_obj = 0;
n_sub_obj = 1;
son_obj = random_object(position1, delta1, r1, pseuod_pos, delta2, r2, n_sub_obj, [Nx, Ny, Nz]);

medium.sound_speed(son_obj > 0) = 2800;      % [m/s]
medium.density(son_obj > 0) = 1850;          % [kg/m^3]

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

% set the record mode such that only the rms and peak values are stored
sensor.record = {'p_rms', 'p_max'};

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
% input_args = {'DisplayMask', transducer.all_elements_mask, ...
%     'PMLInside', false, 'PlotPML', false, ...
%     'DataCast', DATA_CAST, 'DataRecast', true, 'PlotScale', [-1/2, 1/2] * transducer.source_amp};
input_args = {'PMLInside', false, 'PlotPML', false, 'DisplayMask', source.p_mask, ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotScale', [-1/2, 1/2] * transducer.source_amp(1)};

% Plot Masks
% voxelPlot(double(source.p_mask));
voxelPlot(double(source.p_mask | son_obj));

figure;
imagesc(squeeze(sum(t_label, 1)))
colorbar;

% return;

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
% sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor);

% =========================================================================
% COMPUTE THE BEAM PATTERN
% =========================================================================

% reshape the returned rms and max fields to their original position
sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nj]);
sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);

% plot the beam pattern using the pressure maximum
figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_max * 1e-6);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Total Beam Pattern Using Maximum Of Recorded Pressure');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [MPa]');
axis image;

% plot the beam pattern using the pressure rms
figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_rms * 1e-6);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Total Beam Pattern Using RMS Of Recorded Pressure');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [MPa]');
axis image;
    
