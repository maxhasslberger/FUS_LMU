function phase = phase_correction(focus_depth_mm, transducer_id)

if ~strcmp(transducer_id,'CTX500')
    disp("No data available for this transducer");
    return;
end

%% Define grid

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

%% Define medium

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny, Nz);      % [m/s]
medium.density = 1000*ones(Nx, Ny, Nz);          % [kg/m^3]
% medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
% medium.alpha_power = 1.5;
% medium.BonA = 6;

% create the time array
% t_end = 70e-6;                  % [s]
% kgrid.makeTime(medium.sound_speed, [], t_end);
kgrid.makeTime(medium.sound_speed);

%% Define Sensor

% CTX-500 transducer params
transducer.n_elements = 4; % number of elements in the transducer
transducer.Elements_ID_mm = [0, 32.9184, 46.1264, 56.0324];
transducer.Elements_OD_mm = [32.3596, 45.5676, 55.5244, 64.008];
transducer.curv_radius_mm = 63.20; % radius of curvature of the bowl 
transducer.dist_to_plane_mm = 52.38; % distance to the transducer plane from the geometric focus


%% Position + Phase Optimization
t_position = round([Nx/10, Ny/2, Nz/2]);
t_focuspos = round([Nx, Ny/2, Nz/2]);
  
[sensor.mask, t_label] = transducer_setup(transducer, t_position, t_focuspos, [Nx, Ny, Nz], dx*1000);
sensor.record = {'p'};

% Define the sensor area for each element
mask_sensor_p = t_label;
mask_sensor_p(mask_sensor_p(:,:,:) == 0) = [];
mask_sensor_p = reshape(mask_sensor_p,[],1);


%% Define Source

source.p_mask = zeros(Nx, Ny, Nz);
s_pos = round(t_position + [focus_depth_mm/1000 / dx, 0, 0]);
source.p_mask(s_pos(1), s_pos(2), s_pos(3)) = 1;

% define a single time varying sinusoidal source
source_freq = 0.5e6;           % [MHz]
source_mag = 45000;            % [Pa]
source_pressure = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source_pressure);


%% Run Simulation

% input_args = {'PMLInside', false, 'PlotPML', false, 'DisplayMask', source.p_mask, ...
%     'PlotScale', [-1/2, 1/2] * source_mag};
input_args = {'PMLInside', true, 'PlotPML', false, ...
    'PlotScale', [-1/20, 1/20] * source_mag};

% run the simulation
% sensor_signal = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_signal = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

%% Avg Phase correction
% Cut signal before it arrives at sensor
max_idx = size(sensor_signal.p, 2);
cut_samples = focus_depth_mm/1000 / max(medium.sound_speed(:)) / kgrid.dt;
cut_samples = round((max_idx - cut_samples) / 3 + cut_samples);
sensor_signal.p = sensor_signal.p(:, cut_samples:end);

% FFT
Fs = 1/kgrid.dt;
N = size(sensor_signal.p, 2);
fft_signal = fftshift(fft(sensor_signal.p, [], 2));

% Extract fund. frequency
f = -Fs/2:Fs/N:Fs/2;
[~, idx] = min(abs(f - source_freq));

% Extract phases
all_phases = angle(fft_signal(:, idx));
phase = nan(1, transducer.n_elements);
phase_idx = [1, 420, 347, 277];
% phase_idx = [1, 873, 143, 1567];
% phase_idx = [1, 774, 474, 719]; % TODO: configure appropriate param for Nx=256, Ny=128, Nz=128
phase_idx = [1, 1, 1, 1];
for i = 1:transducer.n_elements
%     phase(i) = mean(all_phases(mask_sensor_p == i));
    phase_elements = all_phases(mask_sensor_p == i);
%     [min_val, min_idx] = min(abs(phase_elements - phase(1)));
    phase(i) = phase_elements(phase_idx(i));
end

phase = (phase - phase(1)) /pi*180; % Ref element 1 = 0°

end