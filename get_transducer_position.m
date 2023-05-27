function [bowl_coords, opt_angle] = get_transducer_position(medium, focus_position)

% offsets optimized for M1
phi_var = (0 + (-10:5:10)) /180*pi;
theta_var = (45 + (-20:5:20)) /180*pi;

thr = 1500; % sound speed threshold (=> water)

skull_ids_opt = inf;
out_skull_opt = inf;
opt_angle = inf;
vec_path = zeros(size(medium.sound_speed));
coordinates = [];
sound_speed = [];
coord_vec = [];

for theta = theta_var
for phi = phi_var
    coord_vec_tmp = 0.5 * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    
    % Get all coordinates on vector
    coord_i = focus_position;
    prev_coord_x = coord_i;
    i = 1;
    
    coordinates_tmp = [];
    sound_speed_tmp = [];
    vec_path_tmp = zeros(size(medium.sound_speed));
    while all(coord_i < size(medium.sound_speed) & coord_i > 1)
        coord_i = round(focus_position + i * coord_vec_tmp);
        if sum(prev_coord_x - coord_i) ~= 0
            prev_coord_x = coord_i;
            coordinates_tmp = [coordinates_tmp; coord_i];
            sound_speed_tmp = [sound_speed_tmp; medium.sound_speed(coord_i(1), coord_i(2), coord_i(3))];
            vec_path_tmp(coord_i(1), coord_i(2), coord_i(3)) = 1;
        end
        i = i+1;
    end

    out_skull_tmp = find(sound_speed_tmp > thr, 1, 'last'); % crit. 1: focus near skull
    skull_ids_tmp = sum(sound_speed_tmp > 1500); % crit. 2: skull segment as thin as possible
    if ~isempty(out_skull_tmp) && out_skull_tmp > 0 && skull_ids_tmp > 0

        if out_skull_tmp < out_skull_opt || (out_skull_tmp == out_skull_opt && skull_ids_tmp < skull_ids_opt)
            skull_ids_opt = skull_ids_tmp;
            out_skull_opt = out_skull_tmp;
            coordinates = coordinates_tmp;
            sound_speed = sound_speed_tmp;
            coord_vec = coord_vec_tmp;
            vec_path = vec_path_tmp;
            opt_angle = [theta, phi] /pi*180;
        end

    end

end
end

% Find point near the skull
offset = 32; % 15 - offset points apart from skull (=> mm)

out_skull_idx = find(sound_speed > thr, 1, 'last');

% voxelPlot(double(medium.sound_speed > 1500 | vec_path));

bowl_coords = round(coordinates(out_skull_idx, :) + offset * coord_vec / norm(coord_vec));

bowl_mask = zeros(size(medium.sound_speed));
bowl_mask(bowl_coords(1), bowl_coords(2), bowl_coords(3)) = 1;

% voxelPlot(double(bowl_mask | medium.sound_speed > 1500));

