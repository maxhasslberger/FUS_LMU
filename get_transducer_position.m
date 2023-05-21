function [bowl_coords, opt_angle] = get_transducer_position(medium, focus_position, stim_angle)

stim_angle = stim_angle/180*pi;
stim_angle_var = (25:5:55) /180*pi;
thr = 1500; % sound speed threshold (=> water)

skull_ids_opt = inf;
coordinates = [];
sound_speed = [];
opt_angle = inf;

for stim_angle = stim_angle_var
    coord_vec = [sin(stim_angle), 0, cos(stim_angle)];
    
    % Get all coordinates on vector
    coord_i = focus_position;
    prev_coord_x = coord_i;
    i = 1;
    
    coordinates_tmp = [];
    sound_speed_tmp = [];
    vec_path = zeros(size(medium.sound_speed));
    while all(coord_i < size(medium.sound_speed) & coord_i > 1)
        coord_i = round(focus_position + i * coord_vec);
        if sum(prev_coord_x - coord_i) ~= 0
            prev_coord_x = coord_i;
            coordinates_tmp = [coordinates_tmp; coord_i];
            sound_speed_tmp = [sound_speed_tmp; medium.sound_speed(coord_i(1), coord_i(2), coord_i(3))];
            vec_path(coord_i(1), coord_i(2), coord_i(3)) = 1;
        end
        i = i+1;
    end

    skull_ids = sum(sound_speed_tmp > 1500);
    if skull_ids < skull_ids_opt && skull_ids > 0
        skull_ids_opt = skull_ids;
        coordinates = coordinates_tmp;
        sound_speed = sound_speed_tmp;
        opt_angle = stim_angle/pi*180;
    end

end

% Find point near the skull
offset = 30; % offset points apart from skull

outside_skull_idx = find(sound_speed > thr, 1, 'last');
voxelPlot(double(medium.sound_speed > 1500 | vec_path));

bowl_coords = coordinates(outside_skull_idx + offset, :);

bowl_mask = zeros(size(medium.sound_speed));
bowl_mask(bowl_coords(1), bowl_coords(2), bowl_coords(3)) = 1;
voxelPlot(double(bowl_mask | medium.sound_speed > 1500));

