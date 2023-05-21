function bowl_coords = get_transducer_position(medium, focus_position, stim_angle)

stim_angle = stim_angle/180*pi;
coord_vec = [sin(stim_angle), 0, cos(stim_angle)];

% Get all coordinates on vector
coord_i = focus_position;
prev_coord_x = coord_i;
i = 1;

coordinates = [];
sound_speed = [];
vec_path = zeros(size(medium.sound_speed));
while all(coord_i < size(medium.sound_speed) & coord_i > 1)
    coord_i = round(focus_position + i * coord_vec);
    if sum(prev_coord_x - coord_i) ~= 0
        prev_coord_x = coord_i;
        coordinates = [coordinates; coord_i];
        sound_speed = [sound_speed; medium.sound_speed(coord_i(1), coord_i(2), coord_i(3))];
        vec_path(coord_i(1), coord_i(2), coord_i(3)) = 1;
    end
    i = i+1;
end

% Find point near the skull
thr = 1500; % sound speed threshold (=> water)
offset = 30; % offset points apart from skull

outside_skull_idx = find(sound_speed > thr, 1, 'last');
voxelPlot(double(medium.sound_speed > 1500 | vec_path));

bowl_coords = coordinates(outside_skull_idx + offset, :);

bowl_mask = zeros(size(medium.sound_speed));
bowl_mask(bowl_coords(1), bowl_coords(2), bowl_coords(3)) = 1;
voxelPlot(double(bowl_mask | medium.sound_speed > 1500));

