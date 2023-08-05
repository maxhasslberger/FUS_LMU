function [bowl_coords_rel, opt_angle, pad_offset, focus_depth] = ...
    get_transducer_position(medium, focus_coords_rel, bowl_coord_axis, min_pad_offset, add_offset)

if all(bowl_coord_axis > 0)
    % Manual coordinate axis input
    vec = bowl_coord_axis - focus_coords_rel;
    r_norm = norm(vec);
    phi_var = atan2(vec(2), vec(1));
    theta_var = acos(vec(3) / r_norm);
else
    % offsets optimized for M1
    phi_var = (0 + (-10:5:10)) /180*pi;
    theta_var = (45 + (-20:5:20)) /180*pi;
end

thr = 1500; % sound speed threshold (=> water)

skull_ids_opt = inf;
out_skull_opt = inf;
opt_angle = inf;
% vec_path = zeros(size(medium.sound_speed));
coordinates = [];
sound_speed = [];
coord_vec = [];

for theta = theta_var
for phi = phi_var
    coord_vec_tmp = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    coord_vec_tmp = 0.5 * coord_vec_tmp / norm(coord_vec_tmp);
    
    % Get all coordinates on vector
    coord_i = focus_coords_rel;
    prev_coord_x = coord_i;
    i = 1;
    
    coordinates_tmp = [];
    sound_speed_tmp = [];
    vec_path_tmp = zeros(size(medium.sound_speed));
    while all(coord_i < size(medium.sound_speed) & coord_i > 1)
        coord_i = round(focus_coords_rel + i * coord_vec_tmp);
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
t_face_dis = 13;
hair_offset = 3;
min_offset = t_face_dis + min_pad_offset + hair_offset; % distance plane to transducer face + (hair (3) + min. gel pad thickness (2)) (=> mm)

% 003
% min_offset = t_face_dis + 11.5;% distance plane to transducer face + (hair (3) + min. gel pad thickness (2)) (=> mm)
% add_offset = 13;%1.5, 5 % additional offset as heterogeneous medium deforms focal spot

min_NeuroFUS_fd = t_face_dis + 34;
max_NeuroFUS_fd = t_face_dis + 67;
out_skull_idx = find(sound_speed > thr, 1, 'last');

% voxelPlot(double(medium.sound_speed > 1500 | vec_path));


unit_c_vec = coord_vec / norm(coord_vec);

bowl_coords_rel = round(coordinates(out_skull_idx, :) + min_offset * unit_c_vec);
% bowl_coords_rel_real = round(bowl_coords_rel + add_offset * unit_c_vec);

focus_depth_tmp = norm(focus_coords_rel - bowl_coords_rel);
if min_NeuroFUS_fd > focus_depth_tmp
    pad_offset = min_NeuroFUS_fd - focus_depth_tmp; % add gel pad to yield transducer min. focus depth
    bowl_coords_rel = round(bowl_coords_rel + pad_offset * unit_c_vec);
    pad_offset = pad_offset + min_pad_offset;
    focus_depth = min_NeuroFUS_fd;
% elseif max_NeuroFUS_fd < focus_depth_tmp
%     pad_offset = min_pad_offset;
% %     bowl_coords_rel = round(coordinates(out_skull_idx, :) + min_offset * unit_c_vec);
%     focus_depth = max_NeuroFUS_fd;
else
    pad_offset = min_pad_offset;
    focus_depth = min(round(focus_depth_tmp), max_NeuroFUS_fd);
%     focus_depth = min(max(round(norm(focus_coords_rel - bowl_coords_rel))...
%         , min_NeuroFUS_fd), max_NeuroFUS_fd);
end
pad_offset = pad_offset + add_offset;

% Plot 2D views
img1 = double(squeeze(medium.sound_speed(:, focus_coords_rel(2), :) > thr));
img1(bowl_coords_rel(1), bowl_coords_rel(3)) = 2;
img1(focus_coords_rel(1), focus_coords_rel(3)) = 2;

img2 = double(squeeze(medium.sound_speed(:, :, focus_coords_rel(3)) > thr));
img2(bowl_coords_rel(1), bowl_coords_rel(2)) = 2;
img2(focus_coords_rel(1), focus_coords_rel(2)) = 2;

% bowl_mask1 = zeros(size(medium.sound_speed));
% bowl_mask2 = bowl_mask1;
% bowl_mask1(bowl_coords(1), :, bowl_coords(3)) = 1;
% bowl_mask2(bowl_coords(1), bowl_coords(2), :) = 1;
% 
% focus_mask = zeros(size(medium.sound_speed));
% focus_mask(focus_coords(1), focus_coords(2), focus_coords(3)) = 1;

figure;
imagesc(flipud(imrotate(img1, 90)));
xlabel('x')
ylabel('z')
set(gca,'YDir','normal')

figure;
imagesc(flipud(imrotate(img2, 90)));
xlabel('x')
ylabel('y')
set(gca,'YDir','normal')
colorbar;

% voxelPlot(double(bowl_mask | medium.sound_speed > 1500));

end