function transform = get_transducer_transform(bowl_coords, focus_coords)

% Rotation matrix
unit_v0 = [1 0 0];
vec = focus_coords - bowl_coords;

k = cross(unit_v0, vec);
if(norm(k) > 0)
    k = k / norm(k);
end

theta = acos(dot(unit_v0, vec) / norm(vec));
K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K^2);

% Add translation values and complete transform
transform = [R, bowl_coords'; 0, 0, 0, 1];

