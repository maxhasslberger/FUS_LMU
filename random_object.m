function [son_obj, pos1] = random_object(pos1, del1, r1, pseudo_pos2, del2, r2, n_sub_objs, dims)

% Determine random params
pos1 = pos1{randi(length(pos1))};
r1 = r1(randi(length(r1)));
del1 = del1(randi(length(del1)));

son_obj = makeBall(dims(1), dims(2), dims(3), pos1(1), pos1(2), pos1(3), r1);
son_obj = son_obj - makeBall(dims(1), dims(2), dims(3), pos1(1), pos1(2), pos1(3), r1 - del1);


% Add sub-objects
n_sub_objs = n_sub_objs(randi(length(n_sub_objs)));
for i = 1:n_sub_objs
    r_sub = r2(randi(length(r2)));
    del_sub = del2(randi(length(del2)));

    pos_sub = pseudo_pos2{randi(length(pseudo_pos2))};
    sub_vec = pos_sub - pos1;
    pos_sub = (r1 - del1 / 2) * sub_vec / norm(sub_vec) + pos1; % Put sub-object into main object

    son_obj = son_obj + makeBall(dims(1), dims(2), dims(3), pos_sub(1), pos_sub(2), pos_sub(3), r_sub);
    son_obj = son_obj - makeBall(dims(1), dims(2), dims(3), pos_sub(1), pos_sub(2), pos_sub(3), r_sub - del_sub);
end

son_obj(son_obj > 1) = 1;

end

    

