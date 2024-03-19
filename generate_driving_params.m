% Generate all the files:
% generate_driving_params(13 + (34:67), 'CTX500', 0);
% generate_driving_params(13 + [34, 67], 'CTX500', 0);

function [pressure, phase] = generate_driving_params(focus_depths, transducer, isppa_device, amp)

%% Create driving_param mat files
if nargin < 4
%     amp = [42272, 45000];
    amp = 35000:1250:50000;
end

if ~strcmp(transducer,'CTX500')
    disp("No data available for this transducer");
    return;
end

for focus_depth = focus_depths % 13+34 %47:82
    disp("Processing focus depth " + num2str(focus_depth) + " mm")
    filename = fullfile('driving_params/', strcat('params_dis_', num2str(focus_depth), 'mm.mat'));

    if ~exist(filename, 'file')
        [~, phase] = get_driving_params(focus_depth, transducer); % [Pa, deg]
%         phase = phase_correction(focus_depth, transducer); % [Pa, deg]
    
        for lut_idx = 1:length(amp)
            disp(num2str(lut_idx) + "/" + num2str(length(amp)))
%             get_source_pressure(transducer, phase, focus_depth, 25);
            Isppa = tussim_water_3D(transducer, amp(lut_idx), phase);
    
            lut.isppa(lut_idx) = Isppa;
            lut.pressure(lut_idx) = amp(lut_idx);
        end
    
        lut.phase = phase;
        save(filename, 'lut');
        close all;
    else
        load(filename, 'lut');
        phase = lut.phase;
    end

    [~,idx]=min(abs(lut.isppa-isppa_device));
    pressure = lut.pressure(idx);
    clear lut;
end

end