%% Create driving_param mat files

transducer = 'CTX500';

for focus_depth = 13+50 %47:82
    focus_depth
%     [~, phase] = get_driving_params(focus_depth, transducer); % [Pa, deg]
    phase = phase_correction(focus_depth, transducer); % [Pa, deg]
    get_source_pressure(transducer, phase, focus_depth, 25);
end