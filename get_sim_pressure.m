%% Create driving_param mat files

transducer = 'CTX500';

for focus_depth = 47:82
    [~, phase] = get_driving_params(focus_depth, transducer); % [Pa, deg]
    isppa_filename = get_source_pressure(transducer, phase, focus_depth);
end