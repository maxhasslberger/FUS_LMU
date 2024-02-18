classdef simulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        Label_2                         matlab.ui.control.Label
        StimulationdurationsEditField   matlab.ui.control.NumericEditField
        StimulationdurationsEditFieldLabel  matlab.ui.control.Label
        PulserepfreqHzEditField         matlab.ui.control.NumericEditField
        PulserepfreqHzEditFieldLabel    matlab.ui.control.Label
        PulselengthsEditField           matlab.ui.control.NumericEditField
        PulselengthsEditFieldLabel      matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        TransducerEditField             matlab.ui.control.EditField
        TransducerEditFieldLabel        matlab.ui.control.Label
        RealSham                        matlab.ui.control.Switch
        SubjectNameEditField            matlab.ui.control.EditField
        SubjectNameEditFieldLabel       matlab.ui.control.Label
        T1wpathEditField                matlab.ui.control.EditField
        T1wpathEditFieldLabel           matlab.ui.control.Label
        CTpathEditField                 matlab.ui.control.EditField
        CTpathEditFieldLabel            matlab.ui.control.Label
        LoadsettingsButton              matlab.ui.control.Button
        InitLabel                       matlab.ui.control.Label
        OtherEditField                  matlab.ui.control.EditField
        OtherEditFieldLabel             matlab.ui.control.Label
        SubjectDropDown                 matlab.ui.control.DropDown
        SubjectDropDownLabel            matlab.ui.control.Label
        CenterPanel                     matlab.ui.container.Panel
        ISPPADeviceWcm2EditField        matlab.ui.control.NumericEditField
        ISPPADeviceWcm2EditFieldLabel   matlab.ui.control.Label
        bowlx                           matlab.ui.control.NumericEditField
        xEditField_7Label               matlab.ui.control.Label
        bowly                           matlab.ui.control.NumericEditField
        yEditField_3Label               matlab.ui.control.Label
        bowlz                           matlab.ui.control.NumericEditField
        zEditField_3Label               matlab.ui.control.Label
        BowlaxisLabel                   matlab.ui.control.Label
        AdditionalpadoffsetmmEditField  matlab.ui.control.NumericEditField
        AdditionalpadoffsetmmEditFieldLabel  matlab.ui.control.Label
        MingelpadoffsetmmEditField      matlab.ui.control.NumericEditField
        MingelpadoffsetmmEditFieldLabel  matlab.ui.control.Label
        focusx                          matlab.ui.control.NumericEditField
        xEditField_6Label               matlab.ui.control.Label
        focusy                          matlab.ui.control.NumericEditField
        yEditField_2Label               matlab.ui.control.Label
        focusz                          matlab.ui.control.NumericEditField
        zEditField_2Label               matlab.ui.control.Label
        FocusCoordsLabel                matlab.ui.control.Label
        offsetx                         matlab.ui.control.NumericEditField
        xEditField_5Label               matlab.ui.control.Label
        offsety                         matlab.ui.control.NumericEditField
        yEditFieldLabel                 matlab.ui.control.Label
        offsetz                         matlab.ui.control.NumericEditField
        zEditFieldLabel                 matlab.ui.control.Label
        SubjectoffsetLabel              matlab.ui.control.Label
        StoresettingsButton             matlab.ui.control.Button
        SubjectParamLabel               matlab.ui.control.Label
        RightPanel                      matlab.ui.container.Panel
        ClearoutputwindowButton         matlab.ui.control.Button
        PreparesimulationButton         matlab.ui.control.Button
        RunsimulationButton             matlab.ui.control.Button
        ThermalSimulationCheckBox       matlab.ui.control.CheckBox
        AcousticSimulationCheckBox      matlab.ui.control.CheckBox
        ResultspathEditField            matlab.ui.control.EditField
        ResultspathEditFieldLabel       matlab.ui.control.Label
        LMUMunichFUSSimulationLabel     matlab.ui.control.Label
        TextArea                        matlab.ui.control.TextArea
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    
    properties (Access = public)
        s = [] % subject variable
        save_subject = true % Description
        output_msg = "" % Description
        simulation_ready = false; % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadsettingsButton
        function LoadsettingsButtonPushed(app, event)
            if ~exist('subject_init_params', 'dir')
                mkdir('subject_init_params')
            end

            if app.SubjectDropDown.Value ~= "Other"
                id = app.SubjectDropDown.Value;
            else
                id = app.OtherEditField.Value;
            end

            name = strcat(id, "_", app.RealSham.Value);
            subj_filename = fullfile('subject_init_params/', strcat(name, '.mat'));

            if exist(subj_filename, 'file')
                load(subj_filename, 'subj');
                app.s = subj;

                app.T1wpathEditField.Value = subj.t1_filename;
                app.CTpathEditField.Value = subj.ct_filename;

                app.offsetx.Value = subj.offset(1);
                app.offsety.Value = subj.offset(2);
                app.offsetz.Value = subj.offset(3);

                app.focusx.Value = subj.expl.focus_coords_mm_orig(1);
                app.focusy.Value = subj.expl.focus_coords_mm_orig(2);
                app.focusz.Value = subj.expl.focus_coords_mm_orig(3);


                app.MingelpadoffsetmmEditField.Value = subj.expl.min_pad_offset; % min. pad offset (added to hair offset (~3 mm))
                app.AdditionalpadoffsetmmEditField.Value = subj.expl.add_offset; % additional pad offset after focus_depth computation
                app.ISPPADeviceWcm2EditField.Value = subj.expl.isppa_device; % W/cm^2 - FDA ISPTA limit = 720 mW/cm^2
            
                app.bowlx.Value = subj.expl.bowl_coord_axis(1);
                app.bowly.Value = subj.expl.bowl_coord_axis(2);
                app.bowlz.Value = subj.expl.bowl_coord_axis(3);
            
                app.TransducerEditField.Value = subj.transducer;
                app.PulselengthsEditField.Value = subj.pulse_length;
                app.PulserepfreqHzEditField.Value = subj.pulse_rep_freq;
                app.StimulationdurationsEditField.Value = subj.stim_dur;

                app.SubjectNameEditField.Value = id;

                msg = "Subject " + name + " loaded" + newline + app.output_msg;
            else
                
                msg = "Subject " + name + " not found!" + newline + app.output_msg;
            end

            app.TextArea.Value = msg;
            app.output_msg = msg;
        end

        % Button pushed function: StoresettingsButton
        function StoresettingsButtonPushed(app, event)
             subj.t1_filename = app.T1wpathEditField.Value;
             subj.ct_filename = app.CTpathEditField.Value;

             subj.offset(1) = app.offsetx.Value;
             subj.offset(2) = app.offsety.Value;
             subj.offset(3) = app.offsetz.Value;

             subj.expl.focus_coords_mm_orig(1) = app.focusx.Value;
             subj.expl.focus_coords_mm_orig(2) = app.focusy.Value;
             subj.expl.focus_coords_mm_orig(3) = app.focusz.Value;


             subj.expl.min_pad_offset = app.MingelpadoffsetmmEditField.Value; % min. pad offset (added to hair offset (~3 mm))
             subj.expl.add_offset = app.AdditionalpadoffsetmmEditField.Value; % additional pad offset after focus_depth computation
             subj.expl.isppa_device = app.ISPPADeviceWcm2EditField.Value; % W/cm^2 - FDA ISPTA limit = 720 mW/cm^2
        
             subj.expl.bowl_coord_axis(1) = app.bowlx.Value;
             subj.expl.bowl_coord_axis(2) = app.bowly.Value;
             subj.expl.bowl_coord_axis(3) = app.bowlz.Value;
        
             subj.transducer = app.TransducerEditField.Value;
             subj.pulse_length = app.PulselengthsEditField.Value;
             subj.pulse_rep_freq = app.PulserepfreqHzEditField.Value;
             subj.stim_dur = app.StimulationdurationsEditField.Value;


             % Coordinate transformation into matlab space
            focus_coords_mm = subj.expl.focus_coords_mm_orig + subj.offset;
            
            % Convert into kgrid space
            subj.dxyz = [1.0, 1.0, 1.0] * 1e-3; % m
            subj.focus_coords = focus_coords_mm * 1e-3 ./ subj.dxyz;

            % %% Write transformation matrix into position file for Localite
            % if write_localite_file
            %     position_dir = fullfile('Localite_pos');
            %     t_matrix = get_transducer_transform(bowl_coords*dxyz, focus_coords*dxyz);
            %     
            %     position_filename = fullfile(position_dir, 't_position.txt');
            %     if exist(position_filename, 'file') ~= 2
            %         fopen(position_filename, 'w');
            %     end
            %     dlmwrite(position_filename, t_matrix, 'delimiter', ' ');
            % end

            % Get grid and medium
            close all;
            [medium, focus_coords_rel, ~] = get_medium_param(subj.t1_filename, subj.ct_filename, subj.focus_coords, ~app.save_subject);
        
            %% Tranducer positioning and param
%             close all;
        
            sham_cond = app.RealSham.Value == "sham";
            if sham_cond
                subj_filename_real = fullfile('subject_init_params/', strcat(app.SubjectNameEditField.Value, '_real.mat'));
                subj_real = load(subj_filename_real, 'subj').subj;
        
                subj.expl.bowl_coord_axis = subj_real.bowl_coords_mm;
                subj.expl.isppa_device = subj_real.expl.isppa_device; % W/cm^2 - FDA ISPTA limit = 720 mW/cm^2
                subj.pressure = subj_real.pressure;
                subj.expl.bowl_coord_axis = subj.expl.bowl_coord_axis + subj.offset;
                subj.expl.bowl_coord_axis = subj.expl.bowl_coord_axis - (subj.focus_coords - focus_coords_rel);
            end
        
            [bowl_coords_rel, subj.transducer_angle, subj.pad_offset_mm, focus_depth] ...
                = get_transducer_position(medium, focus_coords_rel, subj.expl.bowl_coord_axis, subj.expl.min_pad_offset, subj.expl.add_offset, ~app.save_subject);
            
            subj.bowl_coords = bowl_coords_rel + (subj.focus_coords - focus_coords_rel);
            
            subj.focus_depth_mm = focus_depth / 1e-3 * subj.dxyz(1); % mm
            
            if ~sham_cond
            %     [subj.pressure, subj.phase] = get_driving_params(focus_depth_mm, transducer); % [Pa, deg]
                [subj.pressure, subj.phase] = generate_driving_params(subj.focus_depth_mm, subj.transducer, subj.expl.isppa_device); % [Pa, deg]
            else
                [~, subj.phase] = get_driving_params(subj.focus_depth_mm, subj.transducer); % [Pa, deg]
            %     [~, subj.phase] = generate_driving_params(focus_depth_mm, transducer, isppa_device); % [Pa, deg]
            end
        
            subj.bowl_coords_mm = ((subj.bowl_coords / 1e-3 .* subj.dxyz) - subj.offset); % mm
            subj.focus_depth_mm_NeuroFUS = subj.focus_depth_mm - 13; % Distance to device face

            subj.id = app.SubjectNameEditField.Value;

            app.s = subj;

            if app.save_subject
                if app.SubjectDropDown.Value ~= "Other"
                    id = app.SubjectDropDown.Value;
                else
                    id = app.OtherEditField.Value;
                end

                name = strcat(id, "_", app.RealSham.Value);
                subj_filename = fullfile('subject_init_params/', strcat(name, '.mat'));
                save(subj_filename, 'subj');

                msg = "Subject " + strcat(app.SubjectNameEditField.Value, "_", app.RealSham.Value) + " saved" + newline + app.output_msg;
                app.TextArea.Value = msg;
                app.output_msg = msg;
            end
        end

        % Button pushed function: PreparesimulationButton
        function PreparesimulationButtonPushed(app, event)
            app.save_subject = false;
            app.StoresettingsButtonPushed();
            app.save_subject = true;

            msg = ...
                 "Bowl_coords_mm: " + num2str(app.s.bowl_coords_mm) + newline + ...
                 "focus_depth_mm_NeuroFUS: " + num2str(app.s.focus_depth_mm_NeuroFUS) + newline + ...
                 "k-wave Pressure: " + num2str(app.s.pressure) + newline + ...
                 "transducer_angle: " + num2str(app.s.transducer_angle) + newline + ...
                 "pad_offset_mm: " + num2str(app.s.pad_offset_mm) + newline + newline + app.output_msg;

            app.TextArea.Value = msg;
            app.output_msg = msg;

            app.simulation_ready = true;
        end

        % Button pushed function: RunsimulationButton
        function RunsimulationButtonPushed(app, event)
            if app.simulation_ready
            
                msg = app.output_msg + newline + ...
                    "Simulation started - Follow Matlab cmd Window" + newline;
                app.TextArea.Value = msg;
                app.output_msg = msg;
                
                tussim_skull_3D(...
                    app.s.id, ...
                    app.s.t1_filename, ...
                    app.s.ct_filename, ...
                    app.ResultspathEditField.Value, ...
                    app.s.focus_coords, ...
                    app.s.bowl_coords, ...
                    app.s.focus_depth_mm, ...
                    app.s.transducer, ...
                    'RunAcousticSim', app.AcousticSimulationCheckBox.Value, ...
                    'RunThermalSim', app.ThermalSimulationCheckBox.Value, ...
                    'PulseLength', app.s.pulse_length, ...
                    'PulseRepFreq', app.s.pulse_rep_freq, ...
                    'StimDuration', app.s.stim_dur, ...
                    'SourcePressure', app.s.pressure, ...
                    'SourcePhase', app.s.phase);
    
                app.s = [];
    
                msg =  newline + "Simulation completed!" + newline + app.output_msg;
                app.TextArea.Value = msg;
                app.output_msg = msg;

                app.simulation_ready = false;
            else
                app.PreparesimulationButtonPushed();
            end

        end

        % Button pushed function: ClearoutputwindowButton
        function ClearoutputwindowButtonPushed(app, event)
            app.TextArea.Value = "";
            app.output_msg = "";
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {480, 480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {195, '1x', 403};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 860 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {195, '1x', 403};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create SubjectDropDownLabel
            app.SubjectDropDownLabel = uilabel(app.LeftPanel);
            app.SubjectDropDownLabel.HorizontalAlignment = 'right';
            app.SubjectDropDownLabel.Position = [16 411 45 22];
            app.SubjectDropDownLabel.Text = 'Subject';

            % Create SubjectDropDown
            app.SubjectDropDown = uidropdown(app.LeftPanel);
            app.SubjectDropDown.Items = {'FUN0001', 'FUN0002', 'FUN0003', 'FUN0004', 'FUN0005', 'FUN0006', 'FUN0007', 'FUN0008', 'FUN0009', 'FUN0010', 'FUN0011', 'FUN0012', 'FUN0013', 'FUN00014', 'FUN0000', 'FUNtest01', 'FUN0017', 'Other'};
            app.SubjectDropDown.Position = [75 411 100 22];
            app.SubjectDropDown.Value = 'FUN0001';

            % Create OtherEditFieldLabel
            app.OtherEditFieldLabel = uilabel(app.LeftPanel);
            app.OtherEditFieldLabel.HorizontalAlignment = 'right';
            app.OtherEditFieldLabel.Position = [26 379 35 22];
            app.OtherEditFieldLabel.Text = 'Other';

            % Create OtherEditField
            app.OtherEditField = uieditfield(app.LeftPanel, 'text');
            app.OtherEditField.Position = [76 379 100 22];

            % Create InitLabel
            app.InitLabel = uilabel(app.LeftPanel);
            app.InitLabel.FontWeight = 'bold';
            app.InitLabel.FontAngle = 'italic';
            app.InitLabel.Position = [71 447 25 22];
            app.InitLabel.Text = 'Init';

            % Create LoadsettingsButton
            app.LoadsettingsButton = uibutton(app.LeftPanel, 'push');
            app.LoadsettingsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadsettingsButtonPushed, true);
            app.LoadsettingsButton.Position = [51 297 100 23];
            app.LoadsettingsButton.Text = 'Load settings';

            % Create CTpathEditFieldLabel
            app.CTpathEditFieldLabel = uilabel(app.LeftPanel);
            app.CTpathEditFieldLabel.HorizontalAlignment = 'right';
            app.CTpathEditFieldLabel.Position = [24 221 47 22];
            app.CTpathEditFieldLabel.Text = 'CT path';

            % Create CTpathEditField
            app.CTpathEditField = uieditfield(app.LeftPanel, 'text');
            app.CTpathEditField.Position = [86 221 85 22];

            % Create T1wpathEditFieldLabel
            app.T1wpathEditFieldLabel = uilabel(app.LeftPanel);
            app.T1wpathEditFieldLabel.HorizontalAlignment = 'right';
            app.T1wpathEditFieldLabel.Position = [16 248 55 22];
            app.T1wpathEditFieldLabel.Text = 'T1w path';

            % Create T1wpathEditField
            app.T1wpathEditField = uieditfield(app.LeftPanel, 'text');
            app.T1wpathEditField.Position = [86 248 85 22];

            % Create SubjectNameEditFieldLabel
            app.SubjectNameEditFieldLabel = uilabel(app.LeftPanel);
            app.SubjectNameEditFieldLabel.HorizontalAlignment = 'right';
            app.SubjectNameEditFieldLabel.Position = [0 21 80 22];
            app.SubjectNameEditFieldLabel.Text = 'Subject Name';

            % Create SubjectNameEditField
            app.SubjectNameEditField = uieditfield(app.LeftPanel, 'text');
            app.SubjectNameEditField.Position = [95 21 90 22];

            % Create RealSham
            app.RealSham = uiswitch(app.LeftPanel, 'slider');
            app.RealSham.Items = {'real', 'sham'};
            app.RealSham.Position = [80 335 45 20];
            app.RealSham.Value = 'real';

            % Create TransducerEditFieldLabel
            app.TransducerEditFieldLabel = uilabel(app.LeftPanel);
            app.TransducerEditFieldLabel.HorizontalAlignment = 'right';
            app.TransducerEditFieldLabel.Position = [6 163 65 22];
            app.TransducerEditFieldLabel.Text = 'Transducer';

            % Create TransducerEditField
            app.TransducerEditField = uieditfield(app.LeftPanel, 'text');
            app.TransducerEditField.Position = [86 163 85 22];
            app.TransducerEditField.Value = 'CTX-500';

            % Create Label
            app.Label = uilabel(app.LeftPanel);
            app.Label.FontWeight = 'bold';
            app.Label.FontAngle = 'italic';
            app.Label.Position = [3 274 194 22];
            app.Label.Text = '-----------------------------------------------';

            % Create PulselengthsEditFieldLabel
            app.PulselengthsEditFieldLabel = uilabel(app.LeftPanel);
            app.PulselengthsEditFieldLabel.HorizontalAlignment = 'right';
            app.PulselengthsEditFieldLabel.Position = [30 120 88 22];
            app.PulselengthsEditFieldLabel.Text = 'Pulse length (s)';

            % Create PulselengthsEditField
            app.PulselengthsEditField = uieditfield(app.LeftPanel, 'numeric');
            app.PulselengthsEditField.Position = [127 120 44 22];
            app.PulselengthsEditField.Value = 0.02;

            % Create PulserepfreqHzEditFieldLabel
            app.PulserepfreqHzEditFieldLabel = uilabel(app.LeftPanel);
            app.PulserepfreqHzEditFieldLabel.HorizontalAlignment = 'right';
            app.PulserepfreqHzEditFieldLabel.Position = [19 91 109 22];
            app.PulserepfreqHzEditFieldLabel.Text = 'Pulse rep. freq (Hz)';

            % Create PulserepfreqHzEditField
            app.PulserepfreqHzEditField = uieditfield(app.LeftPanel, 'numeric');
            app.PulserepfreqHzEditField.Position = [136 91 36 22];
            app.PulserepfreqHzEditField.Value = 5;

            % Create StimulationdurationsEditFieldLabel
            app.StimulationdurationsEditFieldLabel = uilabel(app.LeftPanel);
            app.StimulationdurationsEditFieldLabel.HorizontalAlignment = 'right';
            app.StimulationdurationsEditFieldLabel.Position = [0 62 128 22];
            app.StimulationdurationsEditFieldLabel.Text = 'Stimulation duration (s)';

            % Create StimulationdurationsEditField
            app.StimulationdurationsEditField = uieditfield(app.LeftPanel, 'numeric');
            app.StimulationdurationsEditField.Position = [136 62 36 22];
            app.StimulationdurationsEditField.Value = 80;

            % Create Label_2
            app.Label_2 = uilabel(app.LeftPanel);
            app.Label_2.FontWeight = 'bold';
            app.Label_2.FontAngle = 'italic';
            app.Label_2.Position = [3 41 186 22];
            app.Label_2.Text = '---------------------------------------------';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create SubjectParamLabel
            app.SubjectParamLabel = uilabel(app.CenterPanel);
            app.SubjectParamLabel.FontWeight = 'bold';
            app.SubjectParamLabel.FontAngle = 'italic';
            app.SubjectParamLabel.Position = [79 446 88 22];
            app.SubjectParamLabel.Text = 'Subject Param';

            % Create StoresettingsButton
            app.StoresettingsButton = uibutton(app.CenterPanel, 'push');
            app.StoresettingsButton.ButtonPushedFcn = createCallbackFcn(app, @StoresettingsButtonPushed, true);
            app.StoresettingsButton.Position = [11 21 100 23];
            app.StoresettingsButton.Text = 'Store settings';

            % Create SubjectoffsetLabel
            app.SubjectoffsetLabel = uilabel(app.CenterPanel);
            app.SubjectoffsetLabel.Position = [2 400 81 22];
            app.SubjectoffsetLabel.Text = 'Subject offset:';

            % Create zEditFieldLabel
            app.zEditFieldLabel = uilabel(app.CenterPanel);
            app.zEditFieldLabel.HorizontalAlignment = 'right';
            app.zEditFieldLabel.Position = [187 400 25 22];
            app.zEditFieldLabel.Text = 'z';

            % Create offsetz
            app.offsetz = uieditfield(app.CenterPanel, 'numeric');
            app.offsetz.Position = [220 400 36 22];
            app.offsetz.Value = 128;

            % Create yEditFieldLabel
            app.yEditFieldLabel = uilabel(app.CenterPanel);
            app.yEditFieldLabel.HorizontalAlignment = 'right';
            app.yEditFieldLabel.Position = [125 400 25 22];
            app.yEditFieldLabel.Text = 'y';

            % Create offsety
            app.offsety = uieditfield(app.CenterPanel, 'numeric');
            app.offsety.Position = [158 400 36 22];
            app.offsety.Value = 128;

            % Create xEditField_5Label
            app.xEditField_5Label = uilabel(app.CenterPanel);
            app.xEditField_5Label.HorizontalAlignment = 'right';
            app.xEditField_5Label.Position = [64 400 25 22];
            app.xEditField_5Label.Text = 'x';

            % Create offsetx
            app.offsetx = uieditfield(app.CenterPanel, 'numeric');
            app.offsetx.Position = [97 400 36 22];
            app.offsetx.Value = 96;

            % Create FocusCoordsLabel
            app.FocusCoordsLabel = uilabel(app.CenterPanel);
            app.FocusCoordsLabel.Position = [2 266 83 22];
            app.FocusCoordsLabel.Text = 'Focus Coords:';

            % Create zEditField_2Label
            app.zEditField_2Label = uilabel(app.CenterPanel);
            app.zEditField_2Label.HorizontalAlignment = 'right';
            app.zEditField_2Label.Position = [187 266 25 22];
            app.zEditField_2Label.Text = 'z';

            % Create focusz
            app.focusz = uieditfield(app.CenterPanel, 'numeric');
            app.focusz.Position = [220 266 36 22];

            % Create yEditField_2Label
            app.yEditField_2Label = uilabel(app.CenterPanel);
            app.yEditField_2Label.HorizontalAlignment = 'right';
            app.yEditField_2Label.Position = [125 266 25 22];
            app.yEditField_2Label.Text = 'y';

            % Create focusy
            app.focusy = uieditfield(app.CenterPanel, 'numeric');
            app.focusy.Position = [158 266 36 22];

            % Create xEditField_6Label
            app.xEditField_6Label = uilabel(app.CenterPanel);
            app.xEditField_6Label.HorizontalAlignment = 'right';
            app.xEditField_6Label.Position = [64 266 25 22];
            app.xEditField_6Label.Text = 'x';

            % Create focusx
            app.focusx = uieditfield(app.CenterPanel, 'numeric');
            app.focusx.Position = [97 266 36 22];

            % Create MingelpadoffsetmmEditFieldLabel
            app.MingelpadoffsetmmEditFieldLabel = uilabel(app.CenterPanel);
            app.MingelpadoffsetmmEditFieldLabel.HorizontalAlignment = 'right';
            app.MingelpadoffsetmmEditFieldLabel.Position = [78 225 134 22];
            app.MingelpadoffsetmmEditFieldLabel.Text = 'Min. gel pad offset (mm)';

            % Create MingelpadoffsetmmEditField
            app.MingelpadoffsetmmEditField = uieditfield(app.CenterPanel, 'numeric');
            app.MingelpadoffsetmmEditField.Position = [220 225 36 22];

            % Create AdditionalpadoffsetmmEditFieldLabel
            app.AdditionalpadoffsetmmEditFieldLabel = uilabel(app.CenterPanel);
            app.AdditionalpadoffsetmmEditFieldLabel.HorizontalAlignment = 'right';
            app.AdditionalpadoffsetmmEditFieldLabel.Position = [67 192 145 22];
            app.AdditionalpadoffsetmmEditFieldLabel.Text = 'Additional pad offset (mm)';

            % Create AdditionalpadoffsetmmEditField
            app.AdditionalpadoffsetmmEditField = uieditfield(app.CenterPanel, 'numeric');
            app.AdditionalpadoffsetmmEditField.Position = [220 192 36 22];

            % Create BowlaxisLabel
            app.BowlaxisLabel = uilabel(app.CenterPanel);
            app.BowlaxisLabel.Position = [2 70 59 22];
            app.BowlaxisLabel.Text = 'Bowl axis:';

            % Create zEditField_3Label
            app.zEditField_3Label = uilabel(app.CenterPanel);
            app.zEditField_3Label.HorizontalAlignment = 'right';
            app.zEditField_3Label.Position = [187 70 25 22];
            app.zEditField_3Label.Text = 'z';

            % Create bowlz
            app.bowlz = uieditfield(app.CenterPanel, 'numeric');
            app.bowlz.Position = [220 70 36 22];

            % Create yEditField_3Label
            app.yEditField_3Label = uilabel(app.CenterPanel);
            app.yEditField_3Label.HorizontalAlignment = 'right';
            app.yEditField_3Label.Position = [125 70 25 22];
            app.yEditField_3Label.Text = 'y';

            % Create bowly
            app.bowly = uieditfield(app.CenterPanel, 'numeric');
            app.bowly.Position = [158 70 36 22];

            % Create xEditField_7Label
            app.xEditField_7Label = uilabel(app.CenterPanel);
            app.xEditField_7Label.HorizontalAlignment = 'right';
            app.xEditField_7Label.Position = [64 70 25 22];
            app.xEditField_7Label.Text = 'x';

            % Create bowlx
            app.bowlx = uieditfield(app.CenterPanel, 'numeric');
            app.bowlx.Position = [97 70 36 22];

            % Create ISPPADeviceWcm2EditFieldLabel
            app.ISPPADeviceWcm2EditFieldLabel = uilabel(app.CenterPanel);
            app.ISPPADeviceWcm2EditFieldLabel.HorizontalAlignment = 'right';
            app.ISPPADeviceWcm2EditFieldLabel.Position = [67 163 127 22];
            app.ISPPADeviceWcm2EditFieldLabel.Text = 'ISPPA Device (W/cm2)';

            % Create ISPPADeviceWcm2EditField
            app.ISPPADeviceWcm2EditField = uieditfield(app.CenterPanel, 'numeric');
            app.ISPPADeviceWcm2EditField.Position = [202 163 54 22];
            app.ISPPADeviceWcm2EditField.Value = 20;

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Create TextArea
            app.TextArea = uitextarea(app.RightPanel);
            app.TextArea.Position = [9 13 385 169];

            % Create LMUMunichFUSSimulationLabel
            app.LMUMunichFUSSimulationLabel = uilabel(app.RightPanel);
            app.LMUMunichFUSSimulationLabel.HorizontalAlignment = 'center';
            app.LMUMunichFUSSimulationLabel.FontSize = 36;
            app.LMUMunichFUSSimulationLabel.FontWeight = 'bold';
            app.LMUMunichFUSSimulationLabel.Position = [65 379 274 89];
            app.LMUMunichFUSSimulationLabel.Text = {'LMU Munich '; 'FUS Simulation'};

            % Create ResultspathEditFieldLabel
            app.ResultspathEditFieldLabel = uilabel(app.RightPanel);
            app.ResultspathEditFieldLabel.HorizontalAlignment = 'right';
            app.ResultspathEditFieldLabel.Position = [9 317 71 22];
            app.ResultspathEditFieldLabel.Text = 'Results path';

            % Create ResultspathEditField
            app.ResultspathEditField = uieditfield(app.RightPanel, 'text');
            app.ResultspathEditField.Position = [95 317 224 22];
            app.ResultspathEditField.Value = '..\Results';

            % Create AcousticSimulationCheckBox
            app.AcousticSimulationCheckBox = uicheckbox(app.RightPanel);
            app.AcousticSimulationCheckBox.Text = 'Acoustic Simulation';
            app.AcousticSimulationCheckBox.Position = [95 274 127 22];
            app.AcousticSimulationCheckBox.Value = true;

            % Create ThermalSimulationCheckBox
            app.ThermalSimulationCheckBox = uicheckbox(app.RightPanel);
            app.ThermalSimulationCheckBox.Text = 'Thermal Simulation';
            app.ThermalSimulationCheckBox.Position = [236 274 125 22];

            % Create RunsimulationButton
            app.RunsimulationButton = uibutton(app.RightPanel, 'push');
            app.RunsimulationButton.ButtonPushedFcn = createCallbackFcn(app, @RunsimulationButtonPushed, true);
            app.RunsimulationButton.Position = [95 192 142 43];
            app.RunsimulationButton.Text = 'Run simulation';

            % Create PreparesimulationButton
            app.PreparesimulationButton = uibutton(app.RightPanel, 'push');
            app.PreparesimulationButton.ButtonPushedFcn = createCallbackFcn(app, @PreparesimulationButtonPushed, true);
            app.PreparesimulationButton.Position = [95 240 142 30];
            app.PreparesimulationButton.Text = 'Prepare simulation';

            % Create ClearoutputwindowButton
            app.ClearoutputwindowButton = uibutton(app.RightPanel, 'push');
            app.ClearoutputwindowButton.ButtonPushedFcn = createCallbackFcn(app, @ClearoutputwindowButtonPushed, true);
            app.ClearoutputwindowButton.Position = [252 184 142 22];
            app.ClearoutputwindowButton.Text = 'Clear output window';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = simulationApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end