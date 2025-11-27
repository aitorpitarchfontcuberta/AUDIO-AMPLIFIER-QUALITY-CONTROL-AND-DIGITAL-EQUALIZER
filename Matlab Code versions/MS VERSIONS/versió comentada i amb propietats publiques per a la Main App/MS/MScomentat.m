classdef MScomentat < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        ImportwavButton          matlab.ui.control.Button
        ImportconfigmatButton    matlab.ui.control.Button
        ExportconfigmatButton    matlab.ui.control.Button
        SavewavButton            matlab.ui.control.Button
        PlayButton               matlab.ui.control.Button
        THDEditField             matlab.ui.control.NumericEditField
        THDLabel                 matlab.ui.control.Label
        GainEditField            matlab.ui.control.NumericEditField
        GainLabel                matlab.ui.control.Label
        PowerEditField           matlab.ui.control.NumericEditField
        PowerLabel               matlab.ui.control.Label
        RecordButton             matlab.ui.control.Button
        FFieldLabel              matlab.ui.control.Label
        FField                   matlab.ui.control.NumericEditField
        DurationField            matlab.ui.control.NumericEditField
        DuracisEditFieldLabel    matlab.ui.control.Label
        SignalTypeDropDown       matlab.ui.control.DropDown
        SignalTypeDropDownLabel  matlab.ui.control.Label
        ModeDropDown             matlab.ui.control.DropDown
        ModeDropDownLabel        matlab.ui.control.Label
        FrequencyResponseAxes    matlab.ui.control.UIAxes
        SpectrogramAxes          matlab.ui.control.UIAxes
        WaveformAxes             matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        FS = 44100;                                                         % Sampling frequency
        y                                                                   % Signal used for analysis/plotting
    end  
    
    methods (Access = private)
        
        function generateSignal(app)                                        % Create a synthetic signal based on type and frequency
            t = 0:1/app.FS:app.DurationField.Value - 1/app.FS;              % Signal duration
            F = app.FField.Value;                                           % Frequency of synthetic signal
            type = app.SignalTypeDropDown.Value;                            % Signal type selection
        
            switch type
                case 'sine'
                    x = sin(2*pi*F*t)';
                case 'square'
                    x = sign(sin(2*pi*F*t))';
                case 'triangle'
                     x = (2 * abs(2 * mod(F * t, 1) - 1) - 1)';                
                case 'chirp'
                    durationch     = 5;      
                    Nbits        = 16;      
                    Fch       = 20000;                              
                    samples = durationch*app.FS;
                    fmi     = (0:1:samples-1)'./(samples-1);      
                    phi       = pi*(Fch/app.FS)*(0:samples-1)';
                    x  = sin(phi.*fmi);                  
                case 'saturated'
                    x = max(min(5*sin(2*pi*F*t), 1), -1)';
            end
        
            app.y = x;
            analyzeAndPlot(app);
        end
    
        function analyzeAndPlot(app)
            % Analyze signal: power, gain, THD
            y = app.y;
        
            % Power
            P = (y'*y)/length(y);
            app.PowerEditField.Value = P;
        
            % Gain
            G = NaN;
            try
                rec = audiorecorder(app.FS, 24, 1);
                recordblocking(rec, 1);
                xamp = getaudiodata(rec);
                G = 10*log10(var(xamp)/var(y));
            end
            app.GainEditField.Value = G;
        
            % THD
            Y = fft(y);                                                     % Compute FFT of the signal
            N = length(Y);                                                  % Length of the FFT
            f = (0:N-1)*(app.FS/N);                                         % Frequency vector
            mag = abs(Y(1:floor(N/2)));                                     % Magnitude spectrum (only positive freqs)
            f_half = f(1:floor(N/2));                                       % Corresponding frequencies
        
            [~, idx_fund] = max(mag);                                       % Find fundamental frequency index
            f_fund = f_half(idx_fund);                                      % Fundamental frequency value
            P1 = mag(idx_fund)^2;                                           % Power of the fundamental component
        
            harmonics = zeros(1, 4);
            for k = 2:5                                                     % Loop over 2nd to 5th harmonics
                target_freq = f_fund * k;
                [~, idx] = min(abs(f_half - target_freq));                  % Find closest freq index
                harmonics(k-1) = mag(idx)^2;                                % Power of each harmonic
            end
        
            if P1 == 0
            THD = -Inf;                                                     % Avoid division by zero
            else
                THD = 10*log10(sum(harmonics) / P1);                        % Compute THD in dB
            end
            app.THDEditField.Value = THD;                                   % Display THD in the app
        
            % Color based on THD
            if isinf(THD)                                                   % Handle infinite THD
                THD_norm = 0;
            elseif THD <= -60
                THD_norm = 0;                                               % Clamp low THD to blue
            else
                THD_norm = (THD + 60) / 60;                                 % Normalize THD from -60 dB to 0 dB
            end
            color_signal = (1 - THD_norm)*[0 0 1] + THD_norm*[1 0 0];       % Interpolate color (blue to red)
        
            % Waveform creation
            plot(app.WaveformAxes, (1:length(y))/app.FS, y, 'Color', color_signal, 'LineWidth', 1.5);
            title(app.WaveformAxes, sprintf('Forma d''ona (THD = %.2f dB)', THD));  % Add THD value to title
            xlabel(app.WaveformAxes, 'Temps (s)');
            ylabel(app.WaveformAxes, 'Amplitud');
        
            % Manual spectrogram
            NFFT = 1024;                                                    % FFT size
            overlap = 0.25 * NFFT;                                          % 25% overlap between frames
            hop = NFFT - overlap;                                           % Hop size
            win = 0.54 - 0.46 * cos(2 * pi * (0:NFFT-1)' / (NFFT - 1));     % Hamming window
            win = win / sqrt(sum(win.^2));                                  % Normalize window energy
            
            % If signal is too short, clear axes and show warning
            if length(y) < NFFT
                cla(app.SpectrogramAxes);
                title(app.SpectrogramAxes, 'Senyal massa curta');
                return;
            end
        
            numFrames = floor((length(y) - NFFT) / hop) + 1;                % Number of frames
            S = zeros(NFFT/2, numFrames);                                   % Initialize spectrogram matrix
            T = zeros(1, numFrames);                                        % Time vector for each frame
        
            % Compute FFT for each frame
            for i = 1:numFrames
                startIdx = (i-1)*hop + 1;
                endIdx = startIdx + NFFT - 1;
                frame = y(startIdx:endIdx, 1) .* win;                       % Apply window to frame
                Y_frame = fft(frame, NFFT);                                 % Compute FFT
                S(:,i) = abs(Y_frame(1:NFFT/2));                            % Keep positive frequencies
                T(i) = startIdx / app.FS;                                   % Time in seconds
            end
        
            f = app.FS * (0:(NFFT/2)-1) / NFFT;                             % Frequency vector
        
            % Spectrogram Creation
            surf(app.SpectrogramAxes, T, f, 20*log10(S + eps), 'EdgeColor', 'none');
            axis(app.SpectrogramAxes, 'tight');
            view(app.SpectrogramAxes, [0, 90]);
            xlabel(app.SpectrogramAxes, 'Temps (s)');
            ylabel(app.SpectrogramAxes, 'Freqüència (Hz)');
            title(app.SpectrogramAxes, 'Espectrograma manual');
            colorbar(app.SpectrogramAxes);

            % Frequency Response
            y = app.y;
            Fs = app.FS;
            
            N = length(y);
            n = (0:N-1)';
            w = 0.5 * (1 - cos(2*pi*n/(N-1)));                              % Hann window (manual implementation)
            
            Y = fft(y .* w);                                                % Apply window and compute FFT
            f = (0:N-1)*(Fs/N);                                             % Frequency vector
            
            % Keep only positive frequencies
            half_N = floor(N/2);
            f = f(1:half_N);
            Y = Y(1:half_N);
            
            % Logarithmic scaling (dB)
            Y_mag = abs(Y);
            Y_mag = Y_mag / max(Y_mag);                                     % Normalize magnitude
            Y_dB = 20*log10(Y_mag + eps);                                   % Convert to decibels
            
            % Plot frequency response
            plot(app.FrequencyResponseAxes, f, Y_dB, 'LineWidth', 1.5);
            xlabel(app.FrequencyResponseAxes, 'Freqüència (Hz)');
            ylabel(app.FrequencyResponseAxes, 'Mòdul (dB)');
            title(app.FrequencyResponseAxes, 'Resposta Freqüencial');
            grid(app.FrequencyResponseAxes, 'on');
            xlim(app.FrequencyResponseAxes, [0 Fs/2]);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: ModeDropDown
        function ModeDropdownValueChanged(app, event)
            value = app.ModeDropDown.Value;
            if strcmp(app.ModeDropDown.Value, "synthetic")
                % Show signal type and frequency fields, hide record button
                app.SignalTypeDropDown.Visible = 'on';
                app.SignalTypeDropDownLabel.Visible = 'on';
                app.FField.Visible = 'on';
                app.FFieldLabel.Visible = 'on';
                app.RecordButton.Visible = 'off';
            else
                % Hide signal type and frequency fields, show record button
                app.SignalTypeDropDown.Visible = 'off';
                app.SignalTypeDropDownLabel.Visible = 'off';
                app.FField.Visible = 'off';
                app.FFieldLabel.Visible = 'off';
                app.RecordButton.Visible = 'on';
            end
        end

        % Value changed function: SignalTypeDropDown
        function SignalDropdownValueChanged(app, event)
            value = app.SignalTypeDropDown.Value;                           % Get selected signal type
            generateSignal(app);                                            % Generate signal based on selection
        end

        % Value changed function: DurationField
        function DurationFieldValueChanged(app, event)
            value = app.DurationField.Value;                                % Get new duration value
            minVal = 0.1;                                                   % Minimum allowed duration
        maxVal = 10;                                                        % Maximum allowed duration
            val = app.DurationField.Value;
        
            % Enforce limits on duration value
            if val < minVal
                val = minVal;
            elseif val > maxVal
                val = maxVal;
            end
        
            app.DurationField.Value = val;                                  % Update field with corrected value
            generateSignal(app);                                            % Regenerate signal with new duration
        end

        % Value changed function: FField
        function FFieldValueChanged(app, event)
            value = app.FField.Value;                                       % Get new frequency value
            minVal = 20;                                                    % Minimum allowed frequency (Hz)
            maxVal = 20000;                                                 % Maximum allowed frequency (Hz)
            val = app.FField.Value;
        
            % Enforce limits on frequency value
            if val < minVal
                val = minVal;
            elseif val > maxVal
                val = maxVal;
            end
        
            app.FField.Value = val;                                         % Update field with corrected value
            generateSignal(app);                                            % Regenerate signal with new frequency
        end

        % Button pushed function: RecordButton
        function RecordButtonPushed(app, event)
            duration = app.DurationField.Value;                             % Get recording duration from UI

            % Create temporary label to show countdown
            countdownLabel = uilabel(app.UIFigure);
            countdownLabel.Position = [487 510 120 22];
            countdownLabel.HorizontalAlignment = 'center';
            countdownLabel.FontSize = 14;
            countdownLabel.Text = 'Iniciant gravació en...';
        
            % 3-second countdown before recording
            for i = 3:-1:1
                countdownLabel.Text = sprintf('Gravant en %d...', i);
                pause(1);
            end
        
            countdownLabel.Text = 'Gravant...';
        
            % Start audio recording with specified sample rate, 24 bits, mono
            rec = audiorecorder(app.FS, 24, 1);
            record(rec);
        
            % Update label every second with remaining recording time
            for i = 1:duration
                countdownLabel.Text = sprintf('Temps restant: %d s', duration - i + 1);
                pause(1);
            end
        
            % Stop recording and update label
            stop(rec);
            countdownLabel.Text = 'Gravació completada';
            pause(0.5);
            delete(countdownLabel);                                         % Remove temporary label
        
            % Retrieve recorded audio data and remove initial silence below threshold
            y = getaudiodata(rec);
            threshold = 0.02;
            idx = find(abs(y) > threshold, 1, 'first');
            if ~isempty(idx)
                y = y(idx:end);
            end
            app.y = y;                                                      % Store audio in app property
            analyzeAndPlot(app);                                            % Analyze and update plots
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
            % PlayButton callback: play the stored audio signal
            sound(app.y, app.FS);
        end

        % Button pushed function: SavewavButton
        function SaveButtonPushed(app, event)
            % SavewavButton callback: save current audio to WAV file
            [file, path] = uiputfile('*.wav', 'Desar com');
            if isequal(file, 0), return; end
            audiowrite(fullfile(path, file), app.y, app.FS);
        end

        % Button pushed function: ImportwavButton
        function ImportAudioButtonPushed(app, event)
            % ImportwavButton callback: load audio from WAV file and analyze
            [file, path] = uigetfile('*.wav');
            if isequal(file,0), return; end
            [y, fs] = audioread(fullfile(path, file));
            app.FS = fs;
            app.y = y;
            analyzeAndPlot(app);
        end

        % Button pushed function: ExportconfigmatButton
        function ExportButtonPushed(app, event)
            % ExportconfigmatButton callback: save current app settings to a MAT file
            [file, path] = uiputfile('*.mat');
            if isequal(file,0), return; end
            config = struct();
            config.mode = app.ModeDropDown.Value;
            config.signal = app.SignalTypeDropDown.Value;
            config.duration = app.DurationField.Value;
            config.F = app.FField.Value;
            save(fullfile(path, file), '-struct', 'config');
        end

        % Button pushed function: ImportconfigmatButton
        function ImportButtonPushed(app, event)
            % ImportconfigmatButton callback: load app settings from a MAT file and update UI
            [file, path] = uigetfile('*.mat');
            if isequal(file,0), return; end
            config = load(fullfile(path, file));
            app.ModeDropDown.Value = config.mode;
            app.SignalTypeDropDown.Value = config.signal;
            app.DurationField.Value = config.duration;
            app.FField.Value = config.F;
            generateSignal(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 638 615];
            app.UIFigure.Name = 'MATLAB App';

            % Create WaveformAxes
            app.WaveformAxes = uiaxes(app.UIFigure);
            title(app.WaveformAxes, 'Waveform')
            xlabel(app.WaveformAxes, 'X')
            ylabel(app.WaveformAxes, 'Y')
            zlabel(app.WaveformAxes, 'Z')
            app.WaveformAxes.Position = [1 327 300 185];

            % Create SpectrogramAxes
            app.SpectrogramAxes = uiaxes(app.UIFigure);
            title(app.SpectrogramAxes, 'Spectrogram')
            xlabel(app.SpectrogramAxes, 'X')
            ylabel(app.SpectrogramAxes, 'Y')
            zlabel(app.SpectrogramAxes, 'Z')
            app.SpectrogramAxes.Position = [322 327 300 185];

            % Create FrequencyResponseAxes
            app.FrequencyResponseAxes = uiaxes(app.UIFigure);
            title(app.FrequencyResponseAxes, 'Frequency Response')
            xlabel(app.FrequencyResponseAxes, 'X')
            ylabel(app.FrequencyResponseAxes, 'Y')
            zlabel(app.FrequencyResponseAxes, 'Z')
            app.FrequencyResponseAxes.Position = [170 143 300 185];

            % Create ModeDropDownLabel
            app.ModeDropDownLabel = uilabel(app.UIFigure);
            app.ModeDropDownLabel.HorizontalAlignment = 'right';
            app.ModeDropDownLabel.Position = [66 575 35 22];
            app.ModeDropDownLabel.Text = 'Mode';

            % Create ModeDropDown
            app.ModeDropDown = uidropdown(app.UIFigure);
            app.ModeDropDown.Items = {'synthetic', 'recorded'};
            app.ModeDropDown.ValueChangedFcn = createCallbackFcn(app, @ModeDropdownValueChanged, true);
            app.ModeDropDown.Position = [116 575 100 22];
            app.ModeDropDown.Value = 'synthetic';

            % Create SignalTypeDropDownLabel
            app.SignalTypeDropDownLabel = uilabel(app.UIFigure);
            app.SignalTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SignalTypeDropDownLabel.Position = [33 532 68 22];
            app.SignalTypeDropDownLabel.Text = 'Signal Type';

            % Create SignalTypeDropDown
            app.SignalTypeDropDown = uidropdown(app.UIFigure);
            app.SignalTypeDropDown.Items = {'sine', 'square', 'triangle', 'chirp', 'saturated'};
            app.SignalTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SignalDropdownValueChanged, true);
            app.SignalTypeDropDown.Position = [116 532 100 22];
            app.SignalTypeDropDown.Value = 'sine';

            % Create DuracisEditFieldLabel
            app.DuracisEditFieldLabel = uilabel(app.UIFigure);
            app.DuracisEditFieldLabel.HorizontalAlignment = 'right';
            app.DuracisEditFieldLabel.Position = [259 575 64 22];
            app.DuracisEditFieldLabel.Text = 'Duració (s)';

            % Create DurationField
            app.DurationField = uieditfield(app.UIFigure, 'numeric');
            app.DurationField.ValueChangedFcn = createCallbackFcn(app, @DurationFieldValueChanged, true);
            app.DurationField.Position = [338 575 100 22];
            app.DurationField.Value = 1;

            % Create FField
            app.FField = uieditfield(app.UIFigure, 'numeric');
            app.FField.ValueChangedFcn = createCallbackFcn(app, @FFieldValueChanged, true);
            app.FField.Position = [338 532 100 22];
            app.FField.Value = 440;

            % Create FFieldLabel
            app.FFieldLabel = uilabel(app.UIFigure);
            app.FFieldLabel.Position = [264 532 59 22];
            app.FFieldLabel.Text = 'Freq. (Hz)';

            % Create RecordButton
            app.RecordButton = uibutton(app.UIFigure, 'push');
            app.RecordButton.ButtonPushedFcn = createCallbackFcn(app, @RecordButtonPushed, true);
            app.RecordButton.Visible = 'off';
            app.RecordButton.Position = [487 532 100 65];
            app.RecordButton.Text = 'Record';

            % Create PowerLabel
            app.PowerLabel = uilabel(app.UIFigure);
            app.PowerLabel.HorizontalAlignment = 'right';
            app.PowerLabel.Position = [41 114 39 22];
            app.PowerLabel.Text = 'Power';

            % Create PowerEditField
            app.PowerEditField = uieditfield(app.UIFigure, 'numeric');
            app.PowerEditField.Editable = 'off';
            app.PowerEditField.Position = [95 114 100 22];

            % Create GainLabel
            app.GainLabel = uilabel(app.UIFigure);
            app.GainLabel.HorizontalAlignment = 'right';
            app.GainLabel.Position = [248 114 30 22];
            app.GainLabel.Text = 'Gain';

            % Create GainEditField
            app.GainEditField = uieditfield(app.UIFigure, 'numeric');
            app.GainEditField.Editable = 'off';
            app.GainEditField.Position = [293 114 100 22];

            % Create THDLabel
            app.THDLabel = uilabel(app.UIFigure);
            app.THDLabel.HorizontalAlignment = 'right';
            app.THDLabel.Position = [442 114 30 22];
            app.THDLabel.Text = 'THD';

            % Create THDEditField
            app.THDEditField = uieditfield(app.UIFigure, 'numeric');
            app.THDEditField.Editable = 'off';
            app.THDEditField.Position = [487 114 100 22];

            % Create PlayButton
            app.PlayButton = uibutton(app.UIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.Position = [271 54 100 39];
            app.PlayButton.Text = 'Play';

            % Create SavewavButton
            app.SavewavButton = uibutton(app.UIFigure, 'push');
            app.SavewavButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SavewavButton.Position = [33 63 100 22];
            app.SavewavButton.Text = 'Save .wav';

            % Create ExportconfigmatButton
            app.ExportconfigmatButton = uibutton(app.UIFigure, 'push');
            app.ExportconfigmatButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.ExportconfigmatButton.Position = [386 63 111 22];
            app.ExportconfigmatButton.Text = 'Export config .mat';

            % Create ImportconfigmatButton
            app.ImportconfigmatButton = uibutton(app.UIFigure, 'push');
            app.ImportconfigmatButton.ButtonPushedFcn = createCallbackFcn(app, @ImportButtonPushed, true);
            app.ImportconfigmatButton.Position = [512 63 110 22];
            app.ImportconfigmatButton.Text = 'Import config .mat';

            % Create ImportwavButton
            app.ImportwavButton = uibutton(app.UIFigure, 'push');
            app.ImportwavButton.ButtonPushedFcn = createCallbackFcn(app, @ImportAudioButtonPushed, true);
            app.ImportwavButton.Position = [149 63 100 22];
            app.ImportwavButton.Text = 'Import .wav';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MScomentat

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