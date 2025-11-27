classdef MSDEFINITIU_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        ChangeAmpModeButton           matlab.ui.control.Button
        RefreshButton                 matlab.ui.control.Button
        ChooseamplificationmodePanel  matlab.ui.container.Panel
        btn                           matlab.ui.control.Button
        dd                            matlab.ui.control.DropDown
        ModeAmpliLabel                matlab.ui.control.Label
        ImportwavButton               matlab.ui.control.Button
        ImportconfigmatButton         matlab.ui.control.Button
        ExportconfigmatButton         matlab.ui.control.Button
        SavewavButton                 matlab.ui.control.Button
        PlayButton                    matlab.ui.control.Button
        THDEditField                  matlab.ui.control.NumericEditField
        THDLabel                      matlab.ui.control.Label
        GainEditField                 matlab.ui.control.NumericEditField
        GainLabel                     matlab.ui.control.Label
        PowerEditField                matlab.ui.control.NumericEditField
        PowerLabel                    matlab.ui.control.Label
        RecordButton                  matlab.ui.control.Button
        FFieldLabel                   matlab.ui.control.Label
        FField                        matlab.ui.control.NumericEditField
        DurationField                 matlab.ui.control.NumericEditField
        DuracisEditFieldLabel         matlab.ui.control.Label
        SignalTypeDropDown            matlab.ui.control.DropDown
        SignalTypeDropDownLabel       matlab.ui.control.Label
        ModeDropDown                  matlab.ui.control.DropDown
        ModeDropDownLabel             matlab.ui.control.Label
        AmplifiedSpectrogramAxes      matlab.ui.control.UIAxes
        FrequencyResponseAxes         matlab.ui.control.UIAxes
        SpectrogramAxes               matlab.ui.control.UIAxes
        WaveformAxes                  matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        FS = 44100;                                                         % Sampling frequency
        y                                                                   % Signal used for analysis/plotting
        AmplificationMode; 
        MainApp;
    end  
    
    methods (Access = private)

        function generateSignal(app)                                        % Create a synthetic signal based on type and frequency
            t = 0:1/app.FS:app.DurationField.Value - 1/app.FS;              % Signal duration
            F = app.FField.Value;                                           % Frequency of synthetic signal
            type = app.SignalTypeDropDown.Value;                            % Signal type selection
        
            switch type
                case 'sine'
                    app.y = sin(2*pi*F*t)';
                case 'square'
                    app.y = sign(sin(2*pi*F*t))';
                case 'triangle'
                     app.y = (2 * abs(2 * mod(F * t, 1) - 1) - 1)';                
                case 'chirp'
                    durationch     = 5;      
                    Nbits        = 16;      
                    Fch       = 20000;                              
                    samples = durationch*app.FS;
                    fmi     = (0:1:samples-1)'./(samples-1);      
                    phi       = pi*(Fch/app.FS)*(0:samples-1)';
                    app.y  = sin(phi.*fmi);                  
                case 'saturated'
                    app.y = max(min(5*sin(2*pi*F*t), 1), -1)';
            end
        end
    
        function analyzeAndPlot(app, amplifiedSignal)
            % Analyze signal: power, gain, THD
            y = app.y;
            Fs = app.FS;  % Sampling frequency
        
            % Power
            if isempty(amplifiedSignal)
                P = 0;
                a = 'amplified signal esta empty'
            else
                P = sum(abs(amplifiedSignal.^2)) / length(amplifiedSignal);
            end
            app.PowerEditField.Value = P;
            
            % Determine gain calculation based on mode
            if strcmp(app.AmplificationMode, 'none')
                gainAmplified = 1;
            else
                gainAmplified = 10 * log10(var(amplifiedSignal) / var(app.y));
            end
            
            app.GainEditField.Value = gainAmplified; % Update gain in UI
        
            % THD
            AY = fft(amplifiedSignal);                                       % Compute FFT of the signal
            N = length(AY);                                    % Length of the FFT
            f = (0:N-1)*(app.FS/N);                                         % Frequency vector
            mag = abs(AY(1:floor(N/2)));                       % Magnitude spectrum (only positive freqs)
            f_half = f(1:floor(N/2));                                       % Corresponding frequencies
        
            [~, idx_fund] = max(mag);                                       % Find fundamental frequency index
            f_fund = f_half(idx_fund);                                      % Fundamental frequency value
            P1 = mag(idx_fund).^2;                                           % Power of the fundamental component
        
            harmonics = zeros(1, 4);
            for k = 2:5                                                     % Loop over 2nd to 5th harmonics
                target_freq = f_fund * k;
                [~, idx] = min(abs(f_half - target_freq));                  % Find closest freq index
                idx = idx(1);
                harmonics(k-1) = mag(idx).^2;                                % Power of each harmonic
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
            THD_norm = min(max(THD_norm, 0), 1);
            color_signal = (1 - THD_norm)*[0 0 1] + THD_norm*[1 0 0];       % Interpolate color (blue to red)
        
            % Waveform creation
            L = min(length(app.y), length(amplifiedSignal));
            xlabel(app.WaveformAxes, 'Temps (s)');
            ylabel(app.WaveformAxes, 'Amplitud');
            plot(app.WaveformAxes, (1:L)/Fs, app.y(1:L), 'Color', color_signal, 'LineWidth', 1.5);
            hold(app.WaveformAxes, 'on');
            plot(app.WaveformAxes, (1:L)/Fs, amplifiedSignal(1:L), 'Color', color_signal, 'LineWidth', 1.5);
            hold(app.WaveformAxes, 'off');
            title(app.WaveformAxes, '\color[rgb]{0.94,0.94,0.94}Waveform Comparison');
                    
            % Manual input spectrogram
            NFFT = 1024;                                                    % FFT size
            overlap = 0.25 * NFFT;                                          % 25% overlap between frames
            hop = NFFT - overlap;                                           % Hop size
            win = 0.54 - 0.46 * cos(2 * pi * (0:NFFT-1)' / (NFFT - 1));     % Hamming window
            win = win / sqrt(sum(win.^2));                                  % Normalize window energy
            
            % If signal is too short, clear axes and show warning
            if length(app.y) < NFFT
                cla(app.SpectrogramAxes);
                title(app.SpectrogramAxes, 'Signal too short');
                return;
            end
            if length(amplifiedSignal) < NFFT
                cla(app.AmplifiedSpectrogramAxes);
                title(app.AmplifiedSpectrogramAxes, 'Signal too short');
                return;
            end
        
            numFrames = floor((length(app.y) - NFFT) / hop) + 1;                % Number of frames
            S = zeros(NFFT/2, numFrames);                                   % Initialize spectrogram matrix
            T = zeros(1, numFrames);                                        % Time vector for each frame

            % Compute FFT for each frame
            app.y = app.y(:);
            for i = 1:numFrames
                startIdx = (i-1)*hop + 1;
                endIdx = startIdx + NFFT - 1;
                frame = app.y(startIdx:endIdx, 1) .* win;                   % Apply window to frame
                Y_frame = fft(frame, NFFT);                                 % Compute FFT
                S(:,i) = abs(Y_frame(1:NFFT/2));                            % Keep positive frequencies
                T(i) = startIdx / app.FS;                                   % Time in seconds
            end
            %% Espectrograma de la senyal amplificada amb filtre de silenci
            AnumFrames = floor((length(amplifiedSignal) - NFFT) / hop) + 1;
            AS = zeros(NFFT/2, AnumFrames);
            AT = zeros(1, AnumFrames);

            % Defineix el llindar d'energia per considerar un frame com a silenci
            amplifiedSignal = amplifiedSignal(:);
            % Define an energy threshold for silence (tweak this value as needed)
            energyThreshold = 1e-4;        % Example threshold value
            nonSilentIndices = [];         % To record indices of non-silent frames

            for i = 1:AnumFrames
                AstartIdx = (i-1)*hop + 1;
                AendIdx   = AstartIdx + NFFT - 1;
                Aframe    = amplifiedSignal(AstartIdx:AendIdx, 1) .*win;   % Apply window
                frameEnergy = sum(Aframe.^2);                               % Compute frame energy

                if frameEnergy >= energyThreshold
                    nonSilentIndices = [nonSilentIndices, i];             % Record non-silent frame index
                    AY_frame = fft(Aframe, NFFT);                           % Compute FFT
                    AS(:, i) = abs(AY_frame(1:NFFT/2));                     % Store magnitude of positive frequencies
                    AT(i) = AstartIdx / app.FS;
                else
                    % For silent frames, leave the column as zeros.
                    AS(:, i) = zeros(NFFT/2, 1);
                    AT(i) = AstartIdx / app.FS;
                end
            end

                % Remove silent frames from the spectrogram data
                AS_filtered = AS(:, nonSilentIndices);
                AT_filtered = AT(nonSilentIndices);

            f_vec = Fs * (0:(NFFT/2-1)) / NFFT;                            % Frequency vector

            % Original Spectrogram
            surf(app.SpectrogramAxes, T, f_vec, 20*log10(S+eps), 'EdgeColor', 'none');
            axis(app.SpectrogramAxes, 'tight');
            view(app.SpectrogramAxes, [0 90]);
            xlabel(app.SpectrogramAxes, 'Time (s)');
            ylabel(app.SpectrogramAxes, 'Frequency (Hz)');
            title(app.SpectrogramAxes, '\color[rgb]{0.94,0.94,0.94}Spectrogram - Original Signal');
            colorbar(app.SpectrogramAxes);
        
            % Amplified Signal Spectrogram with Silence Removed
            surf(app.AmplifiedSpectrogramAxes, AT_filtered, f_vec, 20*log10(AS_filtered+eps), 'EdgeColor', 'none');
            axis(app.AmplifiedSpectrogramAxes, 'tight');
            view(app.AmplifiedSpectrogramAxes, [0 90]);
            xlabel(app.AmplifiedSpectrogramAxes, 'Time (s)');
            ylabel(app.AmplifiedSpectrogramAxes, 'Frequency (Hz)');
            title(app.AmplifiedSpectrogramAxes, '\color[rgb]{0.94,0.94,0.94}Spectrogram - Amplified Signal');
            colorbar(app.AmplifiedSpectrogramAxes);

            % Frequency Response
            app.y = app.y(:);

            N = length(app.y);
            n = (0:N-1)';
            w_hann = 0.5*(1-cos(2*pi*n/(N-1)));
            size(w_hann)
            Y_orig = fft(app.y.*w_hann);
            Y_amp = fft(amplifiedSignal.*w_hann); 
            f_fresp = (0:N-1)*(Fs/N);
            half_N = floor(N/2);
            f_fresp = f_fresp(1:half_N);
            Y_orig = Y_orig(1:half_N);
            Y_amp = Y_amp(1:half_N);
            Y_orig_mag = abs(Y_orig)/max(abs(Y_orig));
            Y_amp_mag = abs(Y_amp)/max(abs(Y_amp));
            Y_orig_dB = 20*log10(Y_orig_mag+eps);
            Y_amp_dB = 20*log10(Y_amp_mag+eps);                                   
            
            plot(app.FrequencyResponseAxes, f_fresp, Y_orig_dB, 'b', 'LineWidth', 1.5);
            hold(app.FrequencyResponseAxes, 'on');
            plot(app.FrequencyResponseAxes, f_fresp, Y_amp_dB, 'r', 'LineWidth', 1.5);
            hold(app.FrequencyResponseAxes, 'off');
            title(app.FrequencyResponseAxes, '\color[rgb]{0.94,0.94,0.94}Frequency Response Comparison');
            xlabel(app.FrequencyResponseAxes, 'Frequency (Hz)');
            ylabel(app.FrequencyResponseAxes, 'Magnitude (dB)');
            grid(app.FrequencyResponseAxes, 'on');
            xlim(app.FrequencyResponseAxes, [0 Fs/2]);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, MainApp)
            app.MainApp = MainApp;
            app.generateSignal;
        end

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
            generateSignal(app);                                            % Generate signal based on selection
        end

        % Value changed function: DurationField
        function DurationFieldValueChanged(app, event)
            val = app.DurationField.Value;
            minVal = 0.1;
            maxVal = 10;
            if val < minVal
                val = minVal;
            elseif val > maxVal
                val = maxVal;
            end
            app.DurationField.Value = val;
            generateSignal(app);
        end

        % Value changed function: FField
        function FFieldValueChanged(app, event)
            val = app.FField.Value;
            minVal = 20;
            maxVal = 20000;
            if val < minVal
                val = minVal;
            elseif val > maxVal
                val = maxVal;
            end
            app.FField.Value = val;
            generateSignal(app);
        end

        % Button pushed function: RecordButton
        function RecordButtonPushed(app, event)
            duration = app.DurationField.Value;
            countdownLabel = uilabel(app.UIFigure);
            countdownLabel.Position = [487 510 120 22];
            countdownLabel.HorizontalAlignment = 'center';
            countdownLabel.FontSize = 14;
            countdownLabel.Text = 'Iniciant gravació en...';
            for i = 3:-1:1
                countdownLabel.Text = sprintf('Gravant en %d...', i);
                pause(1);
            end
            countdownLabel.Text = 'Gravant...';
            rec = audiorecorder(app.FS, 24, 1);
            record(rec);
            for i = 1:duration+1
                countdownLabel.Text = sprintf('Temps restant: %d s', duration-i+1);
                pause(1);
            end
            stop(rec);
            countdownLabel.Text = 'Gravació completada';
            pause(0.5);
            delete(countdownLabel);
            y = getaudiodata(rec);
            threshold = 0.02;
            idx = find(abs(y) > threshold, 1, 'first');
            if ~isempty(idx)
                y = y(idx:end);
            end
            app.y = y(1:duration*app.FS);
            % En gravació, utilitzem la senyal importada com original
            analyzeAndPlot(app, app.y);
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
            analyzeAndPlot(app, app.y);
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

        % Callback function
        function AmplificationModeChanged(app, event)

        end

        % Button pushed function: btn
        function ButtonPushedFcn(app, event)
          
            selectedMode = app.dd.Value;
            
            
            app.AmplificationMode = selectedMode; 
            
       
            app.ChooseamplificationmodePanel.Visible = 'off';
        end

        % Button pushed function: RefreshButton
        function RefreshButtonPushed(app, event)
            mode = app.AmplificationMode;
            switch mode
                case 'none'
                    amplifiedSignal = app.y;
                case 'simulated'
                    % Rebem la senyal amplificada d'una altra app (el simulador)
                    amplifiedSignal = app.MainApp.Ampli.y;
                case 'external'
                    sound(app.y, app.FS)
                    rec = audiorecorder(app.FS, 24, 1);
                    durationRecord = max(length(app.y)/app.FS, 1);
                    recordblocking(rec, durationRecord);
                    amplifiedSignal = getaudiodata(rec);
            end
            analyzeAndPlot(app, amplifiedSignal);
        end

        % Button pushed function: ChangeAmpModeButton
        function ChangeAmpModeButtonPushed(app, event)
            app.ChooseamplificationmodePanel.Visible = 'on';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.149 0.149 0.149];
            app.UIFigure.Position = [100 100 638 635];
            app.UIFigure.Name = 'MATLAB App';

            % Create WaveformAxes
            app.WaveformAxes = uiaxes(app.UIFigure);
            title(app.WaveformAxes, '\color[rgb]{0.94,0.94,0.94}Waveform')
            xlabel(app.WaveformAxes, 'X')
            ylabel(app.WaveformAxes, 'Y')
            zlabel(app.WaveformAxes, 'Z')
            app.WaveformAxes.XColor = [0.9412 0.9412 0.9412];
            app.WaveformAxes.YColor = [0.9412 0.9412 0.9412];
            app.WaveformAxes.ZColor = [0.9412 0.9412 0.9412];
            app.WaveformAxes.Color = [0.149 0.149 0.149];
            app.WaveformAxes.Position = [1 343 300 185];

            % Create SpectrogramAxes
            app.SpectrogramAxes = uiaxes(app.UIFigure);
            title(app.SpectrogramAxes, '\color[rgb]{0.94,0.94,0.94}Spectrogram')
            xlabel(app.SpectrogramAxes, 'X')
            ylabel(app.SpectrogramAxes, 'Y')
            zlabel(app.SpectrogramAxes, 'Z')
            app.SpectrogramAxes.XColor = [0.9412 0.9412 0.9412];
            app.SpectrogramAxes.YColor = [0.9412 0.9412 0.9412];
            app.SpectrogramAxes.ZColor = [0.9412 0.9412 0.9412];
            app.SpectrogramAxes.Color = [0.149 0.149 0.149];
            colormap(app.SpectrogramAxes, 'jet')
            app.SpectrogramAxes.Position = [322 347 300 185];

            % Create FrequencyResponseAxes
            app.FrequencyResponseAxes = uiaxes(app.UIFigure);
            title(app.FrequencyResponseAxes, '\color[rgb]{0.94,0.94,0.94}Frequency Response')
            xlabel(app.FrequencyResponseAxes, 'X')
            ylabel(app.FrequencyResponseAxes, 'Y')
            zlabel(app.FrequencyResponseAxes, 'Z')
            app.FrequencyResponseAxes.XColor = [0.9412 0.9412 0.9412];
            app.FrequencyResponseAxes.YColor = [0.9412 0.9412 0.9412];
            app.FrequencyResponseAxes.ZColor = [0.9412 0.9412 0.9412];
            app.FrequencyResponseAxes.Color = [0.149 0.149 0.149];
            app.FrequencyResponseAxes.Position = [2 140 300 185];

            % Create AmplifiedSpectrogramAxes
            app.AmplifiedSpectrogramAxes = uiaxes(app.UIFigure);
            title(app.AmplifiedSpectrogramAxes, '\color[rgb]{0.94,0.94,0.94}Title')
            xlabel(app.AmplifiedSpectrogramAxes, 'X')
            ylabel(app.AmplifiedSpectrogramAxes, 'Y')
            zlabel(app.AmplifiedSpectrogramAxes, 'Z')
            app.AmplifiedSpectrogramAxes.XColor = [0.9412 0.9412 0.9412];
            app.AmplifiedSpectrogramAxes.YColor = [0.9412 0.9412 0.9412];
            app.AmplifiedSpectrogramAxes.ZColor = [0.9412 0.9412 0.9412];
            app.AmplifiedSpectrogramAxes.Color = [0.149 0.149 0.149];
            colormap(app.AmplifiedSpectrogramAxes, 'jet')
            app.AmplifiedSpectrogramAxes.Position = [323 140 300 185];

            % Create ModeDropDownLabel
            app.ModeDropDownLabel = uilabel(app.UIFigure);
            app.ModeDropDownLabel.BackgroundColor = [0.149 0.149 0.149];
            app.ModeDropDownLabel.HorizontalAlignment = 'right';
            app.ModeDropDownLabel.FontColor = [0.9412 0.9412 0.9412];
            app.ModeDropDownLabel.Position = [85 602 35 22];
            app.ModeDropDownLabel.Text = 'Mode';

            % Create ModeDropDown
            app.ModeDropDown = uidropdown(app.UIFigure);
            app.ModeDropDown.Items = {'synthetic', 'recorded'};
            app.ModeDropDown.ValueChangedFcn = createCallbackFcn(app, @ModeDropdownValueChanged, true);
            app.ModeDropDown.FontColor = [0.9412 0.9412 0.9412];
            app.ModeDropDown.BackgroundColor = [0.149 0.149 0.149];
            app.ModeDropDown.Position = [135 602 100 22];
            app.ModeDropDown.Value = 'synthetic';

            % Create SignalTypeDropDownLabel
            app.SignalTypeDropDownLabel = uilabel(app.UIFigure);
            app.SignalTypeDropDownLabel.BackgroundColor = [0.149 0.149 0.149];
            app.SignalTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SignalTypeDropDownLabel.FontColor = [0.9412 0.9412 0.9412];
            app.SignalTypeDropDownLabel.Position = [52 565 68 22];
            app.SignalTypeDropDownLabel.Text = 'Signal Type';

            % Create SignalTypeDropDown
            app.SignalTypeDropDown = uidropdown(app.UIFigure);
            app.SignalTypeDropDown.Items = {'sine', 'square', 'triangle', 'chirp', 'saturated'};
            app.SignalTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SignalDropdownValueChanged, true);
            app.SignalTypeDropDown.FontColor = [0.9412 0.9412 0.9412];
            app.SignalTypeDropDown.BackgroundColor = [0.149 0.149 0.149];
            app.SignalTypeDropDown.Position = [135 565 100 22];
            app.SignalTypeDropDown.Value = 'sine';

            % Create DuracisEditFieldLabel
            app.DuracisEditFieldLabel = uilabel(app.UIFigure);
            app.DuracisEditFieldLabel.BackgroundColor = [0.149 0.149 0.149];
            app.DuracisEditFieldLabel.HorizontalAlignment = 'right';
            app.DuracisEditFieldLabel.FontColor = [0.9412 0.9412 0.9412];
            app.DuracisEditFieldLabel.Position = [259 595 64 22];
            app.DuracisEditFieldLabel.Text = 'Duració (s)';

            % Create DurationField
            app.DurationField = uieditfield(app.UIFigure, 'numeric');
            app.DurationField.ValueChangedFcn = createCallbackFcn(app, @DurationFieldValueChanged, true);
            app.DurationField.FontColor = [0.9412 0.9412 0.9412];
            app.DurationField.BackgroundColor = [0.149 0.149 0.149];
            app.DurationField.Position = [338 595 100 22];
            app.DurationField.Value = 1;

            % Create FField
            app.FField = uieditfield(app.UIFigure, 'numeric');
            app.FField.ValueChangedFcn = createCallbackFcn(app, @FFieldValueChanged, true);
            app.FField.FontColor = [0.9412 0.9412 0.9412];
            app.FField.BackgroundColor = [0.149 0.149 0.149];
            app.FField.Position = [338 552 100 22];
            app.FField.Value = 440;

            % Create FFieldLabel
            app.FFieldLabel = uilabel(app.UIFigure);
            app.FFieldLabel.FontColor = [0.9412 0.9412 0.9412];
            app.FFieldLabel.Position = [264 552 59 22];
            app.FFieldLabel.Text = 'Freq. (Hz)';

            % Create RecordButton
            app.RecordButton = uibutton(app.UIFigure, 'push');
            app.RecordButton.ButtonPushedFcn = createCallbackFcn(app, @RecordButtonPushed, true);
            app.RecordButton.BackgroundColor = [0.149 0.149 0.149];
            app.RecordButton.FontColor = [0.9412 0.9412 0.9412];
            app.RecordButton.Visible = 'off';
            app.RecordButton.Position = [462 552 73 65];
            app.RecordButton.Text = 'Record';

            % Create PowerLabel
            app.PowerLabel = uilabel(app.UIFigure);
            app.PowerLabel.BackgroundColor = [0.149 0.149 0.149];
            app.PowerLabel.HorizontalAlignment = 'right';
            app.PowerLabel.FontColor = [0.9412 0.9412 0.9412];
            app.PowerLabel.Position = [41 101 39 22];
            app.PowerLabel.Text = 'Power';

            % Create PowerEditField
            app.PowerEditField = uieditfield(app.UIFigure, 'numeric');
            app.PowerEditField.Editable = 'off';
            app.PowerEditField.FontColor = [0.9412 0.9412 0.9412];
            app.PowerEditField.BackgroundColor = [0.149 0.149 0.149];
            app.PowerEditField.Position = [95 101 100 22];

            % Create GainLabel
            app.GainLabel = uilabel(app.UIFigure);
            app.GainLabel.BackgroundColor = [0.149 0.149 0.149];
            app.GainLabel.HorizontalAlignment = 'right';
            app.GainLabel.FontColor = [0.9412 0.9412 0.9412];
            app.GainLabel.Position = [248 101 30 22];
            app.GainLabel.Text = 'Gain';

            % Create GainEditField
            app.GainEditField = uieditfield(app.UIFigure, 'numeric');
            app.GainEditField.Editable = 'off';
            app.GainEditField.FontColor = [0.9412 0.9412 0.9412];
            app.GainEditField.BackgroundColor = [0.149 0.149 0.149];
            app.GainEditField.Position = [293 101 100 22];

            % Create THDLabel
            app.THDLabel = uilabel(app.UIFigure);
            app.THDLabel.BackgroundColor = [0.149 0.149 0.149];
            app.THDLabel.HorizontalAlignment = 'right';
            app.THDLabel.FontColor = [0.9412 0.9412 0.9412];
            app.THDLabel.Position = [442 101 30 22];
            app.THDLabel.Text = 'THD';

            % Create THDEditField
            app.THDEditField = uieditfield(app.UIFigure, 'numeric');
            app.THDEditField.Editable = 'off';
            app.THDEditField.FontColor = [0.9412 0.9412 0.9412];
            app.THDEditField.BackgroundColor = [0.149 0.149 0.149];
            app.THDEditField.Position = [487 101 100 22];

            % Create PlayButton
            app.PlayButton = uibutton(app.UIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.BackgroundColor = [0.149 0.149 0.149];
            app.PlayButton.FontColor = [0.9412 0.9412 0.9412];
            app.PlayButton.Position = [271 41 100 39];
            app.PlayButton.Text = 'Play';

            % Create SavewavButton
            app.SavewavButton = uibutton(app.UIFigure, 'push');
            app.SavewavButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SavewavButton.BackgroundColor = [0.149 0.149 0.149];
            app.SavewavButton.FontColor = [0.9412 0.9412 0.9412];
            app.SavewavButton.Position = [33 50 100 22];
            app.SavewavButton.Text = 'Save .wav';

            % Create ExportconfigmatButton
            app.ExportconfigmatButton = uibutton(app.UIFigure, 'push');
            app.ExportconfigmatButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.ExportconfigmatButton.BackgroundColor = [0.149 0.149 0.149];
            app.ExportconfigmatButton.FontColor = [0.9412 0.9412 0.9412];
            app.ExportconfigmatButton.Position = [386 50 111 22];
            app.ExportconfigmatButton.Text = 'Export config .mat';

            % Create ImportconfigmatButton
            app.ImportconfigmatButton = uibutton(app.UIFigure, 'push');
            app.ImportconfigmatButton.ButtonPushedFcn = createCallbackFcn(app, @ImportButtonPushed, true);
            app.ImportconfigmatButton.BackgroundColor = [0.149 0.149 0.149];
            app.ImportconfigmatButton.FontColor = [0.9412 0.9412 0.9412];
            app.ImportconfigmatButton.Position = [512 50 110 22];
            app.ImportconfigmatButton.Text = 'Import config .mat';

            % Create ImportwavButton
            app.ImportwavButton = uibutton(app.UIFigure, 'push');
            app.ImportwavButton.ButtonPushedFcn = createCallbackFcn(app, @ImportAudioButtonPushed, true);
            app.ImportwavButton.BackgroundColor = [0.149 0.149 0.149];
            app.ImportwavButton.FontColor = [0.9412 0.9412 0.9412];
            app.ImportwavButton.Position = [149 50 100 22];
            app.ImportwavButton.Text = 'Import .wav';

            % Create ChooseamplificationmodePanel
            app.ChooseamplificationmodePanel = uipanel(app.UIFigure);
            app.ChooseamplificationmodePanel.ForegroundColor = [1 1 1];
            app.ChooseamplificationmodePanel.Title = 'Choose amplification mode';
            app.ChooseamplificationmodePanel.BackgroundColor = [0.149 0.149 0.149];
            app.ChooseamplificationmodePanel.Position = [194 303 260 137];

            % Create ModeAmpliLabel
            app.ModeAmpliLabel = uilabel(app.ChooseamplificationmodePanel);
            app.ModeAmpliLabel.BackgroundColor = [0.149 0.149 0.149];
            app.ModeAmpliLabel.HorizontalAlignment = 'right';
            app.ModeAmpliLabel.FontColor = [1 1 1];
            app.ModeAmpliLabel.Position = [33 67 68 22];
            app.ModeAmpliLabel.Text = 'Mode Ampli';

            % Create dd
            app.dd = uidropdown(app.ChooseamplificationmodePanel);
            app.dd.Items = {'none', 'simulated', 'external'};
            app.dd.FontColor = [1 1 1];
            app.dd.BackgroundColor = [0.149 0.149 0.149];
            app.dd.Position = [116 67 100 22];
            app.dd.Value = 'none';

            % Create btn
            app.btn = uibutton(app.ChooseamplificationmodePanel, 'push');
            app.btn.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushedFcn, true);
            app.btn.BackgroundColor = [0.149 0.149 0.149];
            app.btn.FontColor = [1 1 1];
            app.btn.Position = [77 23 100 22];
            app.btn.Text = 'confirm';

            % Create RefreshButton
            app.RefreshButton = uibutton(app.UIFigure, 'push');
            app.RefreshButton.ButtonPushedFcn = createCallbackFcn(app, @RefreshButtonPushed, true);
            app.RefreshButton.BackgroundColor = [0.149 0.149 0.149];
            app.RefreshButton.FontColor = [0.9412 0.9412 0.9412];
            app.RefreshButton.Position = [549 552 73 65];
            app.RefreshButton.Text = 'Refresh';

            % Create ChangeAmpModeButton
            app.ChangeAmpModeButton = uibutton(app.UIFigure, 'push');
            app.ChangeAmpModeButton.ButtonPushedFcn = createCallbackFcn(app, @ChangeAmpModeButtonPushed, true);
            app.ChangeAmpModeButton.BackgroundColor = [0.149 0.149 0.149];
            app.ChangeAmpModeButton.FontColor = [0.9412 0.9412 0.9412];
            app.ChangeAmpModeButton.Position = [271 11 100 22];
            app.ChangeAmpModeButton.Text = 'Amplificació';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MSDEFINITIU_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

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