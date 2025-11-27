classdef Eqcodi < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        Toolbar                       matlab.ui.container.Toolbar
        ToggleTool                    matlab.ui.container.toolbar.ToggleTool
        addButton                     matlab.ui.container.toolbar.PushTool
        editButton                    matlab.ui.container.toolbar.PushTool
        deleteButton                  matlab.ui.container.toolbar.PushTool
        saveButton                    matlab.ui.container.toolbar.PushTool
        exportButton                  matlab.ui.container.toolbar.PushTool
        Export                        matlab.ui.container.Panel
        FormatDropDown                matlab.ui.control.DropDown
        TipusLabel                    matlab.ui.control.Label
        XButton_2                     matlab.ui.control.Button
        ExporttoButton                matlab.ui.control.Button
        FileNameEditField             matlab.ui.control.EditField
        FileNameEditFieldLabel        matlab.ui.control.Label
        PauseButton                   matlab.ui.control.Button
        PlayButton                    matlab.ui.control.Button
        Button                        matlab.ui.control.Button
        AudioSlider                   matlab.ui.control.Slider
        FiltreDropDown                matlab.ui.control.DropDown
        colorrgb094094094FiltreLabel  matlab.ui.control.Label
        AfegirfiltrePanel             matlab.ui.container.Panel
        XButton                       matlab.ui.control.Button
        AfegirButton                  matlab.ui.control.Button
        TipusdefiltreDropDown         matlab.ui.control.DropDown
        TipusdefiltreDropDownLabel    matlab.ui.control.Label
        UIAxes                        matlab.ui.control.UIAxes
        UIAxes_2                      matlab.ui.control.UIAxes
        ContextMenu                   matlab.ui.container.ContextMenu
        AfegirfiltreMenu              matlab.ui.container.Menu
    end

    
    properties (Access = private)
        %Eq
        fs = 48000;
        nextID = 1;

        audio;
        audioName;

        audioSegsL;
        audioSegsR;
        audioEQL;
        audioEQR;
        audioWindow = 1;
        audioFiltrat;

        filtersOn = 1;

        ExportMode = 'audio'; % 'audio' o 'param'


        N = 4096;
        H;
        Panell; % app gestó filtre 
    end
    
    properties (Access = public)
        filtres;
        filtresParam;
        player; 
        isChanging = 0;
    end
    
    methods (Access = public)
        
        function actualitza_eq(app)
            
            [s, l]= size(app.filtres);
            f = logspace(log10(10), log10(app.fs/2), app.N);
            app.H=ones(1,app.N);
            H1=ones(1,app.N);
            if (l == app.N)
                for i=1:s
                    H1=H1 .* app.filtres(i,:);
                    area(app.UIAxes, f, 20*log10(abs(app.filtres(i,:))),'FaceAlpha',.2,'EdgeAlpha',.2);
                    hold(app.UIAxes,'on');
                end
            end
            app.H = app.H .*H1;
            semilogx(app.UIAxes, f, 20*log10(abs(app.H)), "w");
            hold(app.UIAxes,'off');

            if(~isplaying(app.player))
                app.actualitza_ft;
                pause(0.05);
            else
                if app.isChanging == 0
                    pause(app.player);
                    app.audioFiltrat = app.filtrar(app.audio(:,1), app.H, app.N*2, app.N);
                    if size(app.audio, 2)==2
                        app.audioFiltrat = [app.audioFiltrat; app.filtrar(app.audio(:,2), app.H, app.N*2, app.N)];
                    end
                    app.player = audioplayer(app.audioFiltrat, app.fs);
                    a = app.audioWindow * app.N/2 + app.N/2;
                    play(app.player,a);
                end
            end
        end

        function actualitza_ft(app)
            f = logspace(log10(10), log10(app.fs/2), app.N);
            
            bloque_filtrado_fft = app.audioEQL(round(app.audioWindow), :).* app.H;
            a = area(app.UIAxes_2, f, 20*log10(abs(bloque_filtrado_fft)*2/app.fs));
            a.EdgeColor = [0.5 0.5 0.5];
            a.FaceColor = [0.15,0.15,0.15];
            a.FaceAlpha = 0.75;
            a.LineWidth = 1.5;
            hold(app.UIAxes_2,'on');

            if size(app.audio, 2)==2
                bloque_filtrado_fft = app.audioEQR(round(app.audioWindow), :).*app.H;
                a = area(app.UIAxes_2, f, 20*log10(abs(bloque_filtrado_fft)*2/app.fs));
                a.EdgeColor = [0.5 0.5 0.5];
                a.FaceColor = [0.15,0.15,0.15];
                a.FaceAlpha = 0.75;
                a.LineWidth = 1.5;
            end
            hold(app.UIAxes_2,'off');
        end

    end
    
    methods (Access = private)
        
        % Enfinestrament temporal
        function w = hamming_custom(~,L)
            if L == 1
                w = 1;
            else
                n = (0:L-1)';
                w = 0.54 - 0.46*cos(2*pi*n/(L-1));
            end
        end
        
        function [bloques, ventanas] = segmentar_hamming(app,s, block_size, overlap)
            
            s = s(:).'; % Assegurar vector fila
            step = block_size - overlap;
            M = length(s);
        
            n_blocks = ceil((M - overlap) / step);
        
            % Zero padding
            padded_len = (n_blocks - 1) * step + block_size;
            s_padded = [s, zeros(1, padded_len - M)];
        
            % Inicializar la matriz de bloques y ventanas
            bloques = zeros(n_blocks, block_size);
            ventanas = zeros(n_blocks, block_size);
        
            w = app.hamming_custom(block_size).'; % Ventana de Hamming como fila
        
            for k = 0:n_blocks-1
                idx_start = k * step + 1;
                idx_end = idx_start + block_size - 1;
                bloque = s_padded(idx_start:idx_end);
                bloques(k+1, :) = bloque .* w;
                ventanas(k+1, :) = w;
            end
        end
        
        function s_rec = ajuntar_hamming(~, bloques, ventanas, block_size, overlap, original_length)
            
            step = block_size - overlap;
            n_blocks = size(bloques, 1);
            total_len = (n_blocks-1)*step + block_size;
            s_rec = zeros(1, total_len);
            norm_factor = zeros(1, total_len);
        
            for k = 0:n_blocks-1
                idx_start = k*step + 1;
                idx_end = idx_start + block_size - 1;
                s_rec(idx_start:idx_end) = s_rec(idx_start:idx_end) + bloques(k+1, :);
                norm_factor(idx_start:idx_end) = norm_factor(idx_start:idx_end) + ventanas(k+1, :);
            end
        
            % Normaliza para corregir la suma en regiones solapadas
            s_rec = s_rec ./ (norm_factor + eps);
        
            % Si se especifica la longitud original, recorta
            if nargin == 5
                s_rec = s_rec(1:original_length);
            end
        end
        
        function s_filtrada = filtrar_blocs(app, s, H, block_size, overlap)
            
        
            s = s(:).'; % Asegura vector fila
            Nfft = length(H);
        
            % Segmenta en bloques solapados y aplica ventana de Hamming
            [bloques, ~] = app.segmentar_hamming(s, block_size, overlap);
            n_blocks = size(bloques, 1);
        
            % Inicializa matriz para los bloques filtrados
            bloques_filtrados = zeros(size(bloques));
        
        
            f = linspace(10, app.fs/2, Nfft);
            f_l = logspace(log10(10), log10(app.fs/2), Nfft);
        
            % Filtra cada bloque en el dominio de la frecuencia
            for k = 1:n_blocks
                bloque = bloques(k, :);
                % FFT con zero-padding si block_size < Nfft
                bloque_fft = fft(bloque, 2*Nfft);
                bloque_log = interp1(f(1:end), bloque_fft(1:block_size), f_l, 'spline', 'extrap');

                % Solo la parte útil (del tamaño del bloque)
                bloques_filtrados(k, :) = bloque_log(1:block_size);
        
            end
        
            % Reconstruye la señal filtrada usando overlap-add y normalización
            s_filtrada = bloques_filtrados;
            
        end

        function s_filtrada = filtrar(app, s, H, block_size, overlap)
               
                s = s(:).'; % Asegura vector fila
                Nfft = length(H);
                [bloques, ventanas] = app.segmentar_hamming(s, block_size, overlap);
                n_blocks = size(bloques, 1);
            
                % Inicializa matriz para los bloques filtrados
                bloques_filtrados = zeros(size(bloques));
                f = linspace(10, app.fs/2, Nfft);
                f_l = logspace(log10(10), log10(app.fs/2), Nfft);

                H1 = interp1(f_l, H, f, 'spline', 0);
            
            for k = 1:n_blocks
                bloque = bloques(k, :);
                bloque_fft = fft(bloque, Nfft*2);
                bloque_filtrado_fft = bloque_fft .* [H1, flip(H1)];
                bloque_filtrado = real(ifft(bloque_filtrado_fft, Nfft*2)); % !!!!
                bloques_filtrados(k, :) = bloque_filtrado;
            end
            
                % Reconstruye la señal filtrada usando overlap-add y normalización
                s_filtrada = app.ajuntar_hamming(bloques_filtrados, ventanas, block_size, overlap, length(s));
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, audio, save)
            
            app.audioName = audio;
            [app.audio, app.fs] = audioread(audio);

            app.player = audioplayer(app.audio, app.fs);

            app.H=ones(1,app.N);

            app.audioSegsL = app.filtrar_blocs(app.audio(:, 1), app.H, app.N, app.N/2);
            
            app.audioEQL = app.audioSegsL;
            app.audioFiltrat = app.audio;

            if size(app.audio, 2)==2
                app.audioSegsR = app.filtrar_blocs(app.audio(:, 2), app.H, app.N, app.N/2);
                app.audioEQR = app.audioSegsR;
            end

            app.AudioSlider.Limits = [1 size(app.audioSegsL, 1)];
            if nargin >2
                p = load(save, 'param');
                param = p.param;
                for i = 1:size(param, 1)
                    if param(i,6)==1
                        app.TipusdefiltreDropDown.Value='Bell';
                    elseif param(i,6)==2
                        app.TipusdefiltreDropDown.Value='Shelf';
                    elseif param(i,6)==3
                        app.TipusdefiltreDropDown.Value='Band';
                    elseif param(i,6)==4
                        app.TipusdefiltreDropDown.Value='HighCut';
                    elseif param(i,6)==5
                        app.TipusdefiltreDropDown.Value='LowCut';
                    end
                    app.AfegirButtonPushed;
                    delete(app.Panell);
                    app.filtresParam(i, :) = param(i, :);
                    app.FiltreDropDownValueChanged;
                end
            end

            app.actualitza_eq;

        end

        % Value changed function: FiltreDropDown
        function FiltreDropDownValueChanged(app, event)
            ID = app.FiltreDropDown.Value;
            if ~isempty(app.Panell)
                if isvalid(app.Panell)
                    app.Panell.valors;
                    delete(app.Panell);
                end
            end
            if ID~=0
                if app.filtresParam(ID,6)==1
                    app.Panell = PanellBell(app, app.filtresParam(ID,1), app.filtresParam(ID,2), double(app.filtresParam(ID,3)), app.filtresParam(ID,4), app.fs, app.N);
                elseif app.filtresParam(ID,6)==2
                    app.Panell = PanellShelf(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.filtresParam(ID,3), app.filtresParam(ID,4), app.filtresParam(ID,5), app.fs, app.N);
                elseif app.filtresParam(ID,6)==3
                    app.Panell = PanellBanda(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.filtresParam(ID,3), app.filtresParam(ID,4), app.filtresParam(ID,5), app.fs, app.N);
                elseif app.filtresParam(ID,6)==4
                    app.Panell = PanellPB(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.fs, app.N);
                elseif app.filtresParam(ID,6)==5
                    app.Panell = PanellPA(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.fs, app.N);
                end
            else
                app.actualitza_eq;
            end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            stop(app.player);
            delete(app.Panell)
            delete(app)
        end

        % Menu selected function: AfegirfiltreMenu
        function AfegirfiltreMenuSelected(app, event)
            app.AfegirfiltrePanel.Visible = 'on';
        end

        % Button pushed function: AfegirButton
        function AfegirButtonPushed(app, event)
            if ~isempty(app.Panell)
                if isvalid(app.Panell)
                    app.Panell.valors;
                    delete(app.Panell)
                end
            end
            filter_type = app.TipusdefiltreDropDown.Value;
            ID = size(app.filtres, 1)+1;
            num = ' (' + string(app.nextID) + ')';

            switch filter_type
                case 'Bell'
                    fc=1000;
                    Q=1;
                    G=3;
                    app.filtresParam = [app.filtresParam; ID, fc, Q, G, 0, 1];
                    app.Panell = PanellBell(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.filtresParam(ID,3), app.filtresParam(ID,4), app.fs, app.N);
                case 'Shelf'
                    fc=676;
                    S=1;
                    G=3;
                    lH = 0;
                    app.filtresParam = [app.filtresParam; ID, fc, S, G, lH, 2];
                    app.Panell = PanellShelf(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.filtresParam(ID,3), app.filtresParam(ID,4), app.filtresParam(ID,5), app.fs, app.N);
                case 'Band'
                    fc=1000;
                    BW=2;
                    G=3;
                    Q=3;
                    app.filtresParam = [app.filtresParam; ID, fc, BW, Q, G, 3];
                    app.Panell = PanellBanda(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.filtresParam(ID,3), app.filtresParam(ID,4), app.filtresParam(ID,5), app.fs, app.N);
                case 'HighCut'
                    fc=20000;
                    app.filtresParam = [app.filtresParam; ID, fc, 0, 0, 0, 4];
                    app.Panell = PanellPB(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.fs, app.N);
                case 'LowCut'
                    fc=20;
                    app.filtresParam = [app.filtresParam; ID, fc, 0, 0, 0, 5];
                    app.Panell = PanellPA(app, app.filtresParam(ID,1), app.filtresParam(ID,2), app.fs, app.N);
            end
            nou =  filter_type + num;
            app.FiltreDropDown.Items = [app.FiltreDropDown.Items nou];
            app.FiltreDropDown.ItemsData = [app.FiltreDropDown.ItemsData ID];
            app.FiltreDropDown.Value = ID;
            app.nextID = app.nextID+1;
            app.AfegirfiltrePanel.Visible = 'off';
        end

        % Value changing function: AudioSlider
        function AudioSliderValueChanging(app, event)
            if isplaying(app.player)
                pause(app.player);
                playing = 1;
            else
                playing = 0;
            end
            
            changingValue = event.Value;
            app.audioWindow = changingValue;

            app.actualitza_ft;
            pause(0.05);

            if playing == 1
                a = app.audioWindow * app.N/2 + app.N/2;
                play(app.player,a);
            end
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
            
            if app.filtersOn == 1
                app.audioFiltrat = app.filtrar(app.audio(:,1), app.H, app.N*2, app.N);
                if size(app.audio, 2)==2
                    app.audioFiltrat= [app.audioFiltrat; app.filtrar(app.audio(:,2), app.H, app.N*2, app.N)];
                end
                app.player = audioplayer(app.audioFiltrat, app.fs);
            else
                app.player = audioplayer(app.audio, app.fs);
            end
            a = app.audioWindow * app.N/2 + app.N/2;
            play(app.player,a);
            
            while(isplaying(app.player))
                newV = round(app.player.CurrentSample*2/app.N);
                if newV~=0
                    if newV < app.AudioSlider.Limits(1, 2)
                        app.AudioSlider.Value = newV;
                    end
                end
                app.audioWindow = app.AudioSlider.Value;
                app.actualitza_ft;
                pause(0.05);
            end
        end

        % Button pushed function: PauseButton
        function PauseButtonPushed(app, event)
            pause(app.player);
        end

        % Button pushed function: XButton
        function XButtonPushed(app, event)
            app.AfegirfiltrePanel.Visible = 'off';
        end

        % Callback function: addButton
        function addButtonClicked(app, event)
            app.AfegirfiltrePanel.Visible = 'on';
        end

        % Callback function: editButton
        function editButtonClicked(app, event)
            app.FiltreDropDownValueChanged;
        end

        % Callback function: deleteButton
        function deleteButtonClicked(app, event)
            ID = app.FiltreDropDown.Value;
            if ~isempty(app.Panell)
                if isvalid(app.Panell)
                    app.Panell.valors;
                    delete(app.Panell);
                end
            end
            if ID~=0
                [s1, s2] = size(app.filtres);
                if s1~=1
                    i=ID;
                    if s1~=ID
                        for i= ID+1:s1
                            app.filtres(i-1,1:s2) = app.filtres(i,1:s2);
                            app.filtresParam(i-1,2:6) = app.filtresParam(i,2:6);
                            app.FiltreDropDown.Items(i) = app.FiltreDropDown.Items(i+1);
                        end
                    end
                    app.FiltreDropDown.ValueIndex = 1;
                    app.FiltreDropDown.Items = app.FiltreDropDown.Items(1:i);
                    app.FiltreDropDown.ItemsData = app.FiltreDropDown.ItemsData(1:i);
                    app.filtres = app.filtres(1:i-1,:);
                    app.filtresParam = app.filtresParam(1:i-1,:);
                    app.actualitza_eq;
                else
                    app.filtres = double.empty;
                    app.filtresParam = double.empty;
                    app.FiltreDropDown.ItemsData=0;
                    app.FiltreDropDown.Items = {'Selecciona un filtre'};
                    app.actualitza_eq;
                    
                end
                
            end
        end

        % Callback function: exportButton
        function exportButtonClicked(app, event)
            app.ExportMode = 'audio';
            app.Export.Title = 'Export audio';
            app.FileNameEditField.Value = app.audioName;
            app.FormatDropDown.Items = {'wav', 'mp3'};
            app.FormatDropDown.Value = 'wav';
            app.Export.Visible = 'on';
        end

        % Button pushed function: ExporttoButton
        function ExporttoButtonPushed(app, event)
            carpeta = uigetdir; 
            if carpeta ~= 0
                nomFitxer = app.FileNameEditField.Value;
                formatFitxer = app.FormatDropDown.Value;
        
                [~,~,ext] = fileparts(nomFitxer);
                if isempty(ext)
                    nomFitxer = [nomFitxer, '.', formatFitxer];
                end
        
                fullPath = fullfile(carpeta, nomFitxer);
        
                switch app.ExportMode
                    case 'audio'
                        audiowrite(fullPath, app.audio, app.fs);
        
                    case 'param'
                        param = app.filtresParam;
                        save(fullPath, 'param');
                end
        
                app.Export.Visible = 'off';
            end
        end

        % Value changed function: FileNameEditField
        function FileNameEditFieldValueChanged(app, event)
            value = app.FileNameEditField.Value;
            app.audioName = value;
        end

        % Button pushed function: XButton_2
        function XButton_2Pushed(app, event)
            app.Export.Visible = 'off';
        end

        % On callback: ToggleTool
        function ToggleToolOn(app, event)
            app.filtersOn = 1;
            if isplaying(app.player)
                sample = app.player.CurrentSample;
                app.player = audioplayer(app.audioFiltrat, app.fs);
                play(app.player, sample)
            end
        end

        % Off callback: ToggleTool
        function ToggleToolOff(app, event)
            app.filtersOn = 0;
            if isplaying(app.player)
                sample = app.player.CurrentSample;
                app.player = audioplayer(app.audio, app.fs);
                play(app.player, sample)
            end
        end

        % Callback function: saveButton
        function saveButtonClicked(app, event)
            app.ExportMode = 'param';
            app.Export.Title = 'Save filter parameters';
            app.FileNameEditField.Value = 'datos';
            app.FormatDropDown.Items = {'mat'}; % Només .mat per ara
            app.FormatDropDown.Value = 'mat';
            app.Export.Visible = 'on';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.149 0.149 0.149];
            app.UIFigure.Position = [180 30 1000 665];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Resize = 'off';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.Pointer = 'crosshair';

            % Create Toolbar
            app.Toolbar = uitoolbar(app.UIFigure);
            app.Toolbar.BackgroundColor = [0.102 0.102 0.102];

            % Create ToggleTool
            app.ToggleTool = uitoggletool(app.Toolbar);
            app.ToggleTool.State = 'on';
            app.ToggleTool.Icon = fullfile(pathToMLAPP, 'Filtrar.png');
            app.ToggleTool.OffCallback = createCallbackFcn(app, @ToggleToolOff, true);
            app.ToggleTool.OnCallback = createCallbackFcn(app, @ToggleToolOn, true);

            % Create addButton
            app.addButton = uipushtool(app.Toolbar);
            app.addButton.ClickedCallback = createCallbackFcn(app, @addButtonClicked, true);
            app.addButton.Icon = fullfile(pathToMLAPP, 'Afegir.png');
            app.addButton.Separator = 'on';

            % Create editButton
            app.editButton = uipushtool(app.Toolbar);
            app.editButton.ClickedCallback = createCallbackFcn(app, @editButtonClicked, true);
            app.editButton.Icon = fullfile(pathToMLAPP, 'Editar.png');

            % Create deleteButton
            app.deleteButton = uipushtool(app.Toolbar);
            app.deleteButton.ClickedCallback = createCallbackFcn(app, @deleteButtonClicked, true);
            app.deleteButton.Icon = fullfile(pathToMLAPP, 'Eliminar.png');

            % Create saveButton
            app.saveButton = uipushtool(app.Toolbar);
            app.saveButton.ClickedCallback = createCallbackFcn(app, @saveButtonClicked, true);
            app.saveButton.Icon = fullfile(pathToMLAPP, 'Guardar.png');
            app.saveButton.Separator = 'on';

            % Create exportButton
            app.exportButton = uipushtool(app.Toolbar);
            app.exportButton.ClickedCallback = createCallbackFcn(app, @exportButtonClicked, true);
            app.exportButton.Icon = fullfile(pathToMLAPP, 'Exportar.png');

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.UIFigure);
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.AmbientLightColor = 'none';
            app.UIAxes_2.XLim = [10 30000];
            app.UIAxes_2.YLim = [-150 0];
            app.UIAxes_2.GridLineWidth = 0.25;
            app.UIAxes_2.MinorGridLineWidth = 0.25;
            app.UIAxes_2.MinorGridLineStyle = '-';
            app.UIAxes_2.XColor = 'none';
            app.UIAxes_2.XTick = [];
            app.UIAxes_2.XScale = 'log';
            app.UIAxes_2.XTickLabel = '';
            app.UIAxes_2.XMinorTick = 'on';
            app.UIAxes_2.YAxisLocation = 'right';
            app.UIAxes_2.YColor = [0.651 0.651 0.651];
            app.UIAxes_2.YTick = [-150 -100 -50 0];
            app.UIAxes_2.ZColor = [0.9412 0.9412 0.9412];
            app.UIAxes_2.Color = [0.251 0.251 0.251];
            app.UIAxes_2.GridColor = [1 1 0.0667];
            app.UIAxes_2.MinorGridColor = [0.502 0.502 0.502];
            app.UIAxes_2.GridAlpha = 0.25;
            app.UIAxes_2.ColorOrder = [0.149019607843137 0.149019607843137 0.149019607843137;0.501960784313725 0.501960784313725 0.501960784313725;0.929411764705882 0.694117647058824 0.125490196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.466666666666667 0.674509803921569 0.188235294117647;0.301960784313725 0.745098039215686 0.933333333333333;0.635294117647059 0.0784313725490196 0.184313725490196];
            colormap(app.UIAxes_2, 'colorcube')
            app.UIAxes_2.Position = [14 60 978 604];

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.AmbientLightColor = [0.9412 0.9412 0.9412];
            app.UIAxes.XLim = [10 24000];
            app.UIAxes.YLim = [-24 12];
            app.UIAxes.GridLineWidth = 0.25;
            app.UIAxes.MinorGridLineWidth = 0.25;
            app.UIAxes.MinorGridLineStyle = '-';
            app.UIAxes.XColor = [0.9412 0.9412 0.9412];
            app.UIAxes.XTick = [20 50 100 200 500 1000 2000 5000 10000 20000];
            app.UIAxes.XScale = 'log';
            app.UIAxes.XTickLabel = {'20'; '50'; '100'; '200'; '500'; '1k'; '2k'; '5k'; '10k'; '20k'};
            app.UIAxes.XMinorTick = 'on';
            app.UIAxes.YAxisLocation = 'right';
            app.UIAxes.YColor = [0.9412 0.9412 0.9412];
            app.UIAxes.YTick = [-24 -21 -18 -15 -12 -9 -6 -3 0 3 6 9 12];
            app.UIAxes.ZColor = [0.9412 0.9412 0.9412];
            app.UIAxes.BoxStyle = 'full';
            app.UIAxes.Color = 'none';
            app.UIAxes.GridColor = [0.9412 0.9412 0.9412];
            app.UIAxes.MinorGridColor = [0.502 0.502 0.502];
            app.UIAxes.GridAlpha = 0.25;
            app.UIAxes.XGrid = 'on';
            app.UIAxes.XMinorGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ColorOrder = [0 1 1;1 0 1;1 1 0;0 1 0;1 0 0;0 0 1];
            colormap(app.UIAxes, 'colorcube')
            app.UIAxes.Position = [14 49 944 615];

            % Create AfegirfiltrePanel
            app.AfegirfiltrePanel = uipanel(app.UIFigure);
            app.AfegirfiltrePanel.ForegroundColor = [0.9412 0.9412 0.9412];
            app.AfegirfiltrePanel.Title = 'Afegir filtre';
            app.AfegirfiltrePanel.Visible = 'off';
            app.AfegirfiltrePanel.BackgroundColor = [0.149 0.149 0.149];
            app.AfegirfiltrePanel.Position = [371 283 260 118];

            % Create TipusdefiltreDropDownLabel
            app.TipusdefiltreDropDownLabel = uilabel(app.AfegirfiltrePanel);
            app.TipusdefiltreDropDownLabel.HorizontalAlignment = 'right';
            app.TipusdefiltreDropDownLabel.FontColor = [0.9412 0.9412 0.9412];
            app.TipusdefiltreDropDownLabel.Position = [28 55 77 22];
            app.TipusdefiltreDropDownLabel.Text = 'Tipus de filtre';

            % Create TipusdefiltreDropDown
            app.TipusdefiltreDropDown = uidropdown(app.AfegirfiltrePanel);
            app.TipusdefiltreDropDown.Items = {'Bell', 'Shelf', 'Band', 'HighCut', 'LowCut'};
            app.TipusdefiltreDropDown.FontColor = [0.9412 0.9412 0.9412];
            app.TipusdefiltreDropDown.BackgroundColor = [0.149 0.149 0.149];
            app.TipusdefiltreDropDown.Position = [120 55 100 22];
            app.TipusdefiltreDropDown.Value = 'Bell';

            % Create AfegirButton
            app.AfegirButton = uibutton(app.AfegirfiltrePanel, 'push');
            app.AfegirButton.ButtonPushedFcn = createCallbackFcn(app, @AfegirButtonPushed, true);
            app.AfegirButton.BackgroundColor = [0.149 0.149 0.149];
            app.AfegirButton.FontColor = [0.9412 0.9412 0.9412];
            app.AfegirButton.Position = [74 19 100 23];
            app.AfegirButton.Text = '+ Afegir';

            % Create XButton
            app.XButton = uibutton(app.AfegirfiltrePanel, 'push');
            app.XButton.ButtonPushedFcn = createCallbackFcn(app, @XButtonPushed, true);
            app.XButton.BackgroundColor = [0.149 0.149 0.149];
            app.XButton.FontColor = [0.9412 0.9412 0.9412];
            app.XButton.Position = [241 99 19 19];
            app.XButton.Text = 'X';

            % Create colorrgb094094094FiltreLabel
            app.colorrgb094094094FiltreLabel = uilabel(app.UIFigure);
            app.colorrgb094094094FiltreLabel.BackgroundColor = [0.149 0.149 0.149];
            app.colorrgb094094094FiltreLabel.HorizontalAlignment = 'right';
            app.colorrgb094094094FiltreLabel.FontColor = [0.9412 0.9412 0.9412];
            app.colorrgb094094094FiltreLabel.Position = [25 619 32 22];
            app.colorrgb094094094FiltreLabel.Text = 'Filtre';

            % Create FiltreDropDown
            app.FiltreDropDown = uidropdown(app.UIFigure);
            app.FiltreDropDown.Items = {'Selecciona un filtre'};
            app.FiltreDropDown.ItemsData = 0;
            app.FiltreDropDown.ValueChangedFcn = createCallbackFcn(app, @FiltreDropDownValueChanged, true);
            app.FiltreDropDown.FontColor = [0.9412 0.9412 0.9412];
            app.FiltreDropDown.BackgroundColor = [0.149 0.149 0.149];
            app.FiltreDropDown.Position = [72 619 100 22];
            app.FiltreDropDown.Value = 0;

            % Create AudioSlider
            app.AudioSlider = uislider(app.UIFigure);
            app.AudioSlider.MajorTicks = [];
            app.AudioSlider.MajorTickLabels = {''};
            app.AudioSlider.ValueChangingFcn = createCallbackFcn(app, @AudioSliderValueChanging, true);
            app.AudioSlider.FontColor = [0.9412 0.9412 0.9412];
            app.AudioSlider.Position = [97 32 887 3];

            % Create Button
            app.Button = uibutton(app.UIFigure, 'push');
            app.Button.Position = [-29 694 2 2];

            % Create PlayButton
            app.PlayButton = uibutton(app.UIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.BackgroundColor = [0.149 0.149 0.149];
            app.PlayButton.FontColor = [0.9412 0.9412 0.9412];
            app.PlayButton.Position = [14 34 59 28];
            app.PlayButton.Text = 'Play';

            % Create PauseButton
            app.PauseButton = uibutton(app.UIFigure, 'push');
            app.PauseButton.ButtonPushedFcn = createCallbackFcn(app, @PauseButtonPushed, true);
            app.PauseButton.BackgroundColor = [0.149 0.149 0.149];
            app.PauseButton.FontColor = [0.9412 0.9412 0.9412];
            app.PauseButton.Position = [14 6 59 26];
            app.PauseButton.Text = 'Pause';

            % Create Export
            app.Export = uipanel(app.UIFigure);
            app.Export.ForegroundColor = [0.9412 0.9412 0.9412];
            app.Export.Title = 'Export audio';
            app.Export.Visible = 'off';
            app.Export.BackgroundColor = [0.149 0.149 0.149];
            app.Export.Position = [370 237 260 163];

            % Create FileNameEditFieldLabel
            app.FileNameEditFieldLabel = uilabel(app.Export);
            app.FileNameEditFieldLabel.BackgroundColor = [0.149 0.149 0.149];
            app.FileNameEditFieldLabel.HorizontalAlignment = 'right';
            app.FileNameEditFieldLabel.FontColor = [0.9412 0.9412 0.9412];
            app.FileNameEditFieldLabel.Position = [25 110 60 22];
            app.FileNameEditFieldLabel.Text = 'File Name';

            % Create FileNameEditField
            app.FileNameEditField = uieditfield(app.Export, 'text');
            app.FileNameEditField.ValueChangedFcn = createCallbackFcn(app, @FileNameEditFieldValueChanged, true);
            app.FileNameEditField.FontColor = [0.9412 0.9412 0.9412];
            app.FileNameEditField.BackgroundColor = [0.149 0.149 0.149];
            app.FileNameEditField.Position = [100 110 121 22];

            % Create ExporttoButton
            app.ExporttoButton = uibutton(app.Export, 'push');
            app.ExporttoButton.ButtonPushedFcn = createCallbackFcn(app, @ExporttoButtonPushed, true);
            app.ExporttoButton.BackgroundColor = [0.149 0.149 0.149];
            app.ExporttoButton.FontColor = [0.9412 0.9412 0.9412];
            app.ExporttoButton.Position = [75 24 100 23];
            app.ExporttoButton.Text = 'Export to...';

            % Create XButton_2
            app.XButton_2 = uibutton(app.Export, 'push');
            app.XButton_2.ButtonPushedFcn = createCallbackFcn(app, @XButton_2Pushed, true);
            app.XButton_2.BackgroundColor = [0.149 0.149 0.149];
            app.XButton_2.FontColor = [0.9412 0.9412 0.9412];
            app.XButton_2.Position = [240 143 19 19];
            app.XButton_2.Text = 'X';

            % Create TipusLabel
            app.TipusLabel = uilabel(app.Export);
            app.TipusLabel.BackgroundColor = [0 0 0];
            app.TipusLabel.HorizontalAlignment = 'right';
            app.TipusLabel.FontColor = [1 1 1];
            app.TipusLabel.Position = [51 82 34 22];
            app.TipusLabel.Text = 'Tipus';

            % Create FormatDropDown
            app.FormatDropDown = uidropdown(app.Export);
            app.FormatDropDown.Items = {'mp3', 'wav'};
            app.FormatDropDown.FontColor = [1 1 1];
            app.FormatDropDown.BackgroundColor = [0 0 0];
            app.FormatDropDown.Position = [105 81 100 22];
            app.FormatDropDown.Value = 'mp3';

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create AfegirfiltreMenu
            app.AfegirfiltreMenu = uimenu(app.ContextMenu);
            app.AfegirfiltreMenu.MenuSelectedFcn = createCallbackFcn(app, @AfegirfiltreMenuSelected, true);
            app.AfegirfiltreMenu.Text = '+ Afegir filtre';
            
            % Assign app.ContextMenu
            app.UIAxes_2.ContextMenu = app.ContextMenu;
            app.UIAxes.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Eqcodi(varargin)

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