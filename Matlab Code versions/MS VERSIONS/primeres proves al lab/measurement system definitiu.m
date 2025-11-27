%% Configuració
clear all

% Mode: 'synthetic' o 'recorded'
mode = 'synthetic';

% Tipus de senyal si és sintètica
signal_type = 'saturated';  % 'sine', 'square', 'triangle', 'chirp', 'saturated'

% Paràmetres de senyal
FS       = 44100;           % Freqüència de mostreig (Hz)
duration = 1.0;             % Durada de la senyal (s)
F        = 440;             % Freqüència base (Hz)

t = 0:(1/FS):duration-(1/FS);    % Vector de temps
n = 0:1:duration*FS-(1/FS);      % Vector d’índexs

%% Generació de senyal
switch signal_type
    case 'sine'
        x = sin(2*pi*F*t)';
    case 'square'
        x = square(2*pi*F*t)';
    case 'triangle'
        x = 2 * abs(2 * mod(F * t, 1) - 1) - 1;
        x = x';
    case 'chirp'
        durationch     = 5;      
        Nbits        = 16;      

        Fch       = 20000;                              
        samples = durationch*FS;
        fmi     = (0:1:samples-1)'./(samples-1);      
        phi       = pi*(Fch/FS)*(0:samples-1)';
        x  = sin(phi.*fmi);  
    case 'saturated'
        x_sat = 5 * sin(2*pi*F*t);  % Amplitud exageradament alta
        % Saturació: retallar a [-1, 1]
        x = max(min(x_sat, 1), -1);
    otherwise
        error('Tipus de senyal no reconegut');
end

%% Enregistrament o assignació segons mode
if strcmp(mode, 'recorded')
    a_r = audiorecorder(FS, 24, 2, 0);
    t_rec = duration;

    disp("Preparat! Gravació en 3 segons...");
    pause(3);

    disp("Recording...");
    recordblocking(a_r, t_rec);
    disp("Recording finished");

    y = getaudiodata(a_r);

    % Convertir a mono si té dos canals
    if size(y, 2) == 2
        y = mean(y, 2);
    end

    % Detecció de llindar
    threshold = 0.02;  % Ajustable
    start_idx = find(abs(y) > threshold, 1, 'first');

    if isempty(start_idx)
        warning("No s'ha detectat activitat a la gravació");
    else
        y = y(start_idx:end);  % Tallar silenci inicial
    end
elseif strcmp(mode, 'synthetic')
    y = x(:);
else
    error("Mode no reconegut. Usa 'synthetic' o 'recorded'.");
end
%% Càlcul Potència
P = pow(y);
fprintf('Potència senyal: %.2f\n', P);
%% Càlcul Guany
G = gainTest(y, FS);
fprintf('Guany senyal: %.2f\n', G);
%% Reproducció de la senyal
kk = audioplayer(y, FS);
play(kk);

%%
% FFT
Y = fft(y);
N = length(Y);
f = (0:N-1)*(FS/N);

% Només meitat positiva
Y_half = Y(1:floor(N/2));
f_half = f(1:floor(N/2));
mag = abs(Y_half);

% Fonamental
[~, idx_fund] = max(mag);
f_fund = f_half(idx_fund);
P1 = mag(idx_fund)^2;

% Buscar harmònics propers a la freq esperada
num_harmonics = 5;
harmonics = zeros(1, num_harmonics - 1);

for k = 2:num_harmonics
    target_freq = f_fund * k;
    [~, idx] = min(abs(f_half - target_freq));
    harmonics(k-1) = mag(idx)^2;
end

% THD en dB
if P1 == 0
    THD_value = -Inf;
else
    THD_value = 10*log10(sum(harmonics) / P1);
end

fprintf("THD: %.2f dB\n", THD_value);

%% Visualització
NFFT    = 1024;  
fc = NFFT / 2;            
win_len = NFFT;           
a = 10;                   
L = NFFT;                 
overlap = 0.25 * NFFT;
hop     = NFFT - overlap;
win = 0.54 - 0.46 * cos(2 * pi * (0:NFFT-1)' / (NFFT - 1));

% Normalització de la finestra per mantenir energia
win = win / sqrt(sum(win.^2));

% Nombre de frames
numFrames = floor((length(y) - NFFT) / hop) + 1;

% Preallocació
S = zeros(NFFT/2, numFrames);  % Magnituds (només part positiva de l'FFT)
T = zeros(1, numFrames);       % Temps

% Calcula espectre per a cada finestra
if length(y) < NFFT
    error('La senyal és massa curta per calcular l’espectrograma amb NFFT = %d', NFFT);
end

for i = 1:numFrames
    startIdx = (i-1)*hop + 1;
    endIdx   = startIdx + NFFT - 1;

    frame = y(startIdx:endIdx, 1) .* win;  % Aplicar finestra
    Y = fft(frame, NFFT);                  % FFT

    S(:,i) = abs(Y(1:NFFT/2));             % Només part positiva
    T(i) = startIdx / FS;                  % Temps del frame
end

% Vector de freqüències
f = FS * (0:(NFFT/2)-1) / NFFT;

% Visualització

% Normalització del THD per color
if isinf(THD_value)
    THD_norm = 0;  % Sense distorsió detectable, posar color blau
elseif THD_value <= -60  % Si el THD és menor o igual a -60 dB, color blau (màxima netedat)
    THD_norm = 0;  % Color blau
else
    % Normalització en el rang [–60 dB, 0 dB] per la distorsió
    THD_norm = (THD_value + 60) / 60;  % Convertir de –60 a 0 en un rang de 0 a 1
end

% Color interpolat entre blau (net) i vermell (distorsionat)
color_signal = (1 - THD_norm)*[0 0 1] + THD_norm*[1 0 0];

% Plot
figure(1)
subplot(2,1,1)
plot((1:length(y))/FS, y, 'Color', color_signal, 'LineWidth', 1.5);
title(sprintf('Forma d’ona (THD = %.2f dB)', THD_value));
xlabel('Temps (s)');
ylabel('Amplitud');




subplot(2,1,2)
surf(T, f, 20*log10(S + eps), 'edgecolor', 'none');  % eps per evitar log(0)
axis tight;
view(0, 90);
xlabel('Temps (s)');
ylabel('Freqüència (Hz)');
title('Espectrograma manual');
colorbar;

%% Funcions

function P = pow(x)
    P = (x'*x)/length(x);
end

function G = gain(x, xamp)
    G = 10*log10(pow(xamp)/pow(x));
end

function GT = gainTest(x, FS)
    recTime = 1;
    recorder = audiorecorder(FS,24,2,0);
    recordblocking(recorder, recTime);
    xamp = getaudiodata(recorder);
    GT = gain(x, xamp);
end
