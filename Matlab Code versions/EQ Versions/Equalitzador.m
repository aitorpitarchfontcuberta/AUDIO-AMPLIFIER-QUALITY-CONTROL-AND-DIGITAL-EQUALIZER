%% Equalitzador
clear all
%%

fs=48000; % Freqüència de mostratge
duration     = 1;              % LIMITACIÓ a 1s

%% -------------------------FASE D'ENTRADA--------------------------------
%% SENYALS DE PROVA
F = 440;                         
F2 = 2.*F;
F3 = 2.*F2;
F4 = 2.*F3;

t = 0:(1/fs):duration-(1/fs);    % eix temporal
n = 0:1:duration*fs-(1/fs);      
so = (sin(2*pi*n*(F/fs))+sin(2*pi*n*(F2/fs))+sin(2*pi*n*(F3/fs))+sin(2*pi*n*(F4/fs)))./4;  

%% SENYAL GRAVADA
a_r =audiorecorder(fs, 16, 1);
disp("Recording...");
recordblocking(a_r, duration);
disp("Recording finished");

so = getaudiodata(a_r)';

%% ARXIU D'AUDIO

so1 = audioread("so1.mp3");
so = so1(1:fs);

%% -----------------------MUNTATGE DE FILTRES-----------------------------

% Coeficients del filtre pass baix chebyshev tipus II
Num =[0.348415248420542,2.739630135368878,9.471519386413812,18.804022197693786,23.447435457214436,18.804022197693786,9.471519386413814,2.739630135368878,0.348415248420542];
Den =[1,5.798820778818033,15.018457076581917,22.629115506603632,21.653557540257474,13.453944119231329,5.294200683421229,1.205120502757305,0.121393185337581];


%% FILTRES GAUSSIANS (PARAMÈTRICS)
% Els paràmetres a passar en ordre son, freq. central, ample de banda,
% amplificació en dB i freq. de mostratge
y0=filtre_gauss(15, 5, 0, fs);
y1=filtre_gauss(29, 10, 0, fs);
y2=filtre_gauss(58, 20, 0, fs);
y3=filtre_gauss(117, 39, 0, fs);
y4=filtre_gauss(234, 78, 0, fs);

y5=filtre_gauss(468, 156, 6, fs);

y6=filtre_gauss(937, 312, 0, fs);
y7=filtre_gauss(1875, 625, 0, fs);
y8=filtre_gauss(3750, 1250, 0, fs);

% freq eliminada : Freqüència central, banda, amplitud, freq de mostratge
y9=filtre_gauss(880, 200, 6, fs);

% banda eliminada : Freqüència 1, 2, banda 1, 2, amplitud, freq de mostratge
y10=filtre_banda_gauss(1760, 3520, 200, 200, 6, fs);


% Unió en un banc de filtres
y=y0.*y1.*y2.*y3.*y4.*y5.*y6.*y7.*y8.*y9.*y10;



% Unió entre filtre passa baix fixe i banc de filtres paramètrics
H = h_filtre_chebyII(Num, Den, fs);
H1=abs(H).*y;

%% Mostra de filtres

x =linspace(0, fs/2, fs/2); 
semilogx(x, 10*log10(abs(H)), "r", x, 10*log10(abs(y5)), "b", x, 10*log10(abs(y9)), "g", x, 10*log10(abs(y10)), "m");
ylim([-60, 30]);
xlim([0, fs/2]);

%% ---------------------FILTRATGE I SORTIDA-------------------------------
% filtrem el senyal
so4 = filter_meu(so, H1, fs);

%El Reproduim i comparem amb l'original
%% Senyal filtrat

sound(so4, fs);

%% Senyal original

sound(so, fs);

%% Funcions

function y=gauss_meu(mu, sigma, a, fs)
% vector de fs/2 posicions
x =linspace(0, fs/2, fs/2);

% finestra gaussiana d'amplitud a, centre a mu i desviació sigma
y = a*(exp(-(x - mu).^2 / (2 * sigma^2)));

end

function feq=filtre_gauss(fc, ftall, A, fs)
    % sigma de la finestra de gauss amb frequència de tall ftall
    sigma = ftall/2; 

    % finestra de gauss corresponent
    f=gauss_meu(fc, sigma, A, fs);

    % equivalent en escala logaritmica
    feq=10.^(f/10);

end

function feq=filtre_gauss_asim(fc, ftall1, ftall2, A, fs)
    % sigma de la finestra de gauss amb frequència de tall ftall
    sigma1 = ftall1/2;
    sigma2 = ftall2/2; 

    % finestra de gauss corresponent
    f=doble_gauss_meu(fc, sigma1, sigma2, A, fs);

    % equivalent en escala logaritmica
    feq=10.^(f/10);

end

function H = h_filtre_chebyII(Num, Den, fs)
    % Vector de freqüències normalitzades (rad/s)
    w = linspace(0, pi, fs/2);

    % càlcul de la funció de transferència
    H = zeros(size(w));
    for k = 1:length(w)
       z = exp(-1j*w(k));        
      H(k) = polyval(Num, z) ./ polyval(Den, z);
    end
end

function so_filtrat=filter_meu(so, H, fs)
% FFT del senyal d'entrada
Fso=fft(so, fs);

% FFT de 0 a fs/2
Fso_meitat=Fso(1:end/2);

% Filtratge del senyal
Fso_filtrat=Fso_meitat.*H;

% IFFT del senyal filtrat
so_filtrat=real(ifft(Fso_filtrat, fs).*2);


% Gràfiques
x =linspace(0, fs/2, fs/2); 
figure;
subplot(2, 2, 1)
plot(so);
subplot(2, 2, 2)

semilogx(10*log10(abs(Fso_meitat)*2/fs));
ylim([-60, 30]);
xlim([0, fs/2]);
subplot(2, 2, 3)
plot(real(so_filtrat));
subplot(2, 2, 4)

semilogx(x, 10*log10(abs(H)), "r", x, 10*log10(abs(Fso_filtrat)*2/fs), "b")
ylim([-60, 30]);
xlim([0, fs/2]);

end

function weq = hamming_meu(fc, win_len, a, L)

    % Crear finestra de Hamming normalitzada
    n = (0:win_len-1);
    h = 0.5 - 0.5 * cos(2*pi*n/(win_len-1));
    h = h * a;  % Escalar l’amplitud

    % Inserir la finestra dins el vector
    w = zeros(1, L);
    start_pos = round(fc - (win_len/2));
    end_pos = start_pos + win_len - 1;

    % Ajustar si surt dels límits
    win_start = 1;
    win_end = win_len;
    if start_pos < 1
        win_start = 2 - start_pos;
        start_pos = 1;
    end
    if end_pos > fs
        win_end = win_len - (end_pos - L);
        end_pos = L;
    end

    w(start_pos:end_pos) = h(win_start:win_end);

    weq=10.^(w/10);
end

function y=doble_gauss_meu(mu, sigma1, sigma2, a, fs)

% Rango de valores x
x =linspace(0, fs/2, fs/2); 

% Función Gaussiana normalizada con amplitud máxima de 1
y1 = a*(exp(-(x - mu).^2 / (2 * sigma1^2)));
y2 = a*(exp(-(x - mu).^2 / (2 * sigma2^2)));
y=[y1(1:mu) y2(mu+1:fs/2)];
end

function f=filtre_banda_gauss(ft1, ft2, bi, bf, a, fs)

b=ft2-ft1;
ya=filtre_gauss_asim(ft1, bi, b, a-a*0.175, fs);

yb=filtre_gauss_asim(ft2, b, bf, a-a*0.175, fs);

f=ya.*yb;
end
