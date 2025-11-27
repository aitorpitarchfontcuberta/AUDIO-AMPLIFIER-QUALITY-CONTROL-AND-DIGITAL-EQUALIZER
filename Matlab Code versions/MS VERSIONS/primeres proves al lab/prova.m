%% 
FS = 48000;
N = 1*FS;
n = (0:N-1);
F = 430;
x = sin(2*pi*(F/FS)*n);
stem(n,x);
kk = audioplayer(x,FS);
play(kk);

%% 
a_r =audiorecorder(FS, 24, 1);
t_rec = 3;
disp("Recording...");
recordblocking(a_r, t_rec);
disp("Recording finished2");

%% 

%play(a_r);
y=getaudiodata(a_r);
kk2=audioplayer(y,FS);
play(kk2);
pause(2);
kk3=audioplayer(y,1.5*FS);
play(kk3);
pause(1);
kk4=audioplayer(y,2*FS);
play(kk4);


%% 
NFFT      = 2^10;
window    = hamming(NFFT);
overlap   = 0.25 * NFFT;

figure(1)

subplot(2,1,1)
plot(y);

[S,FY,T,P1] = spectrogram(y(:,1), hamming(NFFT), overlap, NFFT, Fs);
subplot(2,1,2)
surf(T,FY,10*log10(P1), 'edgecolor', 'none'); axis tight;
view(0, 90);