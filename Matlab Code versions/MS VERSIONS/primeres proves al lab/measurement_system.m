%% 
clear all
FS           = 44100;            % Sampling Rate (check with your computer/sound board doc)
duration     = 2.0;              % Signal duration time in seconds
F = 440;                         % Frequency in Hz. 

t = 0:(1/FS):duration-(1/FS);    % Axis time.
n = 0:1:duration*FS-(1/FS);      % Index vector.
x = sin(2*pi*n*(F/FS))';         % Discrete signal generation.

SNR=snrTest(x,FS)

%% 
a_r =audiorecorder(FS, 24, 1);
t_rec = 3;
disp("Recording...");
recordblocking(a_r, t_rec);
disp("Recording finished2");

%% 
%play(a_r);
y=getaudiodata(a_r);
kk=audioplayer(y,FS);
play(kk);

%% 
NFFT      = 2^10;
window    = hamming(NFFT);
overlap   = 0.25 * NFFT;

figure(1)

subplot(2,1,1)
plot(y);

[S,FY,T,P1] = spectrogram(y(:,1), hamming(NFFT), overlap, NFFT, FS);
subplot(2,1,2)
surf(T,FY,10*log10(P1), 'edgecolor', 'none'); axis tight;
view(0, 90);

%%

audiowrite('sin.wav',x,FS);

y=audioread('sin.wav');

pow(y)

function P=pow(x) 
    P=(x'*x)/length(x);
end

function G=gain(x, xamp)
    G=10*log10(pow(xamp)/pow(x));
end
    
function GT=gainTest(x,FS)
    recTime=1;
    recorder=audiorecorder(FS,24,1);
    sound(x,FS);
    pause(0.1)
    recordblocking(recorder,recTime);
    xamp=getaudiodata(recorder);
    GT=gain(x, xamp);
    
end 

function SNR=snrTest(x,FS)
    recTime=1;
    signalRecorder=audiorecorder(FS,24,1);
    sound(x,FS);
    pause(0.1)
    recordblocking(signalRecorder,recTime);
    xamp=getaudiodata(signalRecorder);
    Ps=pow(xamp);
    
    noiseRecorder=audiorecorder(FS,24,1);
    recordblocking(noiseRecorder,recTime);
    noise=getaudiodata(noiseRecorder);
    Pn=pow(noise);
    
    SNR=Ps/Pn;
end 