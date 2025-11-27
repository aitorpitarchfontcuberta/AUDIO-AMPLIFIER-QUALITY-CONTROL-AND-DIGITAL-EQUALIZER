clear all
Fs           = 44100;            % Sampling Rate (check with your computer/sound board doc)
duration     = 2.0;              % Signal duration time in seconds
F = 440;                         % Frequency in Hz. 


t = 0:(1/Fs):duration-(1/Fs);    % Axis time.
n = 0:1:duration*Fs-(1/Fs);      % Index vector.
x = sin(2*pi*n*(F/Fs))';         % Discrete signal generation.


SNR=snrTest(x,Fs)

%%
audiowrite('sin.wav',x,Fs);

y=audioread('sin.wav');

P=mean(y.^2)



function P=power(x) 
    P=(x'*x)/length(x);
end

function G=gain(x, xamp)
    G=10*log10(power(xamp)/power(x));
end
    
function GT=gainTest(x,Fs)
    recTime=1;
    recorder=audiorecorder(Fs,24,1);
    sound(x,Fs);
    pause(0.1)
    recordblocking(recorder,recTime);
    xamp=getaudiodata(recorder);
    GT=gain(x, xamp);
    
end 

function SNR=snrTest(x,Fs)
    recTime=1;
    signalRecorder=audiorecorder(Fs,24,1);
    sound(x,Fs);
    pause(0.1)
    recordblocking(signalRecorder,recTime);
    xamp=getaudiodata(signalRecorder);
    Ps=power(xamp);
    
    noiseRecorder=audiorecorder(Fs,24,1);
    recordblocking(noiseRecorder,recTime);
    noise=getaudiodata(noiseRecorder);
    Pn=power(noise);
    
    SNR=Ps/Pn;
end 