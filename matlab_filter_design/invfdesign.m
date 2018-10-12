% adesigng44k.m
%


FS=22050
Nb = 6
Na = 6


if FS==22050
    % for 22050 we toss out fitting past 11 kHz
    freq_weights = [15, 20, 25, 32, 50, 63,  125, 250, 500, 1000, 2000, 4000, 8000, 10000, 11000] ;
    weights = [1000, 1000, 1000, 300, 2000, 200, 100, 100, 100, 100, 10, 10, 10, 1, 01];
else
    freq_weights = [15, 20, 25, 32, 50, 63,  125, 250, 500, 1000, 2000, 4000, 8000, 10000, 16000, 18000, 19000, 20000] ;
    weights = [600, 6000, 6000, 3000, 2000, 200, 10, 10, 10, 1, 1, 1, 1, 1, 0.1, 0.1, 1, 1.1];
end



fmax = FS/2;
fs=logspace(1,log10(fmax),N);

%fs = [10, 15, 20, 32, 40, 63, 100, 125, 150, 200, 500, 1000, 2000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 11050]
wn = fs*pi/max(fs);

wt=0*wn;
    % and adjust those weights to ones fixed above
    for J=1:length(freq_weights)
        [~,I] = min(abs(fs-freq_weights(J)));
        wt(I) = weights(J);
    end

%wt = 1./(wn).^2
%wt = []
% A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
% Hz and 738 Hz according to S1.42.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;

NUMa = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENa = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENa = conv(conv(DENa,[1 2*pi*f3]),[1 2*pi*f2]); 



% generate the freq response of the analog filter
ha=freqs(NUMa,DENa,2*pi*fs);
absha=abs(ha);

% C weighting filter has 2 poles at 20 Hz and 12.2 kHz according ANSI S1.42
f1 = 20.598997; 
f4 = 12194.217;
C1000 = 0.0619;

NUMc = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
DENc = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 

% generate the freq response of the analog filter
hc=freqs(NUMc,DENc,2*pi*fs);



[ba,aa] = invfreqz(ha,wn,Nb, Na, wt)
[ba,aa] = invfreqz(ha,wn,Nb, Na, wt, 100, 0.00001)
% [ba,aa] = bilinear(NUMa,DENa,FS, 100)

% generate the freq response of the digital filter at fs
[hd,fd]=freqz(ba,aa,fs,FS);

ampc = 20*log10(abs(ha));
ampd = 20*log10(abs(hd));


ampdiff = ampd - ampc;

subplot(2,1,1)
semilogx(fs, ampc,'b',fd,ampd,'r')

axis([min(fs),max(fs),-60,+10])
legend('Analog','Digital','location','south')

subplot(2,1,2)

semilogx(fs,ampdiff)
axis([min(fs),max(fs),-10,+10])
grid
%axis([20,fmax,-60,+10])
%legend('Analog','Digital','location','south')

