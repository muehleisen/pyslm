
% ADSGN  Design of a A-weighting filter.
%    [B,A] = ADSGN(Fs) designs a digital A-weighting filter for 
%    sampling frequency Fs. Usage: Y = FILTER(B,A,X). 
%    Warning: Fs should normally be higher than 20 kHz. For example, 
%    Fs = 48000 yields a class 1-compliant filter.
%
%    This routine uses invfreqz which does a least squares curve fit to the
%    frequency response instead of the more usual bilinear transform
%    The BLT has problems with the low and high frequency poles at 20 Hz
%    and 12200 Hz.
%
%    Requires the Signal Processing Toolbox. 
%
%    See also ASPEC, CDSGN, CSPEC. 



% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 20, 1997, 10:00am.

% References: 
%    [1] IEC/CD 1672: Electroacoustics-Sound Level Meters, Nov. 1996. 

% set sampling frequency at Fs = 44100 Hz by default


% Definition of analog A-weighting filter according to IEC/CD 1672.
% Analog A-weighting filter coefficients according to ANSI
N=8192;
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;

fawt=[f4 f4 f3 f2 f1 f1]




%K=4*pi^2*12200^2;
A1000 = 1.9997;
K=(2*pi*f4)^2*(10^(A1000/20));

fs=linspace(0,Fs/2,N);
Z=[0 0 0 0]';
P=2*pi*[f4 f4 f3 f2 f1 f1]';
[Ba,Aa]=zp2tf(Z,P,K);
[ha,wa]=freqs(Ba,Aa,2*pi*fs);
fa=wa/2/pi;



% % Zl=1';
% Pl=2*pi*[f3 f2 f1 f1];
% A1000 = 1.9997;
% Kl=(2*pi*f4).^2*(10^(A1000/20));
% [Bl,Al]=zp2tf(Zl,Pl,Kl)
% [hl,wl]=freqs(Bl,Al,1024);
% fl=wl/2/pi;

%pi = 3.14159265358979;
NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);



[hs]=freqs(NUMs,DENs,2*pi*fs);
%[ha,wa]=freqs(Ba,Aa,2*pi*fs);
%fs=ws/2/pi;

%fs=linspace(0,Fs/2,1024);
% Use the bilinear transformation to get the digital filter. 
[numd,dend] = bilinear(NUMs,DENs,Fs);
[hb,fb]=freqz(numd,dend,fs,Fs);

% Use the bilinear transformation to get the digital filter. 
% [n2,d2] = bilinear(NUMs,DENs,Fs,2000)
% [hd,fd]=freqz(n2,d2,1024,Fs);
%[bd,ad]=yulewalk(4,fs/max(fs),abs(hs));

%wn=pi*linspace(0,1,1024);
wt=ones(size(fs));


I=find((fs<1000));
wt(I)=15*I./I;
I=find((fs<300));
wt(I)=200*I./I;
I=find((fs<100));
wt(I)=200*I./I;
I=find((fs<40));
wt(I)=500*I./I;b

I=find((fs>1000));
wt(I)=.1*I./I;
I=find((fs>1000));
wt(I)=0.0*I./I;


%[hd1,wd1]=freqs(NUMs,DENs,fs,Fs);
[bd,ad]=invfreqz(hs,pi*fs/max(fs),4,4,wt,30,0.00001,'trace')


[bd,ad]=bilinear(NUMa,DENa,FS,10000)

[hd,fd]=freqz(bd,ad,fs,Fs);


semilogx(fa,20*log10(abs(ha)),'r',fb,20*log10(abs(hb)),'g',fd,20*log10(abs(hd)),'k')
axis([10,20000,-70,5])
