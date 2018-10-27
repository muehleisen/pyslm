function [bd, ad] = invfreqzw(Fs, weighting, Nab, Npts, lambda)
% inputs:
% Fs = sample frequency (default = 48000)
% weighting = "A" for A-weighting, "B" for B-weighting, and "C" for
% C-weighting (default = "A
% Nab = Na = Nb = number of zeros and poles (default = 5)
% Npts = number of frequencies to use in error calculation
%  
% 

if nargin<5
    lambda = 0.995
end

if nargin < 4
    Npts = 2048
end

if nargin < 3
    Nab = 4
end

if nargin < 2
    weighting = "A"
end

if nargin <1
    Fs = 48000
end

Na = Nab;
Nb = Nab;
weighting = upper(weighting);

fmatch = 1000; % set frequency for matching filter gains
% set frequency for splitting resolution of freq vector
% 100 Hz was found to work by trail and error
fsplit = 100; 

if weighting == "C"
    % C weighting filter has poles at 20.6 and 12194 Hz
    f1 = 20.598997; 
    f4 = 12194.217;
    C1000 = 0.0619;
    bs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
    as = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
      
elseif weighting =="B"
    % B weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole 
    % 158 Hz according to S1.42.
    f1 = 20.598997; 
    f2 = 158.48932;
    f4 = 12194.217;
    B1000 = 0.17480;
    Ks = (2*pi*f4)^2*10^(B1000/20);
    Zs = [0 0 0].';
    Ps = 2*pi*[f4, f4, f2, f1, f1].';
    [bs, as] = zp2tf(Zs, Ps, Ks);

else % assume A weighting filter if none given 
    % A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
    % Hz and 738 Hz according to S1.42.
    f1 = 20.598997; 
    f2 = 107.65265;
    f3 = 737.86223;
    f4 = 12194.217;
    A1000 = 1.9997;
    bs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
    as = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
    as = conv(conv(as,[1 2*pi*f3]),[1 2*pi*f2]); 
end

% define tolerance for Type 0 and 1 filters from S1.43
ftol= [10, 12.5, 16, 20, 25, 31.5,  40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];
up1 = [ 4, 3.5,  3, 2.5,  2,  1.5, 1.5,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,    1,    1,    1,  1.5,  1.5,  1.5,     2,     3,   3.0, 3 ];
low1 =[ 4, 3.5,  3, 2.5,  2,  1.5, 1.5,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,    1,    1,    1,  1.5,    2,    3,     4,     6,   inf, inf ];
up0  = [ 2,  2,  2,   2,1.5,     1,  1,  1,  1,  1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,    1,    1,    1,     2,     2,     2, 2];  
low0 = [ 5,  3,  3,   2,1.5,     1,  1,  1,  1,  1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,    1,  1.5,    2,     3,     3,     3, 3];  

% generate a frequency vector for use in the curve fit make it evenly 
% spaced at low frequencies and high frequencies but
% different spacing in the two regions.  
f1 = linspace(0, fsplit, Npts);
f2 = linspace(fsplit, Fs/2, Npts +1);
f = [f1, f2(2:end)];

wn=2*pi*f/Fs;

% generate the freq response of the analog filter for fitting
hs=freqs(bs,as,wn*Fs);
dbs = 10*log10(abs(hs));

% generate a frequency warping for the ditigal filter.
% frequency warping as described in many papers.  One of the first is in
% A. Oppenheim, D. Johnson, and K. Steiglitz, “Computation of spectra with 
% unequal resolution using the fast Fourier transform,” Proceedings of the 
% IEEE, vol. 59, no. 2, pp. 299–301, Feb. 1971.
wwn=fwarp(wn, lambda);

% convert to a minimum phase version for easier fitting
% use minphasen from B. Bank with matlab code at
% http://home.mit.bme.hu/~bank/parfilt/
% but first, make the wn=0 value of hs be nonzero so the hilbert xform exists
% by making it 1/100 of the value of the 2nd point
hs(1) = hs(2)/100;
hmp = minphasen(abs(hs),wwn);

% fit the ditial filter at the warped frequencies
[bdw, adw] = invfreqz(hmp, wwn, Nb, Na);

% unwarp the coefficients to get back to normal frequencies
[bd, ad] = unwarp(bdw, adw, lambda);

% compute the frequency response of the digital filter over unwarped freq.
hd = freqz(bd, ad, wn);

% find the index of wn closest to 1000 Hz and match |Hd| and |Hs| at that frequency 
[~, Imatch] = min(abs(wn-2*pi*fmatch/Fs));
gainmatch = abs(hs(Imatch))./abs(hd(Imatch));
hd = hd*gainmatch;
bd = bd*gainmatch;

dbs = 10*log10(abs(hs));
dbd = 10*log10(abs(hd));

% Use the bilinear transformation to get the BLT digital filter for comparison. 
[bblt,ablt] = bilinear(bs,as,Fs);
hblt=freqz(bblt,ablt,f,Fs);
dbblt = 10*log10(abs(hblt));

% compute db differences between digital and analog filters
rdiff = dbd(:) - dbs(:);
bdiff = dbblt(:) - dbs(:);

% plot freq response, difference, and pole-zero plot of final filter
subplot(3,1,1)
semilogx(f,dbs , 'k',f,dbd, 'r', f, dbblt, 'g')
axis([1,30000,-70,+3])

subplot(3,1,2)
semilogx(f, rdiff, 'r', f, bdiff, 'g', ftol, up0, 'k--', ftol, -low0, 'k--')
grid
legend('invfreqzw','bilinear','location', 'north')
axis([1,30000,-5,+5])

subplot(3,1,3)
zplane(bd,ad)

[filepath,name,ext] = fileparts(mfilename('fullpath'));
disp(sprintf('# from %s(%d, \"%s\", %d, %d, %1.3f)', name, Fs, weighting, Nab, Npts, lambda));

astring = ['a = ['];
for I=1:(Na+1)
    astring = strcat(astring, sprintf(' %#1.15g,',ad(I)));
end
astring(end)=']';

bstring = ['b = ['];
for I=1:(Nb+1)
    bstring = strcat(bstring, sprintf(' %#1.15g,',bd(I)));
end
bstring(end)=']';

disp(astring)
disp(bstring)

% end
