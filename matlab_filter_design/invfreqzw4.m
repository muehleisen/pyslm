function [bdmin, admin] = invfreqzw(Fs, weighting, Nab, Nerr, lambda)
% inputs:
% Fs = sample frequency (default = 48000)
% weighting = "A" for A-weighting, "B" for B-weighting, and "C" for
% C-weighting (default = "A
% Nab = Na = Nb = number of zeros and poles (default = 5)
% Nerr = number of frequencies to use in error calculation
%  
% 
% Fs =44100;
% Nerr =128;
% Nab = 6;
% weighting = "A";
Nloops=1

if nargin<5
    lambda = 0.995
end

if nargin < 4
    Nerr = 256
end

if nargin < 3
    Nab = 4
end

if nargin < 2
    weighting = "A"
end

if nargin <1
    Fs = 96000
end

Na = Nab;
Nb = Nab;
weighting = upper(weighting);

fmatch = 1000;

if weighting == "C"
    % C weighting filter has poles at 20.6 and 12194 Hz
    f1 = 20.598997; 
    f4 = 12194.217;
    C1000 = 0.0619;
    bs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
    as = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
else
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
    
ferr = logspace(log10(10),log10(20000),Nerr);
werr = 2*pi*ferr;

hserr = freqs(bs, as, werr);
dbserr = 10*log10(abs(hserr));

% f = logspace(0, log10(Fs/2), Nerr);
% f(1) = 0;
% f(end)=Fs/2;

f1 = linspace(0, 100, Nerr+1);
f2 = linspace(100, Fs/2, Nerr);

f = [f1, f2(2:end)];


wn=2*pi*f/Fs;

% generate the freq response of the analog filter for fitting
hs=freqs(bs,as,wn*Fs);
dbs = 10*log10(abs(hs));


% frequency warping as described in many papers including
% B. Bank, “Warped IIR Filter Design with Custom Warping Profiles and 
% Application to Room Response Modeling and Equalization,” presented at 
% the Audio Engineering Society Convention 130, 2011.
% 
% At low freq we use a standard linear mapping of high resolution and
% at high frequencies we use a logarithmic mapping

wwn=fwarp(wn, lambda);
wwnerr = fwarp(werr/Fs, lambda);

wnp = funwarp(wwn,lambda);

% convert to a minimum phase version for easier fitting
% use minphasen from B. Bank as mentioned above and matlab code at
% http://home.mit.bme.hu/~bank/parfilt/

% but first, make the zero value of hs is nonzero so the hilbert xform exists
% by taking 1/100 of the value of the 2nd point
hs(1) = hs(2)/100;
hmp = minphasen(abs(hs),wwn);
hmp = hmp(:).';   %make sure hmp is a row like all the other h

% dbmp = 10*log10(abs(hmp));

% W = ones(size(wwn));
% [~,I10]  =min(abs(wwn*Fs/2/pi - 10)); % find index of wn closest to 10 Hz
% Wmask = ones(size(wwn));
% Wmask(1:I10) = 0;
% errmin = 1000;

[bdwmin, adwmin] = invfreqz(hmp, wwn, Nb, Na);

[bdmin, admin] = unwarp(bdwmin, adwmin, lambda);
hd = freqz(bdmin, admin, wn);


% find the index of wn closest to 1000 Hz to match Hn at 1000 Hz later
[~,Imatch] = min(abs(wn-2*pi*fmatch/Fs));

gainmatch = abs(hs(Imatch))./abs(hd(Imatch));
hd = hd*gainmatch;

bdmin = bdmin*gainmatch;

dbs = 10*log10(abs(hs));
dbd = 10*log10(abs(hd));

% Use the bilinear transformation to get the BLT digital filter for comparison. 
[bblt,ablt] = bilinear(bs,as,Fs);

% get the freq response of the bilinear filter
hblt=freqz(bblt,ablt,f,Fs);
dbblt = 10*log10(abs(hblt));

% tolerance for Type 1 filters from S1.43
ftol= [10, 12.5, 16, 20, 25, 31.5,  40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];
up1 = [ 4, 3.5,  3, 2.5,  2,  1.5, 1.5,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,    1,    1,    1,  1.5,  1.5,  1.5,     2,     3,   3.0, 3 ];
low1 =[ 4, 3.5,  3, 2.5,  2,  1.5, 1.5,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,    1,    1,    1,  1.5,    2,    3,     4,     6,   inf, inf ];
up0  = [ 2,  2,  2,   2,1.5,     1,  1,  1,  1,  1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,    1,    1,    1,     2,     2,     2, 2];  
low0 = [ 5,  3,  3,   2,1.5,     1,  1,  1,  1,  1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,    1,  1.5,    2,     3,     3,     3, 3];  

rdiff = dbd(:) - dbs(:);
bdiff = dbblt(:) - dbs(:);

subplot(3,1,1)
semilogx(f,dbs , 'k',f,dbd, 'r', f, dbblt, 'g')
axis([1,30000,-70,+3])

subplot(3,1,2)
semilogx(f, rdiff, 'r', f, bdiff, 'g', ftol, up0, 'k--', ftol, -low0, 'k--')
grid
legend('invfreqzw','bilinear','location', 'north')
axis([1,30000,-5,+5])



subplot(3,1,3)
zplane(bdmin,admin)

% disp(sprintf('# %s-Weighting',weighting));
[zd, pd, kd]=tf2zp(bdmin, admin);
[zp, pb, kb] = tf2zp(bblt, ablt);




[filepath,name,ext] = fileparts(mfilename('fullpath'));
% disp(sprintf('# From %s.m: Fs = %d, Na = %d, Nb = %d, lambda = %3.3f, Nwt = %3.1f, Nerr=%d', name, Fs, Na, Nb, lambdamin, Nwtmin, Nerr))

disp(sprintf('# from %s(%d, \"%s\", %d, %d, %1.3f)', name, Fs, weighting, Nab, Nerr, lambda));
% disp(sprintf('# solution has lambda = %3.3f',lambdamin, Nwtmin));

astring = ['a = ['];
for I=1:(Na+1)
    astring = strcat(astring, sprintf(' %#1.15g,',admin(I)));
end
astring(end)=']';

bstring = ['b = ['];
for I=1:(Nb+1)
    bstring = strcat(bstring, sprintf(' %#1.15g,',bdmin(I)));
end
bstring(end)=']';
disp(astring)
disp(bstring)

% end






