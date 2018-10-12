
Fs=44100;
Nb = 4;
Na = 4;

f=logspace(log10(20), log10(20000), 512)
% 
% A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
% Hz and 738 Hz according to S1.42.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
% 
% Generate analog filter by using Numerator and Denominator Polynomials
NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]); 
% generate the freq response of the analog filter
% hsa=freqs(NUMs,DENs,2*pi*fs)

% generate the filter using Pole, Zero, Gain formulation
Ka = (2*pi*f4).^2 * (10.^(A1000/20))
Za = [0, 0, 0, 0]'
Pa = 2*pi*[f1, f4, f3, f2, f1, f4]'

[ba,aa] = zp2tf(Za,Pa,Ka)
hsa = freqs(ba, aa, 2*pi*f);
hsa = freqs(NUMs, DENs, 2*pi*f);
% Use the bilinear transformation to get the digital filter coefficients 
% for comparison. 
%[bb,ab] = bilinear(NUMs,DENs,Fs);
[Zbd, Pbd, Kbd] = bilinear(Za, Pa, Ka, Fs);
[bb, ab] = zp2tf(Zbd, Pbd, Kbd)
% get the freq response of the bilinear filter 
[hb,fb]=freqz(bb,ab,f,Fs);


% semilogx(f, 20*log10(abs(hb)), f, 20*log10(abs(hsa)))
% legend('bilinear', 'analog')


%for 44100, Nlist=[5307 5789 6810 6935 9335 10575 14952];
% for 48000 Nlist=14860

errmin=100;
errba=10;

Nmins=[1];

Nmin=32;
Nspace=1;
Nmax=500;
% Nlist = 169 363 897
Nlist = [Nmin:Nspace:Nmax];
%Nlist = [897, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]

freq_weights = [15, 20, 25, 32, 50, 63,  125, 250, 500, 1000, 2000, 4000, 8000, 10000, 16000, 18000, 19000, 20000] 
% works well at lf with  N = 197,  
% weights = [5000, 5000, 500, 200, 200, 200, 200, 200, 200, 200, 200, 200, 20, .01, .001, 0.0001]
% works better at lf with N=139
% weights = [5000, 5000, 5000, 2000, 200, 200, 200, 200, 200, 200, 200, 200, 20, .001, .0001, 0.001]

% N=150 with weights below
weights = [10000, 10000, 10000, 300, 2000, 200, 100, 100, 100, 100, 10, 10, 10, 1, 0.01, 0.01, 0.01, .01];
weights = [6000, 5000, 4000, 3000, 200, 200, 100, 100, 100, 100, 10, 1, 1, .01, 0.001, 0.001, 0.001, 0.001]

%size(freq_weights)
%size(weights)


freq_weights = [15, 20, 25, 32, 50, 63,  125, 250, 500, 1000, 2000, 4000, 8000, 10000, 16000, 18000, 19000, 20000] 
weights = [6000, 5000, 4000, 3000, 200, 200, 100, 100, 100, 100, 10, 1, 1, .01, 0.001, 0.001, 0.001, 0.001]

for N=Nlist
    [N,errba,Nmin,errmin]
    %fs=linspace(10,Fs/2,N);
    if Fs>40000
        fmax=20000;
    else
        fmax=Fs/2;
    end

    % generate a new freq vector of length N
    f2=logspace(1,log10(fmax),N);
    
    f2 = logspace(0.5,log10(Fs/2),N);
    f2(end) = Fs/2;

    wt=1./(f2).^2;
    
    wt(1)=0;
    %wt=zeros(size(fs));
    
    % find the frequencies closest to the weight freq
    % and adjust those weights to ones fixed above
    for J=1:length(freq_weights)
        [~,I] = min(abs(f2-freq_weights(J)));
        wt(I) = weights(J);
    end
    
    % generate the freq response of the analog filter at f2
    
    hs2=freqs(NUMs,DENs,2*pi*f2);
    %hs2=freqs(ba, aa, 2*pi*f2);
    
    
    % use invfreqz to get filter coefficients
    [bd2,ad2]=invfreqz(hs2,2*pi*f2/Fs,Nb,Na,wt,10,0.01);
    % [bd2,ad2]=invfreqz(hs2, pi*f2/max(f2), Nb,Na,wt);
    
    % generate the freq response of the digital filter at fs
    [hd2,fd2]=freqz(bd2,ad2,f,Fs);
     
   
    [~,Imin]=min(abs(f2-20));  % find I closest to 10
    [~,Imax]=min(abs(f2-20000)); % find I closest to 20k
    % ahd=abs(hd4(Imin:Imax));
    % ahs=abs(hs(Imin:Imax));
    ahd=abs(hd2);
    ahs=abs(hsa);
      
    
    diff=20*log10(ahd./ahs);
    % errba=norm(diff)/1000;
    errba = max(abs(diff));
    if errba<errmin
        
        [hdf,fdf]=freqz(bd2,ad2,f,Fs);
        
        fulldiff=20*log10(abs(hdf)./abs(hsa));
        bidiff = 20*log10(abs(hb)./abs(hsa));
        Nmin=N;
        Nmins=[Nmins;Nmin];
        errmin=errba;
        wtmin=wt;
 
        
        subplot(2,1,1)
        semilogx(f,20*log10(abs(hb)),'g',f,20*log10(abs(hsa)),'r',f,20*log10(abs(hdf)),'k')
     	axis([20,max(f),-60,+10])
        legend('bilinear','analog','invfreqz','location','south')
        title(sprintf("N=%d",Nmin))
        
        subplot(2,1,2)
     	semilogx(f,bidiff,'r', f, fulldiff, 'k')
        grid
        
        axis([20,max(f),-2,+5])
        legend('bilinear','invfreqz','location', 'north')
        title(sprintf("N=%d",Nmin))
        

        pause(1)
    end
end
Nmin
disp('From freqz')
disp(sprintf('a_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,ad2))
disp(sprintf('b_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g,]',Fs/1000,bd2))
disp('From bilinear')
disp(sprintf('a_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,ab))
disp(sprintf('b_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,bb))

