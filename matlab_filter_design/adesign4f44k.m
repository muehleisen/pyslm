
Fs=44100;
Norder = 7;

% A weighting filter has 2 poles at 20 Hz and 12.2 kHz and one pole at 108
% Hz and 738 Hz according to S1.42.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;

NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]); 



fs=logspace(log10(20),log10(20000), 1000)

% generate the freq response of the analog filter
hs=freqs(NUMs,DENs,2*pi*fs);

% Use the bilinear transformation to get the digital filter for comparison. 
[bbi,abi] = bilinear(NUMs,DENs,Fs);
% get the freq response of the bilinear filter
[hb,fb]=freqz(bbi,abi,fs,Fs);

%for 44100, Nlist=[5307 5789 6810 6935 9335 10575 14952];
% for 48000 Nlist=14860

errmin=10;
errba=10;

Nmins=[1];

Nmin=32;
Nspace=1;
Nmax=200;
% Nlist = 169 363 897
Nlist = [Nmin:Nspace:Nmax];
%Nlist = [897, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]

freq_weights = [15, 20, 25, 32, 50, 63,  125, 250, 500, 1000, 2000, 4000, 8000, 10000, 16000, 18000, 19000, 20000] 
% works well at lf with  N = 197,  
% weights = [5000, 5000, 500, 200, 200, 200, 200, 200, 200, 200, 200, 200, 20, .01, .001, 0.0001]
% works better at lf with N=139
% weights = [5000, 5000, 5000, 2000, 200, 200, 200, 200, 200, 200, 200, 200, 20, .001, .0001, 0.001]

% N=150 with weights below
weights = [6000, 5000, 4000, 3000, 200, 200, 100, 100, 100, 100, 10, 1, 1, .01, 0.001, 0.001, 0.001, 0.001]


%size(freq_weights)
%size(weights)




for N=Nlist
    [N,errba,Nmin,errmin]
    %fs=linspace(10,Fs/2,N);
    if Fs>40000
        fmax=20000;
    else
        fmax=Fs/2;
    end

    % generate a new freq vector of length N
    fs2=logspace(0.5,log10(Fs/2),N);

    wt=1./(fs2).^2;
    
    wt(1)=0;
    %wt=zeros(size(fs));
    
    % find the frequencies closest to the weight freq
    % and adjust those weights to ones fixed above
    for J=1:length(freq_weights)
        [~,I] = min(abs(fs2-freq_weights(J)));
        wt(I) = weights(J);
    end
    
    % generate the freq response of the analog filter at fs2
    hs2=freqs(NUMs,DENs,2*pi*fs2);
    % use invfreqz to get filter coefficients
    [bd4,ad4]=invfreqz(hs2,pi*fs2/max(fs2),Norder,Norder,wt,10,0.001);
    %[bd4,ad4]=invfreqz(hs2,pi*fs2/max(fs2),Norder,Norder,wt);
    % generate the freq response of the digital filter at fs
    [hd4,fd4]=freqz(bd4,ad4,fs,Fs);
     
    
    
    
    [~,Imin]=min(abs(f2-20));  % find I closest to 10
    [~,Imax]=min(abs(f2-20000)); % find I closest to 20k
    % ahd=abs(hd4(Imin:Imax));
    % ahs=abs(hs(Imin:Imax));
    ahd=abs(hd4);
    ahs=abs(hs);
      
    
    diff=20*log10(ahd./ahs);
    % errba=norm(diff)/1000;
    errba = max(abs(diff));
    if errba<errmin
        
        fulldiff=20*log10(abs(hd4)./abs(hs));
        bidiff = 20*log10(abs(hb)./abs(hs));
        Nmin=N;
        Nmins=[Nmins;Nmin];
        errmin=errba;
        wtmin=wt;
 
        
        subplot(2,1,1)
        semilogx(fs,20*log10(abs(hb)),'g',fb,20*log10(abs(hs)),'r',fd4,20*log10(abs(hd4)),'k')
     	axis([20,fmax,-60,+10])
        legend('bilinear','analog','invfreqz','location','south')
        title(sprintf("N=%d",Nmin))
        
        subplot(2,1,2)
     	semilogx(fs,bidiff,'r', fs, fulldiff, 'k')
        grid
        
        axis([20,fmax,-2,+5])
        legend('bilinear','invfreqz','location', 'north')
        title(sprintf("N=%d",Nmin))
        

        pause(1)
    end
end
Nmin
Fs
disp('From freqz')
disp(sprintf('a = [%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',ad4))
disp(sprintf('b = [%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',bd4))
disp('From bilinear')
disp(sprintf('a_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,abi))
disp(sprintf('b_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,bbi))

