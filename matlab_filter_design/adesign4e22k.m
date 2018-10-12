
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

Fs=22050;

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


%for 44100, Nlist=[5307 5789 6810 6935 9335 10575 14952];
% for 48000 Nlist=14860

errmin=10;
errba=1;

Nmins=[1];

Nmin=10000;
Nspace=1;
Nmax=15000;
% Nlist =  216, 897
Nlist = [897, Nmin:Nspace:Nmax];

for N=Nlist
    [N,errba,Nmin,errmin]
    %fs=linspace(10,Fs/2,N);
    if Fs>40000
        fmax=20000;
    else
        fmax=Fs/2;
    end
    %fs=linspace(10,fmax,N);
    fs=logspace(0,log10(Fs/2),N);
    [hs]=freqs(NUMs,DENs,2*pi*fs);


    % Use the bilinear transformation to get the digital filter for comparison. 
    [bbi,abi] = bilinear(NUMs,DENs,Fs);
    [hb,fb]=freqz(bbi,abi,fs,Fs);

    wt=1./(fs).^2;
    wt(1)=0;
    wtimpt=1;
    [~,I]=min(abs(fs-15));
    wt(I)=wtimpt*5000;
    [~,I]=min(abs(fs-20));
    wt(I)=wtimpt*5000;
    [~,I]=min(abs(fs-25));
    wt(I)=wtimpt*500;    
    [~,I]=min(abs(fs-32));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-50));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-63));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-125));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-250));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-500));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-1000));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-2000));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-4000));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-6000));
    wt(I)=wtimpt*20;
    [~,I]=min(abs(fs-10000));
    wt(I)=wtimpt*0.01;
      [~,I]=min(abs(fs-16000));
    wt(I)=wtimpt*0.001;
    [~,I]=min(abs(fs-18000));
     wt(I)=wtimpt*0.0001;
%       [~,I]=find(abs(fs-19000)<50);
%     wt(I)=wtimpt*0.001;  
%     [~,I]=find(abs(fs-20000)<50);
%     wt(I)=wtimpt*0.001;
    Ilast=I;

    [bd4,ad4]=invfreqz(hs,pi*fs/max(fs),4,4,wt,100,0.00001);
    [hd4,fd4]=freqz(bd4,ad4,fs,Fs);
    
    
    [~,I10]=min(abs(fs-10));  % find I closest to 10
    [~,I20k]=min(abs(fs-20000)); % find I closest to 20k
    ahd=abs(hd4(I10:I20k));
    ahs=abs(hs(I10:I20k));
    diff=20*log10(ahd./ahs);
    errba=norm(diff)/N;
    
    
    
    
    if errba<errmin
        
        
        fulldiff=20*log10(abs(hd4)./abs(hs));
        bidiff = 20*log10(abs(hb)./abs(hs));
        Nmin=N;
        Nmins=[Nmins;Nmin];
        errmin=errba;
        wtmin=wt;
 
        
        subplot(2,1,1)
        semilogx(fs,20*log10(abs(hb)),'g',fb,20*log10(abs(hs)),'r',fd4,20*log10(abs(hd4)),'k')
     	axis([20,fmax,-70,+5])
        legend('bilinear','analog','invfreqz','location','south')
        title(sprintf("N=%d",Nmin))
        
        subplot(2,1,2)
     	semilogx(fs,bidiff,'r', fs, fulldiff, 'k')
       
        axis([20,fmax,-2,+2])
        legend('bilinear','invfreqz','location', 'south')
        title(sprintf("N=%d",Nmin))
        

        pause(1)
    end
end
Nmin
disp('From freqz')
disp(sprintf('a_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,ad4))
disp(sprintf('b_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,bd4))
disp('From bilinear')
disp(sprintf('a_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,abi))
disp(sprintf('b_Awt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,bbi))

