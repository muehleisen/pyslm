
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

Fs=44100;

f1 = 20.598997; 
f4 = 12194.217;
C1000 = 0.0619;
pi = 3.14159265358979;
NUMs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 


%for 44100, Nlist=[5307 5789 6810 6935 9335 10575 14952];
% for 48000 Nlist=14860

errmin=1;
errba=1;

Nmins=[1];

% [ 1045, 1831, 3205
Nmin=3000;
Nspace=1;
Nmax=10000;

Nlist = [1831, Nmin:Nspace:Nmax];
% Nlist=1831
for N=Nlist
    [N,errba,Nmin,errmin]
    %fs=linspace(10,Fs/2,N);
    if Fs>40000
        fmax=20000;
    else
        fmax=Fs/2;
    end
    fs=linspace(0,fmax,N);
    %fs=logspace(1,log10(fmax),N);
    [hs]=freqs(NUMs,DENs,2*pi*fs);


    % Use the bilinear transformation to get the digital filter for comparison. 
    [numd,dend] = bilinear(NUMs,DENs,Fs);
    [hb,fb]=freqz(numd,dend,fs,Fs);

    wt=1./(fs).^2;
    wt(1)=0;
    wtimpt=1;
    [~,I]=min(abs(fs-15));
    wt(I)=wtimpt*600;
    [~,I]=min(abs(fs-20));
    wt(I)=wtimpt*500;
    [~,I]=min(abs(fs-25));
    wt(I)=wtimpt*400;    
    [~,I]=min(abs(fs-32));
    wt(I)=wtimpt*300;
    [~,I]=min(abs(fs-50));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-63));
    wt(I)=wtimpt*200;
    [~,I]=min(abs(fs-125));
    wt(I)=wtimpt*100;
    [~,I]=min(abs(fs-250));
    wt(I)=wtimpt*100;
    [~,I]=min(abs(fs-500));
    wt(I)=wtimpt*100;
    [~,I]=min(abs(fs-1000));
    wt(I)=wtimpt*100;
    [~,I]=min(abs(fs-2000));
    wt(I)=wtimpt*10;
    [~,I]=min(abs(fs-4000));
    wt(I)=wtimpt*1;
    [~,I]=min(abs(fs-8000));
    wt(I)=wtimpt*0.1;
     [~,I]=min(abs(fs-16000));
    wt(I)=wtimpt*0.001;
    [~,I]=find(abs(fs-17000)<50);
    wt(I)=wtimpt*0.0001;
    [~,I]=find(abs(fs-18000)<50);
    wt(I)=wtimpt*0.0001;
      [~,I]=find(abs(fs-19000)<50);
    wt(I)=wtimpt*0.0001;  
    [~,I]=find(abs(fs-20000)<50);
    wt(I)=wtimpt*0.00001;
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
        
        Nmin=N;
        Nmins=[Nmins;Nmin];
        errmin=errba;
        wtmin=wt;
        
        subplot(2,1,1)
        semilogx(fs,20*log10(abs(hb)),'g',fb,20*log10(abs(hs)),'r',fd4,20*log10(abs(hd4)),'k')
     	axis([20,fmax,-70,+5])
        legend('bilinear','analog','invfreqz','location', 'southwest')
        title(sprintf("N=%d",Nmin))
        
        subplot(2,1,2)
        semilogx(fs,unwrap(angle(hb)-2*pi),'g',fb,unwrap(angle(hs)),'r',fd4,unwrap(angle(hd4)),'k')
        axis([20,fmax,-30,+5])
%         subplot(2,2,3)
%         plot(fd4,fulldiff)
%         axis([5,40,-5,20])
%         grid
% 
%         subplot(2,2,4)
%         plot(fd4,fulldiff)
%         axis([10000,Fs/2,-1,2])
%         grid

        pause(1)
    end
end

Nmin
disp(sprintf('a_Cwt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,ad4))
disp(sprintf('b_Cwt_%2.0fk=[%#1.15g,%#1.15g,%#1.15g,%#1.15g,%#1.15g]',Fs/1000,bd4))

