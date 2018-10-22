function ww = fwarp(w, wc)
    % frequency warping as described in many papers including
    % B. Bank, “Warped IIR Filter Design with Custom Warping Profiles and 
    % Application to Room Response Modeling and Equalization,” presented at 
    % the Audio Engineering Society Convention 130, 2011.
    % 
    % We do a scaled linear mapping for w<wc and a log mapping for w>wc
    % 
    % Inputs:
    % w = normalized frequencies i.e. 2*pi*f/Fs or 0 ... pi
    % wc = normalized frequency to split the mapping
    %
    % Outputs:
    % ww = warped frequencies
%     
%     a=pi/(wc*(1+log(pi/wc)));
%     b = exp(1)/wc;
%     logbpi = log(b*pi);
%     
%     ww=zeros(size(w));
%     I = find (w < wc);
%     ww(I) = a*w(I);
%    
%     I = find (w >= wc);
%     ww(I)=pi*log(b*w(I))/logbpi;
%     
    L = wc;
    
    ww = atan2(  (1-L.^2)*sin(w),((1+L.^2)*cos(w)-2*L));

end % function ww

