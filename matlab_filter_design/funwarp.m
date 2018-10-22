function w = invfwarp(ww, wc);
    % frequency unwarping as described in many papers including
    % B. Bank, “Warped IIR Filter Design with Custom Warping Profiles and 
    % Application to Room Response Modeling and Equalization,” presented at 
    % the Audio Engineering Society Convention 130, 2011.
    % 
    % This is the inverse of fwarp
    % 
    % Inputs:
    % ww = normalized warped frequencies i.e. 2*pi*f/Fs or 0 ... pi
    % wc = normalized frequency to split the mapping
    %
    % Outputs:
    % w = unwarped frequencies
    

%     
%     a=pi/(wc*(1+log(pi/wc)));
%     b = exp(1)/wc;;
%     logbpi = log(b*pi);
%     
%     w = zeros(size(ww));
% 
%     
%     % get the warped corner frequency for transition
%     wwc = fwarp(wc, wc);
%     
%     % unwarp the low frequencies (ww < wwc)
%     I = find(ww < wwc);
%     w(I) = ww(I)/a;
%     
%     % unwarp the high frequencies (ww >= wwc)
%     I = find(ww >= wwc);
%     w(I) = exp(ww(I)/pi*logbpi)/b;
    
    L = wc;
    w = atan2( (1-L.^2)*sin(ww),((1+L^2)*cos(ww) + 2*L));
end

