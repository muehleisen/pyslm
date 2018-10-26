function [b,a] = invlogz(D, wn, Nb, Na, Fs )
%invlogz Find IIR coefficients b and a that minimize error to ln(D) 
% 
%
% This function finds the [b, a] that best match a transfer function that
% is minimum phase and has the amplitude defined in D
%
% It uses invfreqz and iterates using a weighting function that is related
% to the error in the fit as defined in the paper
%
% [1]T. Kobayashi and S. Imai, “Design of IIR digital filters with arbitrary 
% log magnitude function by WLS techniques,” IEEE Transactions on Acoustics, 
% Speech, and Signal Processing, vol. 38, no. 2, pp. 247–252, Feb. 1990.
%
% This routine does not do the minimization as described in that paper, 
% just updates the weighting function using the method described in the 
% paper

% first create a minimum phase version of the filter using minphasen
% use magnitude of D in case someone inputs a complex number and not a
% magnitude

% before we do that, make sure D(1) is not zero by interpolating D(2) and D(3)
% or D(2)/100 if D(1) <=0
% use c = polyfit( [x2, x3],[y2, y3], 1) where c(2)is b, the y intercept.
if D(1) == 0
    % fit a line to the 2nd and 3rd pts of D to get slope (C(1) and
    % y intercept C(2) )
    C = polyfit( wn(2:3), D(2:3), 1);
    
    % if the intercept > 0 then set D(1) to be intercept, 
    if C(2) > 0
        D(1) = c(2);
    else % otherwise intercept <=0
        D(1) = D(2)/100;
    end % if c(2) > 0
end
Dmp = minphasen(abs(D));

% get starting attempt at b and a using simple invfreqz


Nw = length(wn)

W = ones(size(wn));
[~,I10]  =min(abs(wn*Fs/2/pi - 10)); % find index of wn closest to 10 Hz
Wmask = ones(size(wn));
Wmask(1:I10) = 0;


errmin = 100
for J=1:100
    [b, a] = invfreqz(Dmp, wn, Nb, Na, W);
    Hd = freqz(b, a, wn); % get transfer function of starting filter
    E = Dmp ./ Hd; % compute the error function
    B = polyval(b, exp(j*wn)); % find b(e*j*w)
    W = Wmask.*abs(log(abs(E))./(E -1))./abs(B); % Eqn 11 in Kobayashi & Imai
    %plot(abs(E.*Wmask))
    err = sum(abs(abs(E.*Wmask) - Wmask))
    if err<errmin
        bmin = b;
        amin = a;
    % semilogx(wn,10*log10(abs(Hd)),wn,10*log10(abs(Dmp)))
    end
end



end

