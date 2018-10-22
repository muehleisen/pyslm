function [bu, au] = unwarp(bw, aw, lambda);

    [zw,pw,kw] = tf2zpk(bw, aw);

    zu = (zw + lambda)./(1 + zw* lambda);
    pu = (pw + lambda)./(1 + pw*lambda);
    
    [bu, au] = zp2tf(zu, pu, kw);


end