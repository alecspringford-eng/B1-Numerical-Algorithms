function [c, V, mu, s] = fit_poly(t, y, m)
% FIT_POLY Least-squares polynomial fit of degree, centered and scaled for
% time. 
%   Backslash avoids squaring the condition number as V'*V gives cond(V)^2
%
%   Inputs:
%     t  - [N x 1 double] sample times
%     y  - [N x 1 double] sample values
%     m  - [scalar int] polynomial degree
%   Outputs:
%     c  - [(m+1) x 1] coefficients in the scaled variable, c(j+1) multiplies tc^j
%     V  - [N x (m+1)] Vandermonde matrix of the scaled times tc
%     mu - [scalar] centre of t used for scaling
%     s  - [scalar] half-width of t used for scaling

    t = t(:);
    y = y(:);
    
    mu = (t(1) + t(end))/2;
    s = (t(end) - t(1))/2;
    tc = (t - mu)/s;
    
    V = tc .^ (0:m);
    c = V \ y;

end 