function [dF, dB, dC] = diff_schemes(y,h)
%   DIFF_SCHEMES Finds the Forward (dF), backward (dB) and central difference
%   (dC) estimates of the derivative of y (y'(t)).
%   y: samples on a uniform grid with spacing h
%   Each output is evaluated at the interior points so that all three
%   schemes can use the same evaluation points.
%
%     forward : (y(i+1) - y(i)) / h
%     backward: (y(i)   - y(i-1)) / h
%     central : (y(i+1) - y(i-1)) / (2h)
%
%   Inputs:
%     y  - [N x 1 double] (or row) samples y(t_i) on a uniform grid
%     h  - [scalar double] grid spacing between consecutive samples
%   Outputs:
%     dF - [(N-2) x 1 double] forward-difference estimate
%     dB - [(N-2) x 1 double] backward-difference estimate
%     dC - [(N-2) x 1 double] central-difference estimate

    y = y(:);
    dF = (y(3:end) - y(2:end-1)) / h;
    dB = (y(2:end-1) - y(1:end-2)) / h;
    dC = (y(3:end) - y(1:end-2)) / (2*h);
end