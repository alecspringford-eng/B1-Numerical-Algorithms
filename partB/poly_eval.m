function [yfit, dyfit] = poly_eval(c, mu, s, t)
% POLY_EVAL Evaluates the fitted polynomial and its numerical derivative at
% t.
%   Inputs:
%     c  - [(m+1) x 1] coefficients from fit_poly
%     mu - [scalar] centre from fit_poly
%     s  - [scalar] half-width from fit_poly
%     t  - [n x 1 double] times at which to evaluate
%   Outputs:
%     yfit  - [n x 1] fitted polynomial values
%     dyfit - [n x 1] analytic derivative d(yfit)/dt

    t = t(:);
    m = numel(c) - 1;
    tc = (t - mu)/s;
    
    yfit = (tc .^ (0:m)) * c;
    
    if m == 0
        dyfit = zeros(size(t));
    else
        dVdtc = [zeros(size(t)), (1:m) .* tc .^ (0:m-1)];
        dyfit = (dVdtc * c) / s;
    end
end