function R = partC_ODE(P)
% PARTC_ODE Part C: Euler and Runge-Kutta (RK4) to solve simple
% harmonic oscillator (SHO), using the given equation and initial
% conditions:
%   y'' + omega^2 y = 0 
%   y(0) = sin(phi), 
%   y'(0) = omega*cos(phi);


    f  = @(t, Y) [Y(2); -P.omega^2 * Y(1)];
    Y0 = [sin(P.phi); P.omega * cos(P.phi)];
    
    hNominal = [0.2 0.1 0.05 0.025 0.0125];
    nh = numel(hNominal);
    
    R.hNom = hNominal;
    R.N = ceil(P.T ./ hNominal);
    R.h = P.T ./ R.N;
    R.euler = zeros(1, nh);
    R.rk4   = zeros(1, nh);
    
    for i = 1:nh
        [t, Yeuler] = euler_solve(f, Y0, P.T, R.N(i));
        [~, Yrk4]   = rk4_solve(f, Y0, P.T, R.N(i));
    
        yExact = sin(P.omega * t + P.phi);
        R.euler(i) = max(abs(Yeuler(1, :) - yExact));
        R.rk4(i)   = max(abs(Yrk4(1, :)   - yExact));
    end
    
    % Observed order = slope of log(error) against log(h), using the actual h
    R.pairEuler = diff(log(R.euler)) ./ diff(log(R.h));
    R.pairRK4   = diff(log(R.rk4))   ./ diff(log(R.h));
    
    R.slopeEuler  = fit_slope(R.h,         R.euler);
    R.slopeRK4    = fit_slope(R.h,         R.rk4);
    R.slopeEuler3 = fit_slope(R.h(end-2:end), R.euler(end-2:end));   
    R.slopeRK43   = fit_slope(R.h(end-2:end), R.rk4(end-2:end));
end 

function slope = fit_slope(h, err)
% FIT_SLOPE  Least-squares slope of log(err) against log(h).
%
%   Inputs:
%     h     - [1 x n] step sizes
%     err   - [1 x n] errors at those steps
%   Outputs:
%     slope - [scalar] least-squares slope of log(err) vs log(h)
    p = polyfit(log(h), log(err), 1);
    slope = p(1);
end