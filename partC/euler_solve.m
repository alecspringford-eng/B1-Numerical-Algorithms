function [t, Y] = euler_solve(f, Y0, T, N)
% EULER_SOLVE  Forward Euler for Y' = f(t, Y) using N steps of size h = T/N.
%
%   Inputs:
%     f  - [function handle] f(t, Y) -> [d x 1] derivative; t scalar, Y [d x 1]
%     Y0 - [d x 1] initial state
%     T  - [scalar] final time
%     N  - [scalar int] number of steps
%   Outputs:
%     t  - [1 x (N+1)] time points
%     Y  - [d x (N+1)] state at each time point

    h = T / N;
    t = (0:N) * h;
    
    Y = zeros(numel(Y0), N + 1);
    Y(:, 1) = Y0(:);
    
    for n = 1:N
        Y(:, n+1) = Y(:, n) + h * f(t(n), Y(:, n));
    end
end
