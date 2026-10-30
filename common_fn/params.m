function P = params()
% PARAMS Fixed parameters shared by every scripts, come from the given
% values in the sheet.

    P.T = 2*pi;      % total time interval
    P.dt = 0.001;     % base sampling step
    P.omega = 2.0;       % angular frequency
    P.phi = pi/6;      % phase
    P.sigma = 0.01;      % noise standard deviation
    
    P.N = floor(P.T/P.dt) + 1;
    P.t = (0:P.N-1)' * P.dt;
    
    % Noise free signal and its derivative to compare against
    P.ytrue = sin(P.omega*P.t + P.phi);
    P.dytrue = P.omega * cos(P.omega * P.t + P.phi);

end