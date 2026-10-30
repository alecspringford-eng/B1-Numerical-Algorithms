function R = partB_poly(P, y)
% PARTB_POLY  Part B: least-squares polynomial fits of degree m = 1..10.
%   Reports conditioning, hold-out validation, and the accuracy of the 
%   analytic derivative of each fit.
%   Validation uses a random 70/30 split; the caller sets rng beforehand.
% 
%   Inputs:
%     P - [struct] from params()
%     y - [N x 1 double] noisy samples on P.t

    trainFraction = 0.7;
    degrees = 1:12;
    nd = numel(degrees);
    N = P.N;
    
    [trainIdx, valIdx] = holdout_split(N, trainFraction);
    inner = 2:N-1;
    mid = round(0.1*N): round(0.9*N);
    
    R.m = degrees;
    R.trainIdx = trainIdx;  R.valIdx = valIdx;
    R.condScaled = zeros(nd, 1);    R.condRaw = zeros(nd, 1);
    R.trainRMS = zeros(nd, 1);  R.valRMS = zeros(nd, 1);    R.resRMS = zeros(nd, 1);
    R.dErrInner = zeros(nd, 1); R.dErrMid = zeros(nd, 1);   R.dErrAll = zeros(nd, 1);
    
    for i = 1:nd
        m = degrees(i);
        [c, V, mu, s] = fit_poly(P.t, y, m);
    
        % Conditioning with/without centring/scaling
        Vraw = P.t .^ (0:m);
        R.condScaled(i) = cond(V' * V);
        R.condRaw(i) = cond(Vraw' * Vraw);
        R.resRMS(i) = rms_error(V * c, y);
    
        % Hold-out: fit on training rows of the same scaled Vandermonde 
        % matrix
        cTrain = V(trainIdx, :) \ y(trainIdx);
        R.trainRMS(i) = rms_error(V(trainIdx, :) * cTrain, y(trainIdx));
        R.valRMS(i) = rms_error(V(valIdx, :) * cTrain, y(valIdx));
    
        % Numerical derivative of the all-data fit against the true 
        % derivative
        [~, dy] = poly_eval(c, mu, s, P.t);
        R.dErrInner(i) = rms_error(dy(inner), P.dytrue(inner));
        R.dErrMid(i) = rms_error(dy(mid), P.dytrue(mid));
        R.dErrAll(i) = rms_error(dy, P.dytrue);
    end

    [~, ~, dC] = diff_schemes(y, P.dt);
    R.dC = dC;                                          
    R.centralBase = rms_error(dC, P.dytrue(inner));
    R.centralMid  = rms_error(dC(mid(1)-1 : mid(end)-1), P.dytrue(mid));

    % Central difference at best sub-sampled step, to be used for
    % comparison plot
    kBest = 128;
    idx = 1:kBest:N;
    [~, ~, R.dC128] = diff_schemes(y(idx), kBest * P.dt);
    R.t128 = P.t(idx(2:end-1));
    R.mChosen = 10;
    [c, ~, mu, s] = fit_poly(P.t, y, R.mChosen);
    [R.yfit, R.dyfit] = poly_eval(c, mu, s, P.t);

end


function [trainIdx, valIdx] = holdout_split(N, trainFraction)
% HOLDOUT_SPLIT Random split into sorted training and validation indices
    perm = randperm(N);
    nTrain = round(trainFraction * N);
    trainIdx = sort(perm(1:nTrain));
    valIdx = sort(perm(nTrain + 1:end));
end 