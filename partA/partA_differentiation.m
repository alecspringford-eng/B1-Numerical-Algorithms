function R = partA_differentiation(P, y)
% PARTA_DIFFERENTIATION Part A: derivative errors against different step
% sizes h.
% Also contains function for noise free signal

    ks = [1 2 4 8 16 32 64 128 256 512];
    nk = numel(ks);
    
    R.k = ks;
    R.h = ks * P.dt;
    R.noisy = zeros(nk, 3);
    R.clean = zeros(nk, 3);
    
    for i = 1:nk
        R.noisy(i, :) = errors(y, ks(i), P);
        R.clean(i, :) = errors(P.ytrue, ks(i), P);
    end
    
    % Noise-only RMSE. Independent errors of std.
    hCol = R.h(:);
    R.noisePred = [sqrt(2)*P.sigma./hCol, sqrt(2)*P.sigma./hCol, P.sigma./(sqrt(2)*hCol)];

    [R.cMin, iBest] = min(R.noisy(:,3)); R.hOptC = R.h(iBest);
    [R.fMin, iBest] = min(R.noisy(:,1)); R.hOptF = R.h(iBest);


    % Seed-to-seed variability of the central-difference error
    R.seedK = [64 256 512];
    R.seedMean = zeros(1,3);
    R.seedStd = zeros(1,3);
    nSeeds = 50;
    rngState = rng;                 % restore after loop so Part B split follows rng(1)
    for j = 1:numel(R.seedK)
        seedErrors = zeros(nSeeds, 1);
        for seed = 1:nSeeds
            rng(seed);
            noisySignal = P.ytrue + P.sigma * randn(P.N, 1);
            e = errors(noisySignal, R.seedK(j), P);
            seedErrors(seed) = e(3);
        end
        R.seedMean(j) = mean(seedErrors);
        R.seedStd(j) = std(seedErrors);
    end

    rng(rngState);

    %Truncation order from noise-free error
    fitRange = ks <= 16;
    R.cleanSlope = zeros(1,3);
    for q = 1:3
        p = polyfit(log(R.h(fitRange)), log(R.clean(fitRange, q)'), 1);
        R.cleanSlope(q) = p(1);
    end
end



function errs = errors(signal, k, P)
% ERRORS samples reference signal and compares with differentiated noisy
% signal to look at errors, for the interior points.

    idx = 1:k:P.N;
    h = k * P.dt;
    
    ref = P.dytrue(idx);
    ref = ref(2:end - 1);       %Only interior points
    
    [dF, dB, dC] = diff_schemes(signal(idx), h);
    errs = [rms_error(dF, ref), rms_error(dB, ref), rms_error(dC, ref)];
end