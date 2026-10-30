function print_results(P, A, B, C)
% PRINT_RESULTS  Title plus results for Parts A, B and C.

fprintf('\nPart A: RMS derivative error (k, h, fwd, bwd, cen)\n');
fprintf('%5d %8.4f %9.4g %9.4g %9.4g\n', [A.k(:), A.h(:), A.noisy(:,1:3)]');

fprintf('\nPart B: polynomial regression (m, train, val, resid)\n');
fprintf('%3d %8.4f %8.4f %8.4f\n', [B.m(:), B.trainRMS(:), B.valRMS(:), B.resRMS(:)]');
fprintf('chosen degree m=%d, sigma=%.3g\n', B.mChosen, P.sigma);

fprintf('\nPart C: max error (h, Euler, RK4)\n');
fprintf('%9.5f %12.4e %12.4e\n', [C.h(:), C.euler(:), C.rk4(:)]');
end
