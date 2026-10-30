function e = rms_error(estimated, reference)
% RMS_ERROR Finds the L2 error, normalized by the number of points, i.e.
% the root mean squared error (RMSE)
%
%   Inputs:
%     est - [array, any shape] estimated values
%     ref - [array, same numel as est] reference (true) values
%   Outputs:
%     e   - [scalar double] root-mean-square of (est - ref)
    e = sqrt(mean((estimated(:) - reference(:)).^2));
end