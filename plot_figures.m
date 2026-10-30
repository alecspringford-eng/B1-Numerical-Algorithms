function plot_figures(P, A, B, C)
% PLOT_FIGURES  Five separate figures, one per result.
%   Uses only stored results; no numerical work is done here.
    plot_derivative_error(A);
    plot_conditioning(B);
    plot_holdout(P, B);
    plot_derivative_estimates(P, B);
    plot_ode_convergence(C);
end

function fig = new_figure(name)
    fig = figure('Name', name, 'Units', 'centimeters', ...
                 'Position', [2 2 13 9], 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', 9);
end

function plot_derivative_error(A)
    fig = new_figure('fig_A_derivative_error');
    col = lines(7);

    loglog(A.h, A.noisy(:,1), 'o-', 'Color', col(1,:), 'MarkerSize', 4); hold on
    loglog(A.h, A.noisy(:,2), 's-', 'Color', col(2,:), 'MarkerSize', 4);
    loglog(A.h, A.noisy(:,3), '^-', 'Color', col(3,:), 'MarkerSize', 4);
    loglog(A.h, A.clean(:,1), '--', 'Color', col(1,:));
    loglog(A.h, A.clean(:,3), '--', 'Color', col(3,:));
    loglog(A.h, A.noisePred(:,3), ':k');
    xline(A.hOptC, ':', 'Color', [.4 .4 .4]);

    grid on
    xlabel('step size h = k\Deltat');
    ylabel('RMS error of y''(t)');
    title('Derivative error vs h');
    legend({'forward', 'backward', 'central', 'forward, noise-free', 'central, noise-free', ...
            'predicted noise error (central)', 'best step size h (central)'}, ...
           'Location', 'southeast', 'FontSize', 7);
end

function plot_conditioning(B)
    fig = new_figure('fig_B_conditioning');
    col = lines(7);

    semilogy(B.m, B.condRaw,    's-', 'Color', col(4,:), 'MarkerSize', 4); hold on
    semilogy(B.m, B.condScaled, 'o-', 'Color', col(5,:), 'MarkerSize', 4);
    yline(1/eps, ':k');

    grid on
    xlabel('polynomial degree m');
    ylabel('cond(V''V)');
    title('Conditioning of the normal equations');
    legend({'raw t', 'centred & scaled t', '1/\epsilon'}, 'Location', 'northwest', 'FontSize', 7);
end

function plot_holdout(P, B)
    fig = new_figure('fig_B_holdout');
    col = lines(7);

    semilogy(B.m, B.trainRMS, 'o-',  'Color', col(1,:), 'MarkerSize', 4); hold on
    semilogy(B.m, B.valRMS,   's--', 'Color', col(2,:), 'MarkerSize', 4);
    yline(P.sigma, ':k');
    ylim([5e-3 1]);

    grid on
    xlabel('polynomial degree m');
    ylabel('RMS residual');
    title('Training and hold-out error (random 70/30 split)');
    legend({'train', 'hold-out', '\sigma'}, 'Location', 'northeast', 'FontSize', 7);
end

function plot_derivative_estimates(P, B)
    fig = new_figure('fig_B_derivative_estimates');
    col = lines(7);

    plot(P.t(2:7:end-1), B.dC(1:7:end), '.', 'Color', [.7 .7 .7], 'MarkerSize', 3); hold on
    plot(B.t128, B.dC128, 'b--', 'LineWidth', 0.8);
    plot(P.t, P.dytrue, 'k-', 'LineWidth', 1);
    plot(P.t, B.dyfit, '-', 'Color', col(3,:), 'LineWidth', 1);

    xlim([0 P.T]);  ylim([-4 6]);  grid on
    xlabel('t');
    ylabel('y''(t)');
    title('Derivative estimates');
    legend({'central, h=\Deltat', 'central, h=128\Deltat', 'true y''', ...
            sprintf('poly m=%d', B.mChosen)}, ...
           'Location', 'north', 'NumColumns', 2, 'FontSize', 7);
end

function plot_ode_convergence(C)
    fig = new_figure('fig_C_ode_convergence');
    col = lines(7);

    loglog(C.h, C.euler, 'o-', 'Color', col(1,:), 'MarkerSize', 4); hold on
    loglog(C.h, C.rk4,   's-', 'Color', col(2,:), 'MarkerSize', 4);
    loglog(C.h, C.euler(end) * (C.h / C.h(end)).^1, '--', 'Color', col(1,:));
    loglog(C.h, C.rk4(end)   * (C.h / C.h(end)).^4, '--', 'Color', col(2,:));

    grid on
    xlabel('step size h');
    ylabel('max |y - y_{exact}|');
    title('Global error, T = 2\pi');
    legend({sprintf('Euler (measured order %.2f, smallest 3 h)', C.slopeEuler3), ...
            sprintf('RK4 (measured order %.2f, smallest 3 h)',   C.slopeRK43), ...
            'Euler expected order 1', 'RK4 expected order 4'}, ...
           'Location', 'northwest', 'FontSize', 7);
end
