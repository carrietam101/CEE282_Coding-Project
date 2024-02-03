function RC_Plot_Errors(load_norm_errors, energy_norm_errors, APRATIOS)

% Plot the load norm and energy norm error indices against applied load ratio

%   load_norm_errors    : row or column vector containing the load norm error index at each increment
%   energy_norm_errors  : row or column vector containing the energy norm error index at each increment
%   APRATIOS            : row or column vector containing the applied load ratio at each increment

%% Check dimensions of arguments
% Check if all arguments are vectors
if ~isvector(load_norm_errors)
    error('load_norm_errors is not a vector');
end

if ~isvector(energy_norm_errors)
    error('energy_norm_errors is not a vector');
end

if ~isvector(APRATIOS)
    error('APRATIOS is not a vector');
end

% Ensure that all vectors have compatible lengths
if length(load_norm_errors) ~= length(energy_norm_errors)
    error('load_norm_errors and energy_norm_errors have unequal lengths');
end

if length(APRATIOS) < length(load_norm_errors)
    error('APRATIOS must have length greater than or equal to the error vectors');
end

% If the analysis is continued from a previous analysis, length of APRATIOS will be larger than the error
% vectors since information about the errors is not retained from previous analyses. In this case, only the
% trailing end of APRATIOS is used.
APRATIOS = APRATIOS((end - length(load_norm_errors) + 1):end);

%% Create the plot
figure;
hold on;
plot(APRATIOS, load_norm_errors, 'b-', 'LineWidth', 1.5);
plot(APRATIOS, energy_norm_errors, 'r--', 'LineWidth', 1.5);
hold off;
set(gca, 'YScale', 'log', 'FontSize', 12);
xlabel('Applied load ratio');
ylabel('Error index');
grid on;

% If Mastan2 is running, color of title and legend text is changed to white, so they need to be explicitly
% set to black
title('Error indices vs. Applied load ratio', 'Color', 'k');
hleg = legend('Load norm error index', 'Energy norm error index', 'Location', 'Best');
set(hleg, 'TextColor', 'k');

end
