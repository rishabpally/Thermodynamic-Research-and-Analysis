close all
clear
clc
format longG
Constants
m_prop = 1; %.06852; % [slugs] = 1kg
g0 = 9.81; %32.2; % [ft/s^2]
Load data
data1 = load ('test1');
Calculate time and frequency
freq = 1.652 * 1000; % [Hz]
% time = (1 / freq) .* linspace(0, length(data1), length(data1))';
% ------------------------------------------ %
% ------------------GOOD TEST--------------- %
% ------------------------------------------ %
fit
% means = movmean(data1(:, 3), 10);
% points = find(diff(means) > .05);
% time_new = time(points);
%
% % Polyfit to clean up line
% c1 = polyfit(time(points), data1(points, 3), 8);
% f1 = c1(1) .* time(points).^8 + c1(2) .* time(points).^7 + c1(3) .*
time(points).^6 + c1(4) .* time(points).^5 + c1(5) .* time(points).^4
+ c1(6) .* time(points).^3 + c1(7) .* time(points) .^2 + c1(8) .*
time(points) + c1(9);
%
% % calc slope

% slope = (data1(points(1), 3) - data1(points(end), 3)) /
(time(points(1)) - time(points(end)));
% y = slope .* time(points) - slope .* time(points(end)) +
data1(points(end), 3);
%
% fx = f1 - y;
%
% fitobject = fit(time_new, fx, 'cubicinterp');
%
% I = trapz(time(points), fx);
% I2 = integrate(fitobject, time_new(end), time_new(1));
%
% I_sp = I / (m_prop * g0);
% I_sp2 = I2 / (m_prop * g0);
Plotting
plot(time(points), data1(points, 3)); hold on plot(time(points), zeros(1, length(points))); plot(time(points),
f1) plot(time(points), y) xlabel('Time [s]') ylabel('Thrust [lbf]')
% ------------------------------------------ %
% ------------------ALL FILES--------------- %
% ------------------------------------------ %
Combine for all files
trials = 16;
N = 0;
for v = 1:trials
% Set to make sure file 12 is not accounted for in calculations
if v == 12
I_sp_mean(12) = NaN;
SEM_I_sp(12) = NaN;
continue
end
%Read in data file v
filename = sprintf('%s%d','LA_Test_FixedMass_Trial',v);
inData = readmatrix(filename);
% Take out first 4 rows
inData = inData(4:end, :);
inData(:,3) = inData(:,3) .* 4.4482216; % Converting to N
time = (1 / freq) .* linspace(0, length(inData), length(inData))';
means = movmean(inData(:, 3), 20);
points = find(diff(means) > .5);
time_new = time(points);
% Fitting line to thrust data

c1 = polyfit(time(points), inData(points, 3), 5);
f1 = c1(1) .* time(points).^5 + c1(2) .* time(points).^4 +
c1(3) .* time(points).^3 + c1(4) ...
.* time(points).^2 + c1(5) .* time(points) + c1(6);
% calculate slope
slope = (inData(points(1), 3) - inData(points(end), 3)) /
(time(points(1)) - time(points(end)));
y = slope .* time(points) - slope .* time(points(end)) +
inData(points(end), 3);
% Zero out line
fx = f1 - y;
% Second integration test
fitobject = fit(time_new, fx, 'cubicinterp');
% Calculate Impulse, I
I(v) = trapz(time(points), fx);
I2(v) = integrate(fitobject, time_new(end), time_new(1));
% Calculate specific impulse
I_sp(v) = I(v) / (m_prop * g0);
I_sp2(v) = I2(v) / (m_prop * g0);
% I_sp Statistical data
I_sp_mean = mean(I_sp);
I_sp_mean2 = mean(I_sp2);
peak_thrust(v) = max(inData(:,3));
thrust_time(v) = time_new(end);
peak_thrust_mean = mean(peak_thrust);
peak_thrust_std = std(peak_thrust);
thrust_time_mean = mean(thrust_time);
thrust_time_std = std(thrust_time);
% Number of entries considered
N = N + length(time_new);
% Standard error of the mean
SEM_I_sp(v) = std(I_sp(1:v)) / N;
% Confidence Interval
CI1 = [I_sp_mean - 1.96 * SEM_I_sp(end), I_sp_mean + 1.96 *
SEM_I_sp(end)]; % 95% CI
CI2 = [I_sp_mean - 2.24 * SEM_I_sp(end), I_sp_mean + 2.24 *
SEM_I_sp(end)]; % 97.5% CI
CI3 = [I_sp_mean - 2.58 * SEM_I_sp(end), I_sp_mean + 2.58 *
SEM_I_sp(end)]; % 99% CI
% Remove 12th Index from plot
for i = 12:length(SEM_I_sp) - 1

SEM_I_sp(i) = SEM_I_sp(i + 1);
end
% Manipulate this section depending on desired plots
% figure()
% %plot(time, inData(:, 3),'k', 'LineWidth', 1.5)
% hold on
% plot(time_new, inData(points, 3), 'LineWidth', 1.5);
% hold on
% % plot(time_new, y, 'LineWidth', 1.5);
% plot(time_new, fx, 'LineWidth', 1.5);
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% legend('Uncondtitioned Thrust Curve', ' Conditioned Thrust
Curve')
% % legend('Thrust Force (Unconditioned)', 'Thrust Force
(Conditioned)')
% title('Representative Force Curve')
% grid on
% hold off
end
dataset = 1:1:15;
% Sem plot
plot(dataset, SEM_I_sp(1:15), '-o', 'MarkerSize', 8, 'LineWidth', 1.5)
grid on
xlabel('Number of Samples (N)')
ylabel('Standard Error of the Mean (SEM) [s]')
title('Standard Error vs Datasets')
legend('SEM')