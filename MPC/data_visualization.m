clc;
clear;

%% Elevation

dataSim = load("ElevationSim.mat");
%dataReal = load("ElevationReal.mat");

time = dataSim.Elevation.time;
valuesRef = dataSim.Elevation.signals.values(1:end,1);
valuesSim = dataSim.Elevation.signals.values(1:end,2);
%valuesReal = dataReal.Elevation.signals.values(1:end,2);

grid on;
hold on;

title("Elevation angle", 'FontSize', 16)
xlabel("t[s]", 'FontSize', 16)
ylabel("\epsilon(t) [deg]", 'FontSize', 16)
xlim([0 45]);
ylim([-10 10]);

plot(time, valuesRef, 'b', 'LineWidth', 2);
plot(time, valuesSim, 'r', 'LineWidth', 2);
%plot(time, valuesReal, 'g', 'LineWidth', 2);

legend("Reference", "MPC sim.", "Real hardware", "se", 'FontSize', 16)


hold off;

%% Pitch

dataSim = load("PitchSim.mat");
%dataReal = load("PitchReal.mat");

time = dataSim.Pitch.time;
valuesRef = dataSim.Pitch.signals.values(1:end,1);
valuesSim = dataSim.Pitch.signals.values(1:end,2);
%valuesReal = dataReal.Pitch.signals.values(1:end,2);

grid on;
hold on;

title("Pitch angle", 'FontSize', 16)
xlabel("t[s]", 'FontSize', 16)
ylabel("\theta(t) [deg]", 'FontSize', 16)
xlim([0 45]);
ylim([-30 30]);

plot(time, valuesRef, 'b', 'LineWidth', 2);
plot(time, valuesSim, 'r', 'LineWidth', 2);
%plot(time, valuesReal, 'g', 'LineWidth', 2);

legend("Reference", "MPC sim.", "Real hardware", "se", 'FontSize', 16)


hold off;



%% Travel

dataSim = load("TravelSim.mat");
%dataReal = load("TravelReal.mat");

time = dataSim.Travel.time;
valuesRef = dataSim.Travel.signals.values(1:end,1);
valuesSim = dataSim.Travel.signals.values(1:end,2);
%valuesReal = dataReal.Travel.signals.values(1:end,2);

grid on;
hold on;

title("Travel angle", 'FontSize', 16)
xlabel("t[s]", 'FontSize', 16)
ylabel("\lambda(t) [deg]", 'FontSize', 16)
xlim([0 45]);
ylim([-40 40]);

plot(time, valuesRef, 'b', 'LineWidth', 2);
plot(time, valuesSim, 'r', 'LineWidth', 2);
%plot(time, valuesReal, 'g', 'LineWidth', 2);

legend("Reference", "MPC sim.", "Real hardware", "se", 'FontSize', 16)


hold off;

%% Input

dataSim = load("InputSim.mat");
%dataReal = load("InputReal.mat");

time = dataSim.Input.time;
valuesSimU1 = dataSim.Input.signals.values(1:end,1);
valuesSimU2 = dataSim.Input.signals.values(1:end,2);
%valuesRealU1 = dataReal.Input.signals.values(1:end,1);
%valuesRealU2 = dataReal.Input.signals.values(1:end,2);

grid on;
hold on;

title("Control Input")
xlabel("t[s]")
ylabel("u [V]")
xlim([0 45]);
ylim([-7 7]);

plot(time, valuesSimU1, 'b', 'LineWidth', 2);
plot(time, valuesSimU2, 'r', 'LineWidth', 2);
%plot(time, valuesRealU1, 'g', 'LineWidth', 2);
%plot(time, valuesRealU2, 'black', 'LineWidth', 2);

legend("Simulation u_1", "Simulation u_2", "Real hardware u_1", "Real hardware u_2", "se")


hold off;