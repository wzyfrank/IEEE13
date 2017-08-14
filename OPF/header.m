clear all;
clc;
close all;

format short g;
% import the loads
loads = load('load_vec.mat');
loads = loads.load;
loads = transpose(loads);
loads_base = loads;

% import the load types
load_type = load('load_type.mat');
load_type = load_type.load_type;
load_type = transpose(load_type);

% import Zbus and Cbus
Zbus = load('Zbus.mat');
Zbus = Zbus.Zbus;
Cbus = load('Cbus.mat');
Cbus = Cbus.Cbus;

% init the voltage profile
[num,txt,raw] = xlsread('loads data.csv');
N_loads = size(num, 1);
Voltage_profile = ones(N_loads, 6);
Voltage_profile(:, 2) = 0;
Voltage_profile(:, 4) = -120;
Voltage_profile(:, 6) = 120;

% base phase voltage (V)
Vbase = 4160 / sqrt(3);


% start iteration
for iter = 1:3
    % output loads
    loads = load_process(Voltage_profile, Vbase);
    [Voltage_profile, phasors] = ieee13_iter(loads, Zbus, Cbus);
    % [Voltage_profile, phasors] = ieee13_yalmip(loads, Zbus, Cbus);
    %[Voltage_profile, phasors] = ieee13_hybrid2(loads, Zbus, Cbus);
    % update the loads
end
