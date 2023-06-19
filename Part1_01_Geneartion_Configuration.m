%% Data Generation Configuration 

disp('<Running Part1_01_Generation_Configuration.mlx>')
%% General

% gravity
g = 9.81; %m/s^2
mg = 0.001*g;

%% file save mode
flg_file_save = 1;

%% simulation
% number of time window
NumSeg = 200; % one event in one window (should >= 20)

% simulation time step
SimStep = 0.01; % for variable-step simulation better no bigger than 0.005

% downsampling setup
hf = 1;
lf = 1;
ds_ratio = lf/hf;

% time window
window = 60; %sec

% number of frames in a time window
seg = window/SimStep; 

% timeseries
TimeSeries = [1:seg]*SimStep;

% figure
flg_img_show = 1;
%% 1 Structure

%% structure
% number of DOF for the structure
nDOF = 6;
%% 2 Event

%% event and fault generation mode

%% event type
NumEventType = 4; % no event / earthquake / impact / strong wind
NumFaultType = 4; % fault free / bias / drift / spike



% 2.1 Ambient Vibration

%% events - ambient vibration - assume to be 0 for simplification purpose
% magnitude of ambient vibration
power = 0.1;% dBW %power of white Gaussian Noise
target_noise_level = 0.5*mg; % rescale factor after normalisation
% 2.2 Earthquake

%% events - earthquake
% earthquake magnitude rescale factor
factoreq = (0.5 + 1.5*rand(1)) * g; %m/s^2
omegag_base = 15;
zetag_base = 0.6;  
S0_base = 0.0049;
fac_time_base = 12.21;
pos_time1_base = 0.1;
%pos_time2_base = 0.5;
% 2.3 Impact

%% events - impact
% define impact intensity baseline
ImpactBaseline = 2*g; % unit:m/s^2
% 2.4 Strong Wind

%% events - strong wind
% standard deviation of strong wind events duration
MaxSig = 3+4*rand(1);

% scale factor of strong wind events 
MaxScale = 1+3*rand(1);
%% 3 Sensor Faults Generation
% 3.1 Bias

% bias range lower boundary and upper boundary defining
bias_lb = 1e-4;
bias_ub = 10;
% 3.2 Drift

% max ratio of drift (absolute value)
k_cap = 5e-2;
% 3.3 Spike

% spike magnitude range lower boundary and upper boundary defining
spike_lb = 1e-2;
spike_ub = 5e1;
%% 4 Measurement Noise

% noise - described by rooted mean square

% snr = 20; 
% rms_measure = 0.5*mg;

ns_level = -40; % noise level (dB) 
%% 
%%