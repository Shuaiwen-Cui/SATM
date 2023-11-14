%% Structure Modelling
% Structure to Simulate - n DOF Shear Building
% 
% Simulation Method - State Space Model
% 
% Theoretical Foundation - Equations of Motion: M*Xdotdot + C*Xdot + K*X = G*P(t)  
% ==>  State-Space Model
% 
% 
% 
% Xdot = Af*X + Bf*F- State Transfer
% 
% Y = Cf*X + Df*F - Measurement
% 
% 
% 
% X - state matrix - 2nx1 - [position; velocity]
% 
% Xdot - state matrix derivative - 2nx1 - [velocity; acceleration]
% 
% Y - measurement matrix - nx1 - [ acceleration]
% 
% 

disp('<Running Part1_02_Structure_Modelling.mlx>')
%% 1 Parameter Configuration

me = 2.9e4; % mass - unit kg
ke = 1.2e9; % stiffness - unit N/m
zeta = 0.02; % damping ratio - actual damp/ critical damp
% 1.1 M - Mass Matrix 

M = me*eye(nDOF);
% 1.2 K - Stiffness Matrix

K = 2*ke*eye(nDOF); K(nDOF,nDOF)=ke;
for i=1:nDOF
    if i<nDOF
        K(i,i+1) = -ke;
    end
    if i>1
        K(i,i-1) = -ke;
    end
end
% 1.3 C - Damping Matrix

% eigen vectors phi and eigen values - # of phi and lambda should all be 10
[phi, lambda] = eig(inv(M)*K);

% identify the characteristic freq
omg = sqrt(diag(lambda));% in rad/s
w = omg/2/pi; % convert in to Hz

% damping matrix
C = inv(phi')*(diag(omg)*2*zeta*me)*inv(phi);
%% 2 Sub Model for Different Events
% consistent in form
% 2.1 General Simulation (used for ambient vibration)

% Aa - state to new state - 2nx2n
% Ba - input to new state - 2nxn
% Ca - state to measurement - [1-3]nxn
% Da - input to% measurement - [1-3]nxn
% Ga - load effect matrix % nxn
% F - input vector - nx1

Ga = eye(nDOF); %varying with different load, by default eyes(nDoF)

Aa = [  zeros(nDOF)     eye(nDOF)
             -M\K           -M\C];

Ba = [ zeros(nDOF);
              M\Ga;];

Ca = [       -M\K           -M\C];  

Da = [        M\Ga;];
% 2.2 Earthquake Simulation

% Ae - state to new state - 2nx2n
% Be - input to new state - 2nxn
% Ce - state to measurement - [1-3]nxn
% De - input to% measurement - [1-3]nxn
% Ge - load effect matrix % nxn
% F - input vector - nx1

Ge = -1*eye(nDOF); %varying with different load, by default eyes(nDoF)

Ae = [  zeros(nDOF)     eye(nDOF)
             -M\K           -M\C];

Be = [ zeros(nDOF);
              M\Ge;];

Ce = [       -M\K           -M\C];  

De = [        M\Ge;];
% 2.3 Impact Simulation

% Ai - state to new state - 2nx2n
% Bi - input to new state - 2nxn
% Ci - state to measurement - [1-3]nxn
% Di - input to% measurement - [1-3]nxn
% Gi - load effect matrix % nxn
% F - input vector - nx1

Gi = eye(nDOF); %varying with different load, by default eyes(nDoF)

Ai = [  zeros(nDOF)     eye(nDOF)
             -M\K           -M\C];

Bi = [ zeros(nDOF);
              M\Gi;];

Ci = [       -M\K           -M\C];  

Di = [        M\Gi;];
% 2.4  Wind Simulation

% Aw - state to new state - 2nx2n
% Bw - input to new state - 2nxn
% Cw - state to measurement - [1-3]nxn
% Dw - input to% measurement - [1-3]nxn
% Gw - load effect matrix % nxn
% F - input vector - nx1

Gw = eye(nDOF); %varying with different load, by default eyes(nDoF)

Aw = [  zeros(nDOF)     eye(nDOF)
             -M\K           -M\C];

Bw = [ zeros(nDOF);
              M\Gw;];

Cw = [       -M\K           -M\C];  

Dw = [        M\Gw;];