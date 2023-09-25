clear all;
close all;
clc;

a = 18;
b = 21;



rE = 6378; % Radius of Earth (km)

%% First orbit
% Initial Conditions
mu = 3.986e5;
r0 = [rE+700,0,0];
v0 = [0,0,sqrt(mu/r0(1))];

s0 = [r0;v0]; % State Vector
n0 = sqrt(mu/r0(1)^3);
tMax0 = (2*pi)/n0;

% Timespans
timeStep0 = 1;
timeSpan0 = 0:1:tMax0; % (seconds) to end of period

%% Second Orbit


r1 = [rE+900,0,0];
v1 = [0,0,sqrt(mu/r1(1))];



s1 = [r1;v1]; % State Vector
n1 = sqrt(mu/r1(1)^3);
tMax1 = (2*pi)/n1;

% Timespans
timeStep1 = 1;
timeSpan1 = 0:1:tMax1; % (seconds) to end of period

%% transfer orbit
h1 = sqrt(2*mu)*sqrt((r0(1)*r1(1))/(r0(1)+r1(1)));
va1 = h1/r0(1);
vb1 = h1/r1(1);

delVa = abs(va1-v0(1,3))
delVb = abs(v1(1,3)-vb1)

delVTot = delVa + delVb

tReq = pi*sqrt(((0.5*(r0(1)+r1(1)))^3)/mu)
timeSpanT = 0:1:tReq;

sT = [r0,[0,0,va1]];


%%


% Output timeSPan and Solutions (sol)
% diffEq used to pass in variables
[~, sol0] = ode45(@(t0,s0)diffEq(t0,s0,mu), timeSpan0, s0); % Need 3 args
[~, sol1] = ode45(@(t1,s1)diffEq(t1,s1,mu), timeSpan1, s1); % Need 3 args
[~, solT] = ode45(@(tT,sT)diffEq(tT,sT,mu), timeSpanT, sT); % Need 3 arg


% Plot Solutions
figure;

hold on; grid on; grid minor;
% Create Earth
[X, Y, Z] = sphere(); 

surf(X*rE, Y*rE, Z*rE)
plot3(sol0(:,1),sol0(:,2),sol0(:,3),'-b','LineWidth',1) % Plot in 3d
plot3(sol1(:,1),sol1(:,2),sol1(:,3),'-r','LineWidth',1) % Plot in 3d
plot3(solT(:,1),solT(:,2),solT(:,3),'-b','LineWidth',1) % Plot in 3d




function sdot = diffEq(t0,s0,mu)

% r was first three elements, v was second three
r_0 = s0(1:3);
v_0 = s0(4:6);

sdot(1:3,1) = v_0; %First three elements are velocity
sdot(4:6,1) = (-mu*r_0)/(norm(r_0))^3; % Actual two body EOM


end
