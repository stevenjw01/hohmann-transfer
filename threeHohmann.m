clear all;
close all;
clc;

format long g

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


r1 = [r0(1)+a,0,0];
v1 = [0,0,sqrt(mu/r1(1))];



s1 = [r1;v1]; % State Vector
n1 = sqrt(mu/r1(1)^3);
tMax1 = (2*pi)/n1;

% Timespans
timeStep1 = 1;
timeSpan1 = 0:1:tMax1; % (seconds) to end of period

%% 1st transfer orbit
h1 = sqrt(2*mu)*sqrt((r0(1)*r1(1))/(r0(1)+r1(1)));
va1 = h1/r0(1);
vb1 = h1/r1(1);

delV1 = abs(va1-v0(1,3))+abs(v1(1,3)-vb1)

tReq = pi*sqrt(((0.5*(r0(1)+r1(1)))^3)/mu)
timeSpanT = 0:1:tReq;

sT = [r0,[0,0,va1]];
%% Third Orbit


r2 = [r0(1)+a+b,0,0];
v2 = [0,0,sqrt(mu/r2(1))];



s2 = [r2;v2]; % State Vector
n2 = sqrt(mu/r2(1)^3);
tMax2 = (2*pi)/n2;

% Timespans
timeStep2 = 1;
timeSpan2 = 0:1:tMax2; % (seconds) to end of period
%% 2nd transfer orbit
h2 = sqrt(2*mu)*sqrt((r1(1)*r2(1))/(r1(1)+r2(1)));
va2 = h2/r1(1);
vb2 = h2/r2(1);

delV2 = abs(va2-v1(1,3))+abs(v2(1,3)-vb2)

tReq2 = pi*sqrt(((0.5*(r1(1)+r2(1)))^3)/mu)
timeSpanT2 = 0:1:tReq2;

sT2 = [r1,[0,0,va2]];

%% Fourth Orbit


r3 = [rE+900,0,0];
v3 = [0,0,sqrt(mu/r3(1))];



s3 = [r3;v3]; % State Vector
n3 = sqrt(mu/r3(1)^3);
tMax3 = (2*pi)/n3;

% Timespans
timeStep3 = 1;
timeSpan3 = 0:1:tMax3; % (seconds) to end of period
%% 3rd transfer orbit
h3 = sqrt(2*mu)*sqrt((r3(1)*r2(1))/(r3(1)+r2(1)));
va3 = h3/r2(1);
vb3 = h3/r3(1);



delV3 = abs(va3-v2(1,3))+abs(v3(1,3)-vb3)

tReq3 = pi*sqrt(((0.5*(r2(1)+r3(1)))^3)/mu)
timeSpanT3 = 0:1:tReq3;

sT3 = [r2,[0,0,va3]];
%%


% Output timeSPan and Solutions (sol)
% diffEq used to pass in variables
[~, sol0] = ode45(@(t0,s0)diffEq(t0,s0,mu), timeSpan0, s0); % Need 3 args
[~, sol1] = ode45(@(t1,s1)diffEq(t1,s1,mu), timeSpan1, s1); % Need 3 args
[~, sol2] = ode45(@(t2,s2)diffEq(t2,s2,mu), timeSpan2, s2); % Need 3 arg
[~, solT] = ode45(@(tT,sT)diffEq(tT,sT,mu), timeSpanT, sT); % Need 3 arg
[~, solT2] = ode45(@(tT2,sT2)diffEq(tT2,sT2,mu), timeSpanT2, sT2); % Need 3 arg

[~, sol3] = ode45(@(t3,s3)diffEq(t3,s3,mu), timeSpan3, s3); % Need 3 arg
[~, solT3] = ode45(@(tT3,sT3)diffEq(tT3,sT3,mu), timeSpanT3, sT3); % Need 3 arg

% Plot Solutions
figure;

hold on; grid on; grid minor;
% Create Earth
[X, Y, Z] = sphere(); 

surf(X*rE, Y*rE, Z*rE)
plot3(sol0(:,1),sol0(:,2),sol0(:,3),'-r','LineWidth',1) % Plot in 3d
plot3(sol1(:,1),sol1(:,2),sol1(:,3),'-r','LineWidth',1) % Plot in 3d
plot3(solT(:,1),solT(:,2),solT(:,3),'-b','LineWidth',1) % Plot in 3d
plot3(solT2(:,1),solT2(:,2),solT2(:,3),'-b','LineWidth',1) % Plot in 3d
plot3(sol2(:,1),sol2(:,2),sol2(:,3),'-r','LineWidth',1) % Plot in 3d
plot3(sol3(:,1),sol3(:,2),sol3(:,3),'-r','LineWidth',1) % Plot in 3d
plot3(solT3(:,1),solT3(:,2),solT3(:,3),'-b','LineWidth',1) % Plot in 3d

%% Output solution
fprintf('Transfer Orbit DelV1: %.6f\n', delV1)
fprintf('Meanuever time: %.3f\n', tReq)

fprintf('Transfer Orbit DelV2: %.6f\n', delV2)
fprintf('Meanuever time: %.3f\n', tReq2)

fprintf('Transfer Orbit DelV3: %.6f\n', delV3)
fprintf('Meanuever time: %.3f\n', tReq3)

ttotal = tReq+tReq2+tReq3;
delVtot = delV1+delV2+delV3;
fprintf('Total Delta V: %.6f\n', delVtot)
fprintf('Total meanuever time: %.3f\n', ttotal)



%%


function sdot = diffEq(t0,s0,mu)

% r was first three elements, v was second three
r_0 = s0(1:3);
v_0 = s0(4:6);

sdot(1:3,1) = v_0; %First three elements are velocity
sdot(4:6,1) = (-mu*r_0)/(norm(r_0))^3; % Actual two body EOM


end


