%% Assignment 3 Monte-Carlo and Finite Difference Method Modelling
% Nathan Lavoy
% 100995612

%% Part One - Uniform Electric Field on Bottleneck
% Part one extends assignment one to add an acceleration due to a uniform
% electric field applied across the simulation. 

% Simulation parameters
m0 = 9.10938356e-31;    % Rest mass of electron(kg)
q0 = 1.60217653e-19;    % electron charge
m = 0.26*m0;            % Effective mass of elctrons
T = 300;                % Temperature(K)
k = 1.38064852e-23;     % Boltzmann constant
tmn = 0.2e-12;          % Mean time between collisions
vth = sqrt(2*k*T/m);    % Thermal Velocity in a 2D space 
lambda = vth*tmn;       % Mean free path
w = 200e-9;             % Simulation width
h = 100e-9;             % Simulation height
partNum = 500;         % Number of particles in simulation
plotNum = 10;           % Number of particles shown in plot
dt = h/vth/100;         % Time step 
iter = 500 ;            % Number of iterations
Vappx = 0.1;            % Applied Voltage X
Vappy = 0;            % Applied Voltage y
eConc = 10e15;

% Simulation state [width, height, xv, yv]
state = zeros(partNum, 4);
traj = zeros(iter,plotNum*2);
temp = zeros(iter,1);

% Random Maxwell-Boltzmann distribution
MB_pdf = makedist('Normal','mu', 0, 'sigma', sqrt(k*T/m));
pScat = 1 - exp(-dt/tmn);
totalCol = 0;    % Total number of collisions

% Initialize the elements
for i=1:partNum
    angle = 2*pi*rand;  % random generated angle
    state(i,:) = [w*rand h*rand random(MB_pdf) random(MB_pdf)];
end

% Electric Field X
Ex = Vappx/w;
Fx = Ex*q0;
ax = Fx/m0;

% Electric Field Y
Ey = Vappy/w;
Fy = Ey*q0;
ay = Fy/m0;

Etot = sqrt(Ex^2 + Ey^2);

fprintf('\nThe electric field is %d V/m.',Etot);
fprintf('\nThe electric force is %d N.\n',Fx);
fprintf('\nThe acceleration is %d m/s^2.\n',ax);
fprintf('Current Density = q0*electron conc*mu*Etot/(h*w)');
Vapp = 1.5;
% Electric Field X
Ex = Vappx/w;
Fx = Ex*q0;
ax = Fx/m0;

% Electric Field Y
Ey = Vappy/w;
Fy = Ey*q0;
ay = Fy/m0;

Etot = sqrt(Ex^2 + Ey^2);
% Main loop
for i = 1:iter
    %Calculate new displacement
    state(:,3) = state(:,3) + ax*dt;
    state(:,4) = state(:,4) + ay*dt;
    dx = dt*state(:,3);
    dy = dt*state(:,4);
    state(:,1) = state(:,1) + dx;
    state(:,2) = state(:,2) + dy;
    
    
    % Width check - should enter the other side with no loss
    % Leaving right
    check = state(:,1) > w;
    state(check,1) = 0;
    % Leaving left
    check = state(:,1) < 0;
    state(check,1) = w;
    % Height check - spectral reflection
    % Top
    check = state(:,2) > h;
    state(check,2) = h;
    state(check,4) = -state(check,4);
    % Bottom
    check = state(:,2) < 0;
    state(check,2) = 0;
    state(check,4) = -state(check,4);
    
    % Scatter particles using exponential scattering probability
    p = rand(partNum, 1) < pScat;
    state(p,3:4) = random(MB_pdf, [sum(p),2]);
    totalCol = totalCol + sum(p(:)==1);
    
    % Record the temperature
    temp(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/partNum;
    
    % Record the trajectories
    for j=1:plotNum
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    % Simulation two plot
    figure(1);
    subplot(2,1,1);
    hold off;
    % Plot current particle positions
    plot(state(1:plotNum,1)./1e-9, state(1:plotNum,2)./1e-9, 'o');
    axis([0 w/1e-9 0 h/1e-9]);
    title(sprintf('Part Two: %d of %d Electrons', plotNum, partNum));
    xlabel('x (nm)');
    ylabel('y (nm)');
    subplot(2,1,2);
    hold off;
    % Plot system temperature
    plot(dt*(0:i-1), temp(1:i));
    axis([0 dt*iter min(temp)*0.98 max(temp)*1.02]);
    title('System Temperature');
    xlabel('Time (s)');
    ylabel('Temperature (K)');
    % System Drift Current
    v = sqrt(state(:,3).^2 + state(:,4).^2);
    mu = sum(v)/partNum/Etot;
    DriftCur(i) = q0*eConc*mu*Etot/(h*w);
end
% Plot trajectories
figure(1);
title(sprintf('%d of %d Electrons Trajectories', plotNum, partNum));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 w/1e-9 0 h/1e-9]);
hold on;
for i=1:plotNum
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
end
hold off;
% Plot temperature
figure(2);
plot(dt*(0:iter-1), temp);
axis([0 dt*iter min(temp)*0.98 max(temp)*1.02]);
title('Semiconductor Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
% Plot Current
figure(4);
plot(dt*(0:iter-1), DriftCur);
axis([0 dt*iter min(DriftCur)*0.98 max(DriftCur)*1.02]);
title('Drift Current');
xlabel('Time (s)');
ylabel('Current (A)');

% Electron Density
density = hist3(state(:,1:2),[30 30])';
figure(5);
imagesc(density);
title('Electron Density');
axis off;
xlabel('x (nm)');
ylabel('y (nm)');
colorbar;

% Temperature calculation
tempDen = zeros(30,30);
partDen = zeros(30,30);
binX = 0:w/30:w;
binY = 0:h/30:h;
%Breakup space into grids and add a particles temp the grid it appears in
for i = 1:partNum
    for n = 1:30
        for p = 1:30
            if (state(i,1) > binX(n) && state(i,1) <= binX(n+1) && state(i,2) > binY(p) && state(i,2) <= binY(p+1))
                tempDen(n,p) = tempDen(n,p) + (state(i,3)^2 + state(i,4)^2)*m/k/2;
                partDen(n,p) = partDen(n,p) + 1;
            end
        end
    end
end
% Electron Density
figure(6);
imagesc((tempDen./partDen)');
title('Temperature Map');
axis off;
xlabel('x (nm)');
ylabel('y (nm)');
colorbar;
%% Part 2 - Distributed Electric Field in Bottle Neck
% Part 2 sets up the distributed electric field using FD method.
xlimleft = 0.8e-7;
xlimright = 1.2e-7;
ylimtop = 0.6e-7;
ylimbot = 0.4e-7;

xlim = 200e-9;
ylim = 100e-9;
nx = 200;
ny = 100;
Vapp = 1.5;
sigma1 = 1;
sigma2 = 0.0001;

[Vmap,Ex,Ey] = SolveBottleNeck(xlim,ylim,nx,ny,sigma1,sigma2,Vapp,xlimleft, xlimright, ylimtop, ylimbot);

figure(7)
subplot(2,1,1),H = surface(Vmap');
title('Voltage Map with Bottleneck')
set(H, 'linestyle', 'none');
view(45,45)
voltages=colorbar;
title(voltages,'Volts')

subplot(2,1,2)
quiver(Ex',Ey');
title('Electric Field from Potential - Vector Plot')
axis([0 nx 0 ny]);
%% Part 3 - MC and FD Combination
% Part 3 implements the distributed electric field in the Monte Carlo
% simulation. The position of the electrons are used to determine the
% acceleration they are influenced by.
% Simulation parameters
m0 = 9.10938356e-31;    % Rest mass of electron(kg).
q0 = 1.60217653e-19;    % electron charge
m = 0.26*m0;            % Effective mass of elctrons
T = 300;                % Temperature(K)
k = 1.38064852e-23;     % Boltzmann constant
tmn = 0.2e-12;          % Mean time between collisions
vth = sqrt(2*k*T/m);    % Thermal Velocity in a 2D space 
lambda = vth*tmn;       % Mean free path
w = 200e-9;             % Simulation width
h = 100e-9;             % Simulation height
partNum = 1000;         % Number of particles in simulation
plotNum = 10;           % Number of particles shown in plot
dt = h/vth/100;         % Time step 
iter = 1000;             % Number of iterations

% Simulation state [width, height, xv, yv]
state = zeros(partNum, 4);
traj = zeros(iter,plotNum*2);
temp = zeros(iter,1);

% Random Maxwell-Boltzmann distribution
MB_pdf = makedist('Normal','mu', 0, 'sigma', sqrt(k*T/m));
pScat = 1 - exp(-dt/tmn);
totalCol = 0;    % Total number of collisions

% Defines walls
boxes = 1e-9.*[80 120 0 40; 80 120 60 100];
% Fill simulation
for i = 1:partNum
    angle = rand*2*pi;
    state(i,:) = [w*rand h*rand random(MB_pdf) random(MB_pdf)];
    
    % remap illegal particle locations
    while(invalidPos(state(i,1), state(i,2), boxes))
        state(i,1:2) = [w*rand h*rand];
    end
end

% Electric Field X
Fx = Ex*q0;
ax = Fx/m0;

% Electric Field Y
Fy = Ey*q0;
ay = Fy/m0;
% Main loop
for i = 1:iter
    %Calculate new displacement
    dx = dt*state(:,3);
    dy = dt*state(:,4);
    state(:,1) = state(:,1) + dx;
    state(:,2) = state(:,2) + dy;
    % Top
    check = state(:,2) > h;
    state(check,2) = h;
    state(check,4) = -state(check,4);
    % Bottom
    check = state(:,2) < 0;
    state(check,2) = 0;
    state(check,4) = -state(check,4);
    % Left
    check = state(:,1) < 0;
    state(check,1) = 0;
    state(check,3) = -state(check,3);
    % Right
    check = state(:,1) > w;
    state(check,1) = w;
    state(check,3) = -state(check,3);
    % Walls
    for j = 1:partNum
        % if in box, determine which wall and reset particle
        if (invalidPos(state(j,1), state(j,2), boxes))
            xDist1 = abs(boxes(1,1)-state(j,1));    % Left wall
            xDist2 = abs(boxes(1,2)-state(j,1));    % Right wall
            yDist1 = abs(boxes(1,4)-state(j,2));    % Lower wall
            yDist2 = abs(boxes(2,3)-state(j,2));    % Upper wall
            wall = min([xDist1 xDist2 yDist1 yDist2]);
            if (xDist1 == wall)
                state(j,1) = boxes(1,1);
                state(j,3) = -state(j,3);
            elseif (xDist2 == wall)
                state(j,1) = boxes(1,2);
                state(j,3) = -state(j,3);
            elseif (yDist1 == wall)
                state(j,2) = boxes(1,4);
                state(j,4) = -state(j,4);
            elseif (yDist2 == wall)
                state(j,2) = boxes(2,3);
                state(j,4) = -state(j,4);   
            end
        end
        xpos = round(round(max(state(j,1)),9)*1e9);
        if (xpos == 0) xpos = xpos+1; end
        if (xpos > 200) xpos = xpos -1; end
        ypos = round(round(max(state(j,2)),9)*1e9);
        if (ypos == 0) ypos = ypos+1; end
        if (ypos > 100) ypos = ypos -1; end
        state(j,3) = state(j,3) + ax(xpos,ypos)*dt;
        state(j,4) = state(j,4) + ay(xpos,ypos)*dt;
    end
    % Scatter particles
    j = rand(partNum, 1) < pScat;
    state(j,3:4) = random(MB_pdf, [sum(j),2]);
    
    % Record temperatures
    temp(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/partNum;
    % Record positions for subset of particles that will be graphed
    for j=1:plotNum
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    % Simulation 3
    figure(8);
    hold off;
    plot(state(1:plotNum,1)./1e-9, state(1:plotNum,2)./1e-9, 'o');
    hold on;
    % Plot the boxes
    for j=1:size(boxes,1)
       plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,...
           [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
    end
    axis([0 w/1e-9 0 h/1e-9]);
    title(sprintf('Part Three: %d of %d Electrons', plotNum, partNum));
    xlabel('x (nm)');
    ylabel('y (nm)');
end
% Show trajectories
figure(8);
title(sprintf('Part Three: %d of %d Electrons Trajectories', plotNum, partNum));
xlabel('X (nm)');
ylabel('Y (nm)');
axis([0 w/1e-9 0 h/1e-9]);
hold on;
for i=1:plotNum
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
end

% Electron Density
density = hist3(state(:,1:2),[30 30])';
figure(9);
imagesc(density);
title('Electron Density');
axis off;
xlabel('x (nm)');
ylabel('y (nm)');
colorbar;

% Used to determine if particle is in the boxes
function val = invalidPos(x, y, boxes)
    val = false;
    if (x > boxes(1,1) && x < boxes(1,2) && y > boxes(1,3) && y < boxes(1,4))
        val = true;
    elseif (x > boxes(2,1) && x < boxes(2,2) && y > boxes(2,3) && y < boxes(2,4))
        val = true;
    end
end
%%
% The results of the simulation appear to support the theory. Electrons are
% shuffled through the gap and then into the right wall as directed by the
% electric field. The next step is to implement forces of electrons on each
% other due to charges rather than just collisions like standard particles. 