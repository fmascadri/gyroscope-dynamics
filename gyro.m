% Gyroscope dynamics

% This script calculates the decoupled equations of motion for the presented
% gyroscope system. The decoupled equations are converted into state-space form
% solved using an ODE solver for the parameters and initial conditions
% given and the results are plotted and animations generated.

% Further details on the modelling are contained in the discussion document.

clear all

%% Variable definitions

% Define symbolic variables
syms t alpha(t) beta(t) gamma(t) delta(t) % angles between frames

% angular velocities
alpha_dot = diff(alpha(t), t);
beta_dot = diff(beta(t), t);
gamma_dot = diff(gamma(t), t);
delta_dot = diff(delta(t), t);

% angular accelerations
alpha_ddot = diff(alpha_dot, t);
beta_ddot = diff(beta_dot, t);
gamma_ddot = diff(gamma_dot, t);
delta_ddot = diff(delta_dot, t);

% other symbols
syms r H h Ro Ri L  % geometry of gyroscope
syms m1 m2 m3 m4    % masses of four body components
syms g              % gravity acceleration

%% Rotation matrices
% Define rotation matrices - where R_a_b is the rotation from frame a to frame b
R_0_1 = [cos(alpha(t)) 0 -sin(alpha(t)); ...
         0 1 0; ...
         sin(alpha(t)) 0 cos(alpha(t))]; % neg rotation about y_0
R_1_2 = [1 0 0; ...
         0 cos(beta(t)) -sin(beta(t)); ...
         0 sin(beta(t)) cos(beta(t))];     % pos rotation about x_1
R_2_3 = [cos(gamma(t)) -sin(gamma(t)) 0; ...
         sin(gamma(t)) cos(gamma(t)) 0; ...
         0 0 1]; % pos rotation about z_2
R_3_4 = [cos(delta(t)) -sin(delta(t)) 0; ...
         sin(delta(t)) cos(delta(t)) 0; ...
         0 0 1]; % pos rotation about z_3

%% Tensors of inertia 
% Define tensors of interia for each body

% Frame part a - Thin central cylindrical rod OC in frame 3 about centre
% of gravity of frame G
I_3_a_G = m1*[(3*r^2+H^2)/12 0 0; ...
              0 (3*r^2+H^2)/12 0; ...
              0 0 0.5*r^2];

% Frame part b - Torus 'vertical ring' in frame 3 about centre of gravity of frame G
I_3_b_G = m2*[(3/4)*r^2 + Ro^2 0 0; ...
              0 (5/8)*r^2 + 0.5*Ro^2 0; ...
              0 0 (3/4)*r^2 + Ro^2];

% Frame part c - Torus 'horiztonal ring' in frame 3 about centre of gravity of frame G
I_3_c_G = m3*[(3/4)*r^2 + Ro^2 0 0; ...
              0 (3/4)*r^2 + Ro^2 0; ...
              0 0 (5/8)*r^2 + 0.5*Ro^2];

% Frame total
I_3_frame_G = I_3_a_G + I_3_b_G + I_3_c_G;

% Rotor solid cylinder in frame 4 about centre of gravity of rotor G
I_4_rotor_G = m4*[(3*Ri^2+h^2)/12 0 0; ...
                  0 (3*Ri^2+h^2)/12 0; ...
                  0 0 0.5*Ri^2];

%% Gravity forces
G_frame_0 = [0; 0; -(m1+m2+m3)*g];              % Weight force of Frame in frame 0
G_rotor_0 = [0; 0; -m4*g];                      % Weight force of Rotor in frame 0

%% Kinematic quantities
% Kinematic quantities for Frame of gyroscope
% -- Frame 1 --
omega_1_1 = [0; -alpha_dot; 0];
omega_1_1_dot = diff(omega_1_1);

r_1_O1_G1 = [0; 0; L];
r_1_O1_G1_dot = deriv(t,r_1_O1_G1,omega_1_1);
r_1_O1_G1_ddot = deriv(t,r_1_O1_G1_dot,omega_1_1);

% All zero valued in this case
r_1_O1_O2 = [0; 0; 0]; % same axes
r_1_O1_O2_dot = deriv(t,r_1_O1_O2,omega_1_1);
r_1_O1_O2_ddot = deriv(t,r_1_O1_O2_dot,omega_1_1);

% -- Frame 2 --
omega_2_21 = [beta_dot; 0; 0]; % i.e. velocity of frame 2, relative to frame 1
omega_2_2 = transpose(R_1_2)*omega_1_1 + omega_2_21;
omega_2_2_dot = transpose(R_1_2)*omega_1_1_dot + deriv(t,omega_2_21,omega_2_2);

r_2_O1_O2_ddot = transpose(R_1_2)*r_1_O1_O2_ddot;

r_2_O2_G2 = [0;0;L];
r_2_O2_G2_dot = deriv(t,r_2_O2_G2,omega_2_2);
r_2_O2_G2_ddot = deriv(t,r_2_O2_G2_dot,omega_2_2);

r_2_O1_G2_ddot = r_2_O1_O2_ddot + r_2_O2_G2_ddot;

% All zero valued in this case
r_2_O2_O3 = [0; 0; 0]; % same axes
r_2_O2_O3_dot = deriv(t,r_2_O2_O3,omega_2_2);
r_2_O2_O3_ddot = deriv(t,r_2_O2_O3_dot,omega_2_2);

% -- Frame 3 --
omega_3_32 = [0; 0; gamma_dot];
omega_3_3 = transpose(R_2_3)*omega_2_2 + omega_3_32;
omega_3_3_dot = transpose(R_2_3)*omega_2_2_dot + deriv(t,omega_3_32,omega_3_3);

r_3_O1_O3_ddot = transpose(R_2_3)*r_2_O1_O2_ddot;

r_3_O3_G3 = [0;0;L];
r_3_O3_G3_dot = deriv(t,r_3_O3_G3,omega_3_3);
r_3_O3_G3_ddot = deriv(t,r_3_O3_G3_dot,omega_3_3);

r_3_O1_G3_ddot = r_3_O1_O3_ddot + r_3_O3_G3_ddot;

r_3_O3_O4 = [0; 0; L]; % origin of frame 4, rotor, is at L above ground
r_3_O3_O4_dot = deriv(t,r_3_O3_O4,omega_3_3);
r_3_O3_O4_ddot = deriv(t,r_3_O3_O4_dot,omega_3_3);

% -- Frame 4 --
omega_4_43 = [0; 0; delta_dot];
omega_4_4 = transpose(R_3_4)*omega_3_3 + omega_4_43;
omega_4_4_dot = transpose(R_3_4)*omega_3_3_dot + deriv(t,omega_4_43,omega_4_4);

r_4_O1_O4_ddot = transpose(R_3_4)*r_3_O1_O3_ddot;

r_4_O4_G4 = [0;0;0]; % origin of frame 4 is coincident with G4
r_4_O4_G4_dot = deriv(t,r_4_O4_G4,omega_4_4);
r_4_O4_G4_ddot = deriv(t,r_4_O4_G4_dot,omega_4_4);

r_4_O1_G4_ddot = r_4_O1_O4_ddot + r_4_O4_G4_ddot;

%% Newton-Euler equations

% NE equations starting from Rotor
F_4_4 = m4*r_4_O1_G4_ddot - transpose(R_0_1*R_1_2*R_2_3*R_3_4)*G_rotor_0;
M_4_4 = cross(r_4_O4_G4, F_4_4) + I_4_rotor_G*omega_4_4_dot + cross(omega_4_4, I_4_rotor_G*omega_4_4);

% NE equations of Frame
F_3_3 = (m1+m2+m3)*r_3_O1_G3_ddot + R_3_4*F_4_4 - transpose(R_0_1*R_1_2*R_2_3)*G_frame_0;
M_3_3 = R_3_4*M_4_4 + cross(r_3_O3_G3, F_3_3) + ...
    I_3_frame_G*omega_3_3_dot + cross(omega_3_3, I_3_frame_G*omega_3_3);

%% EOM

T_rotor_z = M_4_4(3); % Torque from rotor acting in the z_4 axis
T_frame_x = M_3_3(1); % Torque from frame acting in the x_3 axis
T_frame_y = M_3_3(2); % Torque from frame acting in the y_3 axis
T_frame_z = M_3_3(3); % Torque from frame acting in the z_3 axis

%% Decouple

% Neaten (makes easier to visually inspect the terms in these large equations)

syms a b c d a_d b_d c_d d_d a_dd b_dd c_dd d_dd

% Note from R2022b, must specify higher order derivative terms first so
% that they aren't evaluated to 0 because of the lower order terms' symbols 
% being substituted first: (i.e. diff(alpha(t), t) => diff(a, t) = 0
old_symbols = [diff(alpha(t),t,t), diff(beta(t),t,t), diff(gamma(t),t,t), diff(delta(t),t,t),...
    diff(alpha(t),t), diff(beta(t),t), diff(gamma(t),t), diff(delta(t),t),...
    alpha(t), beta(t), gamma(t), delta(t)];

new_symbols = [a_dd, b_dd, c_dd, d_dd,...
    a_d, b_d, c_d, d_d,...
    a, b, c, d];

T_rotor_z_neat = subs(T_rotor_z, old_symbols, new_symbols);
T_frame_x_neat = subs(T_frame_x, old_symbols, new_symbols);
T_frame_y_neat = subs(T_frame_y, old_symbols, new_symbols);
T_frame_z_neat = subs(T_frame_z, old_symbols, new_symbols);

% Decouple
E1 = T_rotor_z_neat;
E2 = T_frame_x_neat;
E3 = T_frame_y_neat;
E4 = T_frame_z_neat;

[A,B] = equationsToMatrix([E1,E2,E3,E4],[a_dd,b_dd,c_dd,d_dd]);

decoupled = A\B;

a_dd = decoupled(1);
b_dd = decoupled(2);
c_dd = decoupled(3);
d_dd = decoupled(4);


%% Evaluate parameters

% Geometry of gyroscope
r = 0.002;  % cross-sectional radius of frame components [m]
H = 0.08;   % height of cylindrical rod [m]
h = 0.007;  % height of rotor [m]
Ro = 0.035; % major radius of frame torus [m];
Ri = 0.030; % radius of rotor [m]
L = 0.040;  % distance from O to G [m]

% Masses
m1 = 0.010;      % mass of cylindrical rod [kg]
m2 = 0.015;    % mass of torus [kg]
m3 = 0.015;    % mass of second torus [kg]
m4 = 0.07;     % mass of rotor [kg] 

% Gravitational acceleration
g = 9.81;   % [m/s^2]

% Substitute parameter values into the equations 
eq1 = subs(a_dd);
eq2 = subs(b_dd);
eq3 = subs(c_dd);
eq4 = subs(d_dd);

% Convert to variable precision from symbolic precision
eq1 = vpa(eq1);
eq2 = vpa(eq2);
eq3 = vpa(eq3);
eq4 = vpa(eq4);

%% Sub in x_n variables for state-space matrix form

eq1_xvars = subs(eq1, {a,b,c,d,...
                       a_d, b_d, c_d, d_d},...
                       {'x1','x2','x3','x4',...
                        'x5','x6','x7','x8'});
eq2_xvars = subs(eq2, {a,b,c,d,...
                       a_d, b_d, c_d, d_d},...
                       {'x1','x2','x3','x4',...
                        'x5','x6','x7','x8'});
eq3_xvars = subs(eq3, {a,b,c,d,...
                       a_d, b_d, c_d, d_d},...
                       {'x1','x2','x3','x4',...
                        'x5','x6','x7','x8'});
eq4_xvars = subs(eq4, {a,b,c,d,...
                       a_d, b_d, c_d, d_d},...
                       {'x1','x2','x3','x4',...
                        'x5','x6','x7','x8'});

% These four equations are substituted into the xdot state matrix for the
% functions of alpha_ddot, beta_ddot, gamma_ddot, delta_ddot


%% ODE solving

% Initial conditions
alpha_0 = deg2rad(5);
beta_0 = deg2rad(0);
gamma_0 = deg2rad(0);
delta_0 = 0;
alpha_dot_0 = deg2rad(0);
beta_dot_0 = deg2rad(0);
gamma_dot_0 = deg2rad(0);
delta_dot_0 = 300;  % rad/s
x0 = [alpha_0; beta_0; gamma_0; delta_0; alpha_dot_0; beta_dot_0; gamma_dot_0;
     delta_dot_0];

%  Time 
dt = 0.01; % seconds
tspan = [0 10]; % seconds
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% Calculation
x5_dotf = matlabFunction(eq1_xvars);
x6_dotf = matlabFunction(eq2_xvars);
x7_dotf = matlabFunction(eq3_xvars);
x8_dotf = matlabFunction(eq4_xvars);

sol = ode45(@(t,X) xdot(t,X,x5_dotf,x6_dotf,x7_dotf,x8_dotf), tspan, x0, options);
t = tspan(1):dt:tspan(2);
X = deval(sol,t);

%% Plot

PLOTTING = 0; % (1 if yes, 0 if no)

if PLOTTING
    figure
    
    subplot(8,1,1)
    plot(t,X(1, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\dot{\alpha}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,2)
    plot(t,X(2, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\dot{\beta}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,3)
    plot(t,X(3, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\dot{\gamma}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,4)
    plot(t,X(4, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\dot{\delta}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,5)
    plot(t,X(5, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\ddot{\alpha}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,6)
    plot(t,X(6, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\ddot{\beta}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,7)
    plot(t,X(7, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\ddot{\gamma}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
    
    subplot(8,1,8)
    plot(t,X(8, 1:end));
    xlabel('Time (t) [s]');
    ylabel('$\ddot{\delta}$', 'interpreter', 'latex', 'Rotation', 0, 'FontSize', 18);
end


%% Setup animation
% Set up animation objects

% Create thin central rod
    r = r;
    [x1,y1,z1] = cylinder(r,30);
    z1 = z1*H;

% Create first ring
    R = Ro; % outer radius of torus
    r = r; % inner tube radius
    th=linspace(0,2*pi,36); % e.g. 36 partitions along perimeter of the tube
    phi=linspace(0,2*pi,18); % e.g. 18 partitions along azimuth of torus
    % we convert our vectors phi and th to [n x n] matrices with meshgrid command:
    [Phi,Th]=meshgrid(phi,th);
    % now we generate n x n matrices for x,y,z according to eqn of torus
    x2=(R+r.*cos(Th)).*cos(Phi);
    y2=(R+r.*cos(Th)).*sin(Phi);
    z2=r.*sin(Th);
    z2 = z2 + L; % shift z coord up along L

% Create second ring
    R = Ro; % outer radius of torus
    r = r; % inner tube radius
    th=linspace(0,2*pi,36); % e.g. 36 partitions along perimeter of the tube
    phi=linspace(0,2*pi,18); % e.g. 18 partitions along azimuth of torus
    % we convert our vectors phi and th to [n x n] matrices with meshgrid command:
    [Phi,Th]=meshgrid(phi,th);
    % now we generate n x n matrices for x,y,z according to eqn of torus
    z3=(R+r.*cos(Th)).*cos(Phi);
    y3=(R+r.*cos(Th)).*sin(Phi);
    x3=r.*sin(Th);
    z3 = z3 + L; % shift z coord up along L

% Create rotor
    % rotor cylinder (only creates cylinder 'side wall', not end caps)
    R = Ri; % radius of rotor
    [x4,y4,z4] = cylinder(R,30);
    z4 = z4*h; % z4 * height of rotor
    z4 = z4+(L-h/2); % shift z coord up along L
   
    % bottom cap surface coords
    x5 = x4(1,:);
    y5 = y4(1,:);
    z5 = z4(1,:);
    [px5,py5,pz5] = create_pinwheel(x5,y5,z5);

    % top cap surface coords
    x6 = x4(end,:);
    y6 = y4(end,:);
    z6 = z4(end,:);
    [px6,py6,pz6] = create_pinwheel(x6,y6,z6);



%% Create animation

% Record video? (0 No, 1 Yes)
VIDEO = 1;

if VIDEO
    fps = 1/dt;
    MyVideo = VideoWriter(['gryo_', num2str(delta_dot_0)],'MPEG-4');
    MyVideo.FrameRate = fps;
    open(MyVideo);
end

handle = figure;
set(gcf, 'Position',  [100, 100, 1080, 1080])
hold on;

% top of gyro trace
top_x = ones(1,120)*100;
top_y = ones(1,120)*100;
top_z = ones(1,120)*100;

% Main animation loop
for i = 1:length(t)
    cla

    MANUAL_TESTING = 0;
    if MANUAL_TESTING
        % Manual inputs to test gyro configuration
        alpha = deg2rad(0);
        beta = deg2rad(5);
        gamma = deg2rad(0);
        delta = deg2rad(0);
    else
        % Values of angles from ODE solution
        alpha = X(1,i);
        beta = X(2,i);
        gamma = X(3,i);
        delta = X(4, i);
    end

    % Evaluate rotation matrices using angle values
    R_0_1 = [cos(alpha), 0, -sin(alpha);
            0, 1, 0;
            sin(alpha), 0, cos(alpha)];
    R_1_2 = [1, 0, 0;
            0, cos(beta), -sin(beta);
            0, sin(beta), cos(beta)];
    R_2_3 = [cos(gamma), -sin(gamma), 0;
            sin(gamma), cos(gamma), 0;
            0, 0, 1];
    R_3_4 = [cos(delta), -sin(delta), 0;
            sin(delta), cos(delta), 0;
            0, 0, 1];

    % Rotate gyro objects coordinates
    [x1_rotated,y1_rotated,z1_rotated] = rotation(x1,y1,z1,R_0_1*R_1_2*R_2_3);
    [x2_rotated,y2_rotated,z2_rotated] = rotation(x2,y2,z2,R_0_1*R_1_2*R_2_3);
    [x3_rotated,y3_rotated,z3_rotated] = rotation(x3,y3,z3,R_0_1*R_1_2*R_2_3);
    [x4_rotated,y4_rotated,z4_rotated] = rotation(x4,y4,z4,R_0_1*R_1_2*R_2_3*R_3_4);

    % Plot surfaces from rotated coordinates
    surf(x1_rotated,y1_rotated,z1_rotated, 'FaceColor', 'red', 'EdgeColor', 'none');
    surf(x2_rotated,y2_rotated,z2_rotated, 'FaceColor', 'blue', 'EdgeColor', 'none');
    surf(x3_rotated,y3_rotated,z3_rotated, 'FaceColor', 'green', 'EdgeColor', 'none');
    surf(x4_rotated,y4_rotated,z4_rotated, 'FaceColor','black', 'EdgeColor', 'none');

    % Rotor top and bottom surfaces (using patch)
    % colormap to create 'pinwheel' effect
    pinwheel_colormap = [ones(3,5) ones(3,5)*2 ones(3,5)*3 ones(3,5)*4 ones(3,5)*5 ones(3,5)*6];

    %bottom rotor surface:
    [px5_rotated,py5_rotated,pz5_rotated] = rotation(px5,py5,pz5,R_0_1*R_1_2*R_2_3*R_3_4);
    patch(px5_rotated,py5_rotated,pz5_rotated,pinwheel_colormap,'FaceColor', 'flat', 'EdgeColor', 'none');
    
    %top rotor surface:
    [px6_rotated,py6_rotated,pz6_rotated] = rotation(px6,py6,pz6,R_0_1*R_1_2*R_2_3*R_3_4);
    patch(px6_rotated,py6_rotated,pz6_rotated,pinwheel_colormap,'FaceColor', 'flat', 'EdgeColor', 'none');

    %top of gyro trace
    z_offset_trace = 0;
    top_x = [x3_rotated(end,end) top_x(1:end-1)];
    top_y = [y3_rotated(end,end) top_y(1:end-1)];
    top_z = [z3_rotated(end,end)+z_offset_trace top_z(1:end-1)];

    condition = top_x < 100; %filter out dummy values (for initial iterations of loop before enough of a trace is built up)

    plot3(top_x(condition),top_y(condition),top_z(condition),'LineWidth',2,'Color',[0,0,0,0.5]);

    % Final plot settings
    view(3)
    axis(0.06*[-1 1 -1 1 0 2])
    axis square
    set(gca,'position',[0.01 0.01 0.99 0.99])
    set(gca,'XTickLabel',"");
    set(gca,'YTickLabel',"");
    set(gca,'ZTickLabel',"");
    grid on;
    %title(['Gyroscope with rotor at ', num2str(delta_dot_0), ' rad/s'], 'FontSize', 18)

    if VIDEO
        writeVideo(MyVideo, getframe(handle));
    else
        pause(dt);
    end
end 

if VIDEO
    close(MyVideo);
end


%% Functions
function xdot = xdot(t,X,x5_dot,x6_dot,x7_dot,x8_dot)
    x1 = X(1);
    x2 = X(2);
    x3 = X(3);
    x4 = X(4);
    x5 = X(5);
    x6 = X(6);
    x7 = X(7);
    x8 = X(8);

    x1_dot = x5;
    x2_dot = x6;
    x3_dot = x7;
    x4_dot = x8;

    xdot =  [x1_dot; ...
             x2_dot; ...
             x3_dot; ...
             x4_dot; ...
             x5_dot(x1,x2,x3,x4,x5,x6,x7,x8); ... 
             x6_dot(x1,x2,x3,x4,x5,x6,x7,x8); ...
             x7_dot(x1,x2,x3,x4,x5,x6,x7,x8); ...
             x8_dot(x2,x3,x5,x6)];
end

function [Xf,Yf,Zf] = rotation(Xi,Yi,Zi,R)
    I=size(Xi,1);
    J=size(Xi,2);

    Xf=zeros(I,J);
    Yf=zeros(I,J);
    Zf=zeros(I,J);

    for ii=1:I
        for jj=1:J
            vector=[Xi(ii,jj);Yi(ii,jj);Zi(ii,jj)];
            vector=R*vector;

            Xf(ii,jj)=vector(1);
            Yf(ii,jj)=vector(2);
            Zf(ii,jj)=vector(3);
        end
    end
end

function v_dot = deriv(t, v_f, omega)
    v_dot = diff(v_f, t) + cross(omega,v_f); 
end

function [X,Y,Z] = create_pinwheel(x,y,z)
    % Takes the array of points along the circle of the rotor and returns
    % new arrays with points added at the origin of the circle to create
    % the pinwheel wedges which can then be colored in using a colormap.
    
    X = zeros(3,size(x,2)-1);
    Y = zeros(3,size(x,2)-1);
    Z = zeros(3,size(x,2)-1);
    C = zeros(1,size(x,2)-1);

    for i = 1:(size(x,2)-1)
        X(1,i) = x(i);
        X(2,i) = x(i+1);
        X(3,i) = 0;     % adds point at origin to create the pinwheel wedge

        Y(1,i) = y(i);
        Y(2,i) = y(i+1);
        Y(3,i) = 0;     % adds point at origin to create the pinwheel wedge

        % Where this is called in the code, 
        % the z dimension is constant. This would need
        % to be updated if that changed
        Z(1,i) = z(i);
        Z(2,i) = z(i);
        Z(3,i) = z(i);  
    end    
end
