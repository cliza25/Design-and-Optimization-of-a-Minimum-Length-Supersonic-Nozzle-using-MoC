%% Minimum Length Nozzle Design using Method of Characteristics
% Converted from Python to MATLAB
% Gamma = 1.3
% Exit Mach Number = 3.0

clear;
clc;
close all;

%% ---------------- Problem Data ----------------

gamma = 1.3;
Me = 3.0;

n = 20;                 % Number of characteristic lines
ht = 20.0;              % Throat height (mm)
ya0 = ht/2;

%% ------------ Prandtl-Meyer quantities ------------

nu_e = nu_of_M(Me,gamma);

theta_max = nu_e/2;

dtheta = theta_max/n;

fprintf('(a) nu_e = %.4f deg\n',nu_e);
fprintf('(b) theta_max = %.4f deg\n',theta_max);
fprintf('dtheta = %.4f deg\n',dtheta);

%% ------------ Total mesh points ------------

Ntot = n*(n+3)/2;

fprintf('Total mesh points = %d\n',Ntot);

%% ------------ Allocate Arrays ------------

theta = zeros(Ntot,1);
nu = zeros(Ntot,1);

Mach = zeros(Ntot,1);
mu = zeros(Ntot,1);

Kminus = zeros(Ntot,1);
Kplus = zeros(Ntot,1);

x = zeros(Ntot,1);
y = zeros(Ntot,1);

rowOf = zeros(Ntot,1);
lineOf = zeros(Ntot,1);

onaxis = false(Ntot,1);
onwall = false(Ntot,1);

idx = -ones(n+1,n+1);

wallidx = -ones(n+1,1);

%% =====================================================
%% Build theta, nu, compatibility constants
%% =====================================================

pt = 0;

for r = 1:n

    if r==1
        Lrange = 1:n;
    else
        Lrange = r:n;
    end

    for L = Lrange

        pt = pt + 1;

        if r==1

            th = L*dtheta;
            nuv = L*dtheta;

        elseif L==r

            th = 0;

            nuv = 2*r*dtheta;

            onaxis(pt)=true;

        else

            th = (L-r)*dtheta;

            nuv = (L+r)*dtheta;

        end

        theta(pt)=th;
        nu(pt)=nuv;

        Kminus(pt)=th+nuv;
        Kplus(pt)=th-nuv;

        rowOf(pt)=r;
        lineOf(pt)=L;

        idx(L,r)=pt;

    end

    pt=pt+1;

    lastpt=idx(n,r);

    theta(pt)=theta(lastpt);

    nu(pt)=nu(lastpt);

    Kminus(pt)=Kminus(lastpt);

    Kplus(pt)=Kplus(lastpt);

    onwall(pt)=true;

    rowOf(pt)=r;

    lineOf(pt)=n;

    wallidx(r)=pt;

end

%% =====================================================
%% Compute Mach number and Mach angle
%% =====================================================

for i=1:Ntot

    Mach(i)=M_of_nu(nu(i),gamma);

    mu(i)=asind(1/Mach(i));

end

%% =====================================================
%% Geometry
%% =====================================================

xa=0;

ya=ya0;

theta_a_flow=0;

mu_a_flow=90;

theta_a_wall=theta_max;

for r=1:n

    if r==1
        Lrange=1:n;
    else
        Lrange=r:n;
    end

    for L=Lrange

        p=idx(L,r);

        %% Previous point on same characteristic

        if r==1

            xA=xa;

            yA=ya;

            thA=theta_a_flow;

            muA=mu_a_flow;

        else

            A=idx(L,r-1);

            xA=x(A);

            yA=y(A);

            thA=theta(A);

            muA=mu(A);

        end

        %% Axis point

        if L==r

            thC=theta(p);

            muC=mu(p);

            slope=tand(0.5*((thA-muA)+(thC-muC)));

            xC=xA-yA/slope;

            yC=0;

        else

            %% Interior point

            B=idx(L-1,r);

            xB=x(B);

            yB=y(B);

            thB=theta(B);

            muB=mu(B);

            thC=theta(p);

            muC=mu(p);

            slopeA=tand(0.5*((thA-muA)+(thC-muC)));

            slopeB=tand(0.5*((thB+muB)+(thC+muC)));

            % Avoid division by zero
            if abs(slopeA-slopeB)<1e-12
                slopeA=slopeA+1e-12;
            end

            xC=((yB-slopeB*xB)-(yA-slopeA*xA))/(slopeA-slopeB);

            yC=yA+slopeA*(xC-xA);

        end

        x(p)=xC;

        y(p)=yC;

    end
    %% =====================================================
%% Compute Wall Points
%% =====================================================

    wpt = wallidx(r);
    lastpt = idx(n,r);

    % Previous wall point
    if r==1

        xW = xa;
        yW = ya;
        thW = theta_a_wall;

    else

        xW = x(wallidx(r-1));
        yW = y(wallidx(r-1));
        thW = theta(wallidx(r-1));

    end

    % Last interior point
    xA = x(lastpt);
    yA = y(lastpt);
    thA = theta(lastpt);
    muA = mu(lastpt);

    % Wall slope
    slopeWall = tand(0.5*(thW+thA));

    % C+ characteristic slope
    slopeC = tand(thA+muA);

    if abs(slopeWall-slopeC) < 1e-12
        slopeWall = slopeWall + 1e-12;
    end

    xN = (yA-yW-slopeC*xA+slopeWall*xW)/(slopeWall-slopeC);

    yN = yW + slopeWall*(xN-xW);

    x(wpt)=xN;
    y(wpt)=yN;

end

%% =====================================================
%% Exit Geometry
%% =====================================================

L_nozzle = x(wallidx(n));

he = y(wallidx(n));

fprintf('\n');
fprintf('Nozzle Length = %.5f mm\n',L_nozzle);

fprintf('Exit Half Height = %.5f mm\n',he);

fprintf('Exit Full Height = %.5f mm\n',2*he);

%% =====================================================
%% Isentropic Area Ratio
%% =====================================================

AeAstar = (1/Me)*...
    ((2/(gamma+1))*...
    (1+((gamma-1)/2)*Me^2))^...
    ((gamma+1)/(2*(gamma-1)));

fprintf('Area Ratio (Isentropic) = %.5f\n',AeAstar);

fprintf('Area Ratio (Geometry)   = %.5f\n',he/ya0);

%% =====================================================
%% Export Wall Coordinates
%% =====================================================

wallData = zeros(n+1,5);

wallData(1,:) = [0 xa ya theta_a_wall 1];

for r=1:n

    w = wallidx(r);

    wallData(r+1,:) = ...
        [r ...
         x(w) ...
         y(w) ...
         theta(w) ...
         Mach(w)];

end

header = {'Rowid','x_mm','y_mm','Theta_deg','Mach'};

T = array2table(wallData,...
    'VariableNames',header);

writetable(T,'wall_table.csv');

fprintf('\nWall coordinates saved to wall_table.csv\n');

%% =====================================================
%% Plot Characteristic Net
%% =====================================================

figure('Color','w');

hold on;

%---------------------------------------
% C- Characteristics
%---------------------------------------

for L=1:n

    xs = xa;
    ys = ya;

    for r=1:min(L,n)

        p = idx(L,r);

        xs = [xs x(p)];
        ys = [ys y(p)];

    end

    plot(xs,ys,'b','LineWidth',0.6);

end

%---------------------------------------
% C+ Characteristics
%---------------------------------------

for r=1:n

    xs = [];
    ys = [];

    for L=r:n

        p = idx(L,r);

        xs = [xs x(p)];
        ys = [ys y(p)];

    end

    xs = [xs x(wallidx(r))];
    ys = [ys y(wallidx(r))];

    plot(xs,ys,'r','LineWidth',0.6);

end

%% =====================================================
%% Plot Nozzle Wall
%% =====================================================

wx = xa;
wy = ya;

for r=1:n

    wx = [wx x(wallidx(r))];
    wy = [wy y(wallidx(r))];

end

plot(wx,wy,'k','LineWidth',2);

plot(wx,-wy,'k','LineWidth',2);

yline(0,'k--');

xlabel('x (mm)');

ylabel('y (mm)');

title('Characteristic Net');

axis equal;

grid on;

saveas(gcf,'characteristic_net.png');
%% =====================================================
%% Nozzle Contour Plot
%% =====================================================

figure('Color','w');

hold on;

plot(wx,wy,'k','LineWidth',2);
plot(wx,-wy,'k','LineWidth',2);

% Fill nozzle interior
fill([wx fliplr(wx)],...
     [wy fliplr(-wy)],...
     [0.85 0.92 1.00],...
     'EdgeColor','none');

% Draw contour again
plot(wx,wy,'k','LineWidth',2);
plot(wx,-wy,'k','LineWidth',2);

yline(0,'k--');

xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Minimum Length Nozzle (M_e = %.1f)',Me));

axis equal;
grid on;
box on;

saveas(gcf,'nozzle_contour.png');

fprintf('\n');
fprintf('Nozzle contour saved as nozzle_contour.png\n');

%% =====================================================
%% Print Wall Coordinate Table
%% =====================================================

fprintf('\n');
fprintf('-----------------------------------------------\n');
fprintf(' Wall Point Coordinates\n');
fprintf('-----------------------------------------------\n');
fprintf(' Row      x(mm)      y(mm)   Theta(deg)     Mach\n');
fprintf('-----------------------------------------------\n');

fprintf('%4d %12.5f %12.5f %12.5f %10.5f\n',...
        0,xa,ya,theta_a_wall,1.0);

for r = 1:n

    w = wallidx(r);

    fprintf('%4d %12.5f %12.5f %12.5f %10.5f\n',...
        r,...
        x(w),...
        y(w),...
        theta(w),...
        Mach(w));

end

fprintf('-----------------------------------------------\n');

%% =====================================================
%% Summary
%% =====================================================

fprintf('\n');
fprintf('=============== RESULTS ===============\n');
fprintf('Gamma                 : %.3f\n',gamma);
fprintf('Exit Mach Number      : %.3f\n',Me);
fprintf('Characteristic Lines  : %d\n',n);
fprintf('Maximum Wall Angle    : %.4f deg\n',theta_max);
fprintf('Prandtl-Meyer Angle   : %.4f deg\n',nu_e);
fprintf('Nozzle Length         : %.4f mm\n',L_nozzle);
fprintf('Exit Half Height      : %.4f mm\n',he);
fprintf('Exit Full Height      : %.4f mm\n',2*he);
fprintf('Geometric Area Ratio  : %.5f\n',he/ya0);
fprintf('Isentropic Area Ratio : %.5f\n',AeAstar);
fprintf('=======================================\n');

disp(' ');
disp('Program Completed Successfully.');