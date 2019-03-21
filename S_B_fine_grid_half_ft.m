%% This is a code for 2D non-isothermal reservoir simulation with compressible single phase fluid
% Using old values of pressure and temprature to calculate paramters 
% All the units are in field units 
close all
clear 
clc
global Rs gamma_g gamma_o API 
%% reservoir and fluid data
dx=30;                                       % length by ft of each grid in x direction  
dy=[256*ones(1,3*150) 64*ones(1,150) 16*ones(1,150) 0.5*ones(1,64*150) 16*ones(1,150) 64*ones(1,150) 256*ones(1,3*150)]'; % length by ft of each grid in y direction
dy_cave=0.5;                                    % grid cell width inside cavity
dz=15;                                        % length by ft of each grid in z direction
lx=150;                                       % number of grids in x direction 
N=length(dy);                                 % number of grids in y direction
API=10;                                       % API gravity of fluid
gamma_o=141.5/(API+131.5);                    % liquid specific gravity
rho_Liquid=62.4*gamma_o;                      % liquid density ib/cf 
Rs=30;                                        % Gas solubility SCF/STB
gamma_g=0.65;
rho_Solid=2.71*62.4;                          % density of rock by lb/ft3
Cp_Solid=0.2;                                 % specific heat of rock by BTU/lb.f
KL_cond=8.016;                                % fluid thermal conductivity by BTU/day.ft.R
KS_cond=61.537;                               % rock thermal conductivity by BTU/day.ft.R
ly=N/lx;                                      % total number of grids
K=[333*ones(((N-lx*48)/2),1); inf*ones(lx*48,1); 333*ones(((N-lx*48)/2),1)]; % horizontal permeability which is equal in x and y direction
phi=[0.18*ones(((N-lx*48)/2),1); ones(lx*48,1); 0.18*ones(((N-lx*48)/2),1)]; % porosity
Bc=1.127e-3; alpha=5.615; B2c = 1.0624E-14;  % conversion factors in DE
%% vertical wells calculations 
rw=0.33;            % raduis of vertical wells
iwell=[2 149];
jwell = [9 66];
mwell = iwell+(jwell-1)*lx;
req=0.14*(dx^2+(dy(mwell)).^2).^(0.5);    % equivalent raduis of vertical wells 
qsc=zeros(N,1);                              
qsc(mwell,1) = [-700 700]; % assining standard flow rate for each vertical completion
Tem_inj=670;   % injection temprature by R
%% initial conditions and time step 
P_initial=6000;       % initial pressure
Tem_ini=550;    % initial temprature
Fx_initial = 0;   % initial inter-grid mass rate in X-direction
Fy_initial = 0;   % initial inter-grid mass rate in Y-direction
tmax=500;      % simulation time days
dt=2;         % time step
AP = spalloc(N+(lx-1)*ly+(ly-1)*lx,N+(lx-1)*ly+(ly-1)*lx,9*N + 7*ly*(lx-1) + 7*lx*(ly-1)); % pressure rate matrix
AT = spalloc(N,N,5*N); % temprature matrix
%% Calculations
[L1,L2]=size(AP);
Vol=dx.*dy.*dz./(alpha*dt);
Ax=dy.*dz; % area in X-direction
Ay=dx.*dz; % area in Y-direction
Ax_x=Ax./dx; % area over length
Ay_y=Ay./dy;
RHS_P=zeros(N+(lx-1)*ly+(ly-1)*lx,1);
RHS_q=qsc;    % assigning values of wells rate inside RHS vector
Sol=zeros(2*N+(lx-1)*ly+(ly-1)*lx,(tmax/dt)+1); % solution matrix
Sol(:,1)=[P_initial*ones(N,1); Fx_initial*ones((lx-1)*ly,1); Fy_initial*ones((ly-1)*lx,1); Tem_ini*ones(N,1)];
Pwf=zeros((tmax/dt),length(mwell)); % bottom hole pressure 
for n = 2:((tmax/dt)+1)    % time loop
    P_old=Sol(1:N,n-1);     % grid old pressures
    T_old=Sol(L1+1:end,n-1); % grid old tempratures
    [phi_old,Der_phi] = phi_Der_phi(P_old,phi); % porosity at old tome level
    [Bo_old,mu_old,Cp_v_old] = mu_Bo_Cp(T_old,P_old);% calculating viscosity, formation volume factor and specific heat at old time level
    [der_Bo_P] = Der_Bo(T_old,P_old);
    bo_old=1./Bo_old; % shrinkage factor
    der_bo_P=-der_Bo_P.*(bo_old.^2); % derivative of shrinkage factor
    Cm=Vol.*(phi_old.*der_bo_P+bo_old.*Der_phi);  % Cm factor in conservation of mass
    Tx=Bc*K.*Ax_x.*bo_old./(mu_old); % transmissiblity in X-direction
    Ty=Bc*K.*Ay_y.*bo_old./(mu_old); % transmissiblity in Y-direction
    Cx=Bc*Ax_x.*bo_old./(B2c.*mu_old); % shear term in X-direction
    Cy=Bc*Ay_y.*bo_old./(B2c.*mu_old); % shear term in Y-direction
    for m=1:N
        j=ceil(m/lx);
        i=m-(j-1)*lx;
        if i > 1
            AP(m,N+(i-1)+(j-1)*(lx-1))=-1;  % Conservation of mass flow rate in from m-1 to m
        end
        if i<lx
            F=N+i+(j-1)*(lx-1);
            Cx_avg=0.5*((1/Cx(m))+(1/Cx(m+1))); % average harmonic shear term
            AP(m,F)=1;       % Conservation of mass flow rate out from m to m+1 
            AP(F,m)=-1;            % Conservation of momentuem pressure at grid m
            AP(F,m+1)=1;           % Conservation of momentuem pressure at grid m+1
            AP(F,F)=0.5*((1/Tx(m))+(1/Tx(m+1)))+2*Cx_avg*((1/(dx^2))+(1/(dy_cave^2)));%+(4/(dz^2))); % Conservation of momentuem Dary's part  
            if i<(lx-1)
                AP(F,F+1)=-Cx_avg/(dx^2);
            end
            if i>1
                AP(F,F-1)=-Cx_avg/(dx^2);
            end
            if j==1
                AP(F,F)=AP(F,F)+2*Cx_avg*(1/(dy_cave^2));
                AP(F,F+(lx-1))=-(4/3)*Cx_avg*(1/(dy_cave^2));
            elseif j==ly
                AP(F,F)=AP(F,F)+2*Cx_avg*(1/(dy_cave^2));
                AP(F,F-(lx-1))=-(4/3)*Cx_avg*(1/(dy_cave^2));
            else
                AP(F,F+(lx-1))=-Cx_avg*(1/(dy_cave^2));
                AP(F,F-(lx-1))=-Cx_avg*(1/(dy_cave^2));
            end
        end
        if j>1
            AP(m,N+ly*(lx-1)+m-lx)=-1;
        end
        if j<ly
            F=N+ly*(lx-1)+m;
            Cy_avg=0.5*((1/Cy(m))+(1/Cy(m+lx))); % average harmonic shear term
            AP(m,F)=1;
            AP(F,m)=-1;
            AP(F,m+lx)=1;
            AP(F,F)=0.5*((1/Ty(m))+(1/Ty(m+lx)))+2*Cy_avg*((1/(dx^2))+(1/(dy_cave^2)));%+(4/(dz^2)));
            if j<(ly-1)
                AP(F,F+lx)=-Cy_avg/(dy_cave^2);
            end
            if j>1
                AP(F,F-lx)=-Cy_avg/(dy_cave^2);
            end
            if i==1
                AP(F,F+1)=-(4/3)*Cy_avg*(1/(dx^2));
                AP(F,F)=AP(F,F)+2*Cy_avg*(1/(dx^2));
            elseif i==lx
                AP(F,F-1)=-(4/3)*Cy_avg*(1/(dx^2));
                AP(F,F)=AP(F,F)+2*Cy_avg*(1/(dx^2));
            else
                AP(F,F+1)=-Cy_avg*(1/(dx^2));
                AP(F,F-1)=-Cy_avg*(1/(dx^2));
            end
        end
        AP(m,m)=Cm(m);
    end
    RHS_P(1:N)=RHS_q+Cm.*P_old; % RHS for pressure system
    Sol(1:L1,n) = AP\RHS_P;     % solution for pressure and mass rate
    P_new=Sol(1:N,n);     % grid new pressures
    Fx=Sol(N+1:N+(lx-1)*ly,n);     % inter-grid new mass rate X-direction
    Fy=Sol(N+(lx-1)*ly+1:L1,n);   % inter-grid new mass rate Y-direction
    phi_new=phi_Der_phi(P_new,phi); % porosity at new pressure
    Bo_new=mu_Bo_Cp(T_old,P_new);% calculating formation volume factor new pressure
    bo_new=1./Bo_new; % shrinkage factor at new pressure
    Vol_T_new=Vol.*(phi_new.*bo_new.*Cp_v_old+(1-phi_new).*(rho_Solid/rho_Liquid)*Cp_Solid); % accumulation energy term at new pressure
    Vol_T_old=Vol.*(phi_old.*bo_old.*Cp_v_old+(1-phi_old).*(rho_Solid/rho_Liquid)*Cp_Solid); % accumulation energy term at old pressure
    Total_KT=(KL_cond.^phi_new).*(KS_cond.^(1-phi_new));  % Total thermal conductivity
    Thermal_KT_X=(Total_KT.*Ax_x./(alpha*rho_Liquid)); % Thermal transmissibility in X-direction
    Thermal_KT_Y=(Total_KT.*Ay_y./(alpha*rho_Liquid)); % Thermal transmissibility in Y-direction
    RHS_T=Vol_T_old.*T_old;
    for m=1:N % New temprature calculations
        j=ceil(m/lx);
        i=m-(j-1)*lx;
        AT(m,m)=Vol_T_new(m); % accumulation part
        if i>1
            F=m-1;
            F1=(i-1)+(j-1)*(lx-1);
            Tx_T_avg=2/((1/Thermal_KT_X(m))+(1/Thermal_KT_X(F)));
            AT(m,m)=AT(m,m)+Tx_T_avg; % conduction part
            AT(m,F)=-Tx_T_avg;
            if Fx(F1)<0 % Up winding % advection part
                Cp_avg=Cp_v_old(m);
                AT(m,m)=AT(m,m)-Cp_avg*Fx(F1);
            else
                Cp_avg=Cp_v_old(F);
                AT(m,F)=AT(m,F)-Cp_avg*Fx(F1);
            end
        end
        if i<lx
            F=m+1;
            F1=i+(j-1)*(lx-1);
            Tx_T_avg=2/((1/Thermal_KT_X(m))+(1/Thermal_KT_X(F)));
            AT(m,m)=AT(m,m)+Tx_T_avg; % conduction part
            AT(m,F)=-Tx_T_avg;
            if Fx(F1)>0
                Cp_avg=Cp_v_old(m);
                AT(m,m)=AT(m,m)+Cp_avg*Fx(F1);
            else
                Cp_avg=Cp_v_old(F);
                AT(m,F)=AT(m,F)+Cp_avg*Fx(F1);
            end
        end
        if j>1
            F=m-lx;
            F1=m-lx;
            Ty_T_avg=2/((1/Thermal_KT_Y(m))+(1/Thermal_KT_Y(F)));
            AT(m,m)=AT(m,m)+Ty_T_avg; % conduction part
            AT(m,F)=-Ty_T_avg;
            if Fy(F1)<0
                Cp_avg=Cp_v_old(m);
                AT(m,m)=AT(m,m)-Cp_avg*Fy(F1);
            else
                Cp_avg=Cp_v_old(F);
                AT(m,F)=AT(m,F)-Cp_avg*Fy(F1);
            end
        end
        if j<ly
            F=m+lx;
            F1=m;
            Ty_T_avg=2/((1/Thermal_KT_Y(m))+(1/Thermal_KT_Y(F)));
            AT(m,m)=AT(m,m)+Ty_T_avg;
            AT(m,F)=-Ty_T_avg;
            if Fy(F1)>0
                Cp_avg=Cp_v_old(m);
                AT(m,m)=AT(m,m)+Cp_avg*Fy(F1);
            else
                Cp_avg=Cp_v_old(F);
                AT(m,F)=AT(m,F)+Cp_avg*Fy(F1);
            end
        end
        for w=1:length(mwell)
            if m==mwell(w) && qsc(m)<0   % producers
                AT(m,m)=AT(m,m)-Cp_v_old(m)*qsc(m);
            elseif m==mwell(w) && qsc(m)>0   % injectors
                [~,~,Cp_inj] = mu_Bo_Cp(T_old(m),P_new(m));% calculating specific heat of injected fluid
                RHS_T(m)=RHS_T(m)+Cp_inj*qsc(m)*Tem_inj;
            end
        end
    end
    AT=AT./Tem_inj;
    RHS_T=RHS_T./Tem_inj;
    New_T = AT\RHS_T;     % solution for temprature
    New_T(New_T<Tem_ini)=Tem_ini;   % New_T(New_T<Tem_ini)=Tem_inj;
    Sol(L1+1:end,n)=New_T;
    
    [Bo_well,mu_well] = mu_Bo_Cp(New_T(mwell),P_new(mwell));% calculating viscosity, formation volume factor and specific heat at old time level
    J_well = qsc(mwell).*mu_well.*Bo_well.*(log(req/rw))./((7.08E-3) .* K(mwell) .* dz);
    Pwf(n-1,:)=P_new(mwell) + J_well ;    % wells bottom hole flowing pressure 
    
%     figure(n-1)
%     surf(vec2mat(New_T(:,1),lx));view(2);colormap jet; 
%     shading interp 
%     drawnow
end
%% ploting temprature distribution at certain time
x=0.5:1:(lx-0.50); y=[0.5:1:4.5 5.0625:0.125:12.9375 13.5:1:17.5];
[X,Y]=meshgrid(x,y);
figure(n)  % at 100 day
surf(X,Y,vec2mat(Sol(L1+1:end,(100/dt)+1),lx));view(2);colormap jet; 
shading interp 
figure(n+1)  % at 200 day
surf(X,Y,vec2mat(Sol(L1+1:end,(200/dt)+1),lx));view(2);colormap jet; 
shading interp 
%% ploting producer termprature
figure(n+2)
time=0:dt:tmax; Temprature=Sol(L1+mwell(1),:);
plot(time,Temprature)