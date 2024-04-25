%Reading data and defining the function
dat = readtable('Data_FluxTunable.csv');
%dat = readtable('Data_Ideal.csv');
phi = dat.(1);
I_s = dat.(2);

Ic = max(I_s);

Bc = 10; %Value chosen for ease of observation
RC = .01;% Aribitrary value

h_cross = 1.05457182E-34;
e_ = 1.60217662E-19;

s_ = @(phi_in) interp1(phi,I_s,mod(phi_in,2*pi),"spline") / Ic;
%The above function uses interpolation to define the function at points
%not in the data file, and a modulo operator is used to ensure that the
%input variable is periodic with period 2*pi. And to be easy to use w.r.t
%the differential equation, we have made our definition such that s_(phi) =
%Is(phi)/Ic < 1.

%Now, we define our parameters which are the parts of the differential
%equation.

%Matrices A and B
A = [1 0 ; 0 Bc];
B = [0 1 ; 0 -1];

%A_inv = inv(A); This line is discarded since we can directly use A\

I_range = [linspace(0,4E-7,1000),linspace((4E-7)-(4E-10),-4E-7,1999),linspace(-(4E-7)+(4E-10),0,999)];

%Choosing a time-step: Note that this can depend on speed and computational
%resource availability, and hence a low value was used. 
t_span = [0,1e4]; %Note: The actual variable represented here is tau which is non-dimensionalized time.

dphis = [];

zeta_ini = [0,0]; %Initial conditions
zeta_c = zeta_ini;

for I = I_range
    diffz = @(t,z) A\(B*z(:,1) - [0 ; s_(z(1,1))] + [0 ; I/Ic]);
    [t_ , z_ ] = ode23s(diffz,t_span,zeta_c);
    zeta_c = z_(length(z_),:);
    dphis = [dphis zeta_c(2)];
end
V = dphis*h_cross*Bc/(2*e_*RC);

%Calculation of efficiency
[VcPlus,IcPlusInd] = max(diff(V));
IcPlus =  I_range(IcPlusInd);
[VcMinus,IcMinusInd] = min(diff(V));
IcMinus =  I_range(IcMinusInd);
eta = abs(IcPlus + IcMinus)/(abs(IcPlus) + abs(IcMinus));

%I_V Characteristics obtained, lets plot!
plot(V,I_range)
hold on
set(gca,'FontSize',20, 'FontName', 'Courier')
title('I-V Characterists - RCSJ Model')
xlabel('Voltage V (V)')
ylabel('Current I (A)')
str = "Efficiency = " + num2str(eta);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold off
