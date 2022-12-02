%Script that approximates the side and drag forces, the average rotor 
%torque acting on a rotor disk, and the power requirement to drive the
%rotor for a rotor disk of user-specified geometry travelling at a 
%user-specified air velocity with vertical and forward components.
%Written by Jonathan De Sousa
%Credit to Dr Errikos Levis for cubicEval.m function
%Date: 08-04-2019

%Clear workspace and screen and close all open figures
clear
clc
close all

%Prompt user to input operating conditions
omega=input('Rotor angular velocity (rad/s): ');
i_c=input('Collective blade pitch (degrees): '); 
V_inf=input('Forward airspeed (knots): '); 
W_inf=input('Rate of descent (ft/min): '); 
altitude=input('Altitude (m): ');%geopotential height above sea level 
[t,a,p,rho]=atmosisa(altitude); %obtain density using ISA model

clc %clear screen

%Prompt user to define rotor geometry
N=input('Number of blades: '); 
D=input('Rotor diameter (m): ');
R_0=input('First aerofoil section (m): ');
c_root=input('Root chord length (m): ');
c_tip=input('Tip chord length (m): ');
epsilon=input('Blade twist angle (degrees): ');

clc %clear screen

%Prompt user for full file name of ASCII file containing lift and drag 
%coefficients at different angles of attack.
filename=input('Full name of data file (e.g. N0012.dat): ','s');

%Load file
fID=fopen(char(filename),'r');

%Read file
fscanf(fID,'%s',3); %skip first three lines (table headings)
data=fscanf(fID,'%f',[3,inf]); %extract data from file

fclose(fID); %close file

%Transpose data array and separate data into categories
data=data';
AoA=data(:,1); %angles of attack
CL=data(:,2); %corresponding coefficients of lift
CD=data(:,3); %corresponding coefficients of drag

%Change all angles into radians, and velocities into m/s
AoA=deg2rad(AoA);
i_c=deg2rad(i_c);
epsilon=deg2rad(epsilon);
V_inf=convvel(V_inf,'kts','m/s');
W_inf=convvel(W_inf,'ft/min','m/s');

%Discretize rotor disk
m_R=200; %number of discrete rotor sections 
R=(R_0:(D/2-R_0)/(m_R-1):D/2); %distances of each section from rotor disc centre (m)
m_azimuth=360; %number of discretized azimuth angles
azimuth=(0:2*pi/(m_azimuth-1):2*pi); %azimuth angles of the blade (rad)

%Calculate local chord length for each discretized section using linear
%interpolation
c_local=(((c_tip-c_root)/(D/2-R_0)).*(R-R_0)+c_root)';

%Set constants
Fn_old=0; %initial thrust force of rotor disc (N)
diff=99; %difference between successive calculated total thrust values
const=N/(2*pi); %constant used in all force and moment calculations

%Pre-allocate arrays
V_e=zeros(m_R,m_azimuth);
delta_alpha=zeros(m_R,m_azimuth);
CL_InterpVal=zeros(m_R,1);
CD_InterpVal=zeros(m_R,1);

%Create meshgrid for azimuth and R values
[azimuth,R]=meshgrid(azimuth,R);

%Calculate polynomial coefficients of interpolant function to be used for
%linear interpolation for lift and drag coefficients inside while loop
polycoeff_CL=coeff(AoA,CL);
polycoeff_CD=coeff(AoA,CD);

%Calculate tangential velocity to blade in rotor disc plane (m/s)
V_T=omega*R+V_inf*sin(azimuth);

angle=i_c+((R-R_0)/(D/2-R_0))*epsilon; %constant angle repeatedly used in loop

%Loop iterative calculation of total thrust force, Fn, until Fn converges
%to within an precision of 10^-5 
while diff>10^-5
    %Calculate total air velocity normal to the rotor plane (m/s)
            if W_inf<=0
                W=(W_inf/2)-0.5*sqrt(W_inf^2+8*Fn_old/(rho*pi*D^2));
            elseif W_inf>sqrt(8*Fn_old/(rho*pi*D^2))
                W=(W_inf/2)+0.5*sqrt(W_inf^2-8*Fn_old/(rho*pi*D^2));
            else
                %Display warning if rotorcraft is in Vortex Ring State,
                %i.e   0 < W_inf <= sqrt((8*Fn)/(rho*pi*D^2))
                disp('WARNING: ROTORCRAFT IS IN DANGEROUS VORTEX RING STATE')
            end
            
            %Calculate local airspeed in plane of blade element (m/s)
            V_e=sqrt(V_T.^2+W^2);
            %Calculate local rotation of the airflow relative to the plane of
            %the disk (rad)
            delta_alpha=atan(W./V_T);
            %Calculate angle of attack (rad) of blade element
            alpha_e=angle+delta_alpha;
            %Use cubic spline interpolation to find coefficients of drag and
            %lift at each section using cubicEval.m and coeff.m functions used
            %in Cubic Spline Interpolation question
            for i=1:m_R
                CL_InterpVal(i,1)=cubicEval(AoA,polycoeff_CL,alpha_e(i,1)); %interpolated CL value
                CD_InterpVal(i,1)=cubicEval(AoA,polycoeff_CD,alpha_e(i,1)); %interpolated CD value
            end
            
            %Calculate total thrust force Fn_new of the rotor disc
            dFn=0.5*rho*V_e.^2.*c_local.*(CL_InterpVal.*cos(delta_alpha)+CD_InterpVal.*sin(delta_alpha));
            Fn_new=const*trapz((azimuth(1,:))',trapz(R(:,1),dFn));
            
            %Calculate difference between successive thrust force values
            diff=abs(Fn_new-Fn_old);
            Fn_old=Fn_new;
end

%Calculate the force generated by each section in the plane of the rotor disk
dFp=0.5*rho*V_e.^2.*c_local.*(CD_InterpVal.*cos(delta_alpha)-CL_InterpVal.*sin(delta_alpha));

%Calculate Forces, Moments and Power required to drive rotor disk
Fx=const*trapz((azimuth(1,:)),trapz(R(:,1),dFp.*sin(azimuth))); %drag force
Fy=-const*trapz((azimuth(1,:)),trapz(R(:,1),dFp.*cos(azimuth))); %side force
Mx=-const*trapz((azimuth(1,:)),trapz(R(:,1),dFn.*R.*cos(azimuth))); %rolling moment
My=-const*trapz((azimuth(1,:)),trapz(R(:,1),dFn.*R.*sin(azimuth))); %pitching moment
T=const*trapz((azimuth(1,:)),trapz(R(:,1),dFp.*R));%average torque over the rotor
P=T*omega; %power P required to drive the rotor

clc %clear screen

%Display results
fprintf(1,'Total thrust force.. = %.3e N \n',Fn_new);
fprintf(1,'Power to drive rotor = %.3e W \n',P);
fprintf(1,'Drag force.......... = %.3e N \n',Fx);
fprintf(1,'Side force.......... = %.3e N \n',Fy);
fprintf(1,'Rolling moment...... = %.3e Nm \n',Mx);
fprintf(1,'Pitching moment..... = %.3e Nm \n',My);
fprintf(1,'Rotor Torque........ = %.3e Nm \n',T);

%end of script