function Bioregen_Model (t_max, dt, N_radial_intervals, surface_diffusion_coeff, NH4_initial, Na_initial,NO2_initial,NO3_initial)

% To Run type: Bioregen_Model(200,0.01,50,0.0019,80.13,24.7,0.27,0.007)
%  t_max -- time period for the simulation runs
%  dt -- time step used in finite difference method
%  N_radial_intervals -- spatial discretization of chabazite
%  surface_diffusion_coeff -- surface diffusion coefficient 
%  NH4_initial -- initial ammonium concentration
%  Na_initial -- initial sodium concentration
%  NO2_initial -- initial nitrite concentration
%  NO3_initial -- initial nitrate concentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Physical parameters -- from experimental isotherm data   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Q = 1.81;               % maximum sorption capacity of the chabazite [meq/g]
K = 2.92;               % isotherm parameter, value from equilibrium data
M_over_V = 150;         % g chabazite per L of solution  
d_grain = 1.5;          % grain diameter[mm]
r_grain = d_grain / 2;  % grain radius[mm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Biokinetic parameters              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
qmax_1=0.81/14;        %max specific growth rate AOB/yield coefficient AOB (fitted)
qmax_2=1.53/14;        %max specific growth rate NOB/yield coefficient NOB (fitted)
KS_1=5.0/14;           %half sat. constant ammonium (mg/L) (literature)
KS_2=0.9/14;           %half sat. constant nitrite (mg/L) (literature)
K_I1=1123.0/14;        %inhibition coefficient AOB (fitted)
K_I2=122/14;           %inhibition coefficient NOB (fitted) 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Short-cut notation
%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = r_grain / N_radial_intervals;
g = surface_diffusion_coeff * dt / dr^2 ;
h = (3../8.) * M_over_V * (1/N_radial_intervals)^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Build matrix of coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros (N_radial_intervals + 5);         % 2 additional nodes for NO2- and NO3- 
i = 1;                                      % corresponds to the very center of the chabazite
    A(i, i) = 1;
    A(i, i+1) = -1;
     
for i = 2 : N_radial_intervals              % interior chabazite nodes
    A(i, i-1) = g/(2*(i-1)) - g/2;
    A(i, i) = 1 + g;
    A(i, i+1) = -g/(2*(i-1)) - g/2;
end

i = N_radial_intervals + 1;                 % surface of the chabazite
    A(i, i) = 1;
    
i = N_radial_intervals + 2;                 % aqueous [NH4+]
    A(i, i-2) = -3*M_over_V*g/(2*N_radial_intervals);
    A(i, i-1) =  3*M_over_V*g/(2*N_radial_intervals);
    A(i, i) = 1;

i = N_radial_intervals + 3;                 % aqueous [Na+]
    j = 1;
        A(i, j) = h;
    for j = 2:N_radial_intervals
        A(i, j) = h * ((2*j-1)^2 + (2*j-3)^2);
    end
    j = N_radial_intervals + 1;
        A(i, j) = h * (2*j-3)^2;
        
    j = N_radial_intervals + 3;
        A(i, j) = -1;
      
i = N_radial_intervals + 4;              % aqueous [NO2-]
    A(i, i) = 1;        
        
i = N_radial_intervals + 5;              % aqueous [NO3-]
    A(i, i) = 1;        
        
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Initialize system
%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;
C_new = zeros([N_radial_intervals+5, 1]);
C_new(N_radial_intervals+2) = NH4_initial;                  
C_new(N_radial_intervals+3) = Na_initial;
C_new(N_radial_intervals+4) = NO2_initial;                  
C_new(N_radial_intervals+5) = NO3_initial;

counter = 1;
NH4_history(counter)  = C_new(N_radial_intervals + 2);
Na_history(counter)   = C_new(N_radial_intervals + 3);
NO2_history (counter) = C_new(N_radial_intervals + 4); 
NO3_history (counter) = C_new(N_radial_intervals + 5); 


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    March through time
%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros([N_radial_intervals+3, 1]);

while t < t_max
    
    t = t + dt;
    counter = counter + 1;
    C_old = C_new;
    
    %%%%    PREDICTOR LOOP
    i = 1;                                      % center of grain
        b(i) = 0;
        
    for i = 2:N_radial_intervals                % interior of grain
        b(i) = C_old(i-1) * (-g/(2*(i-1)) + g/2);
        b(i) = b(i) + C_old(i) * (1 - g);
        b(i) = b(i) + C_old(i+1) * (g/(2*(i-1)) + g/2);
    end
    
    i = N_radial_intervals + 1;                 % surface of grain
        NH4_conc = C_old(N_radial_intervals + 2);
        Na_conc = C_old(N_radial_intervals + 3);
        b(i) = Q*K*NH4_conc / (Na_conc + K*NH4_conc);
        
    i = N_radial_intervals + 2;                 % aqueous NH4+ conc.
        b(i) = C_old(i-2) * (3*M_over_V*g/(2*N_radial_intervals));
        b(i) = b(i) + C_old(i-1) * (-3*M_over_V*g/(2*N_radial_intervals));
        b(i) = b(i) + C_old(i) * (1);
        b(i) = b(i)-((dt*qmax_1*C_old(i))/((C_old(i)^2/K_I1)+ KS_1 + C_old(i)));
        
    i = N_radial_intervals + 3;                 % aqueous Na+ conc.
        b(i) = -Na_initial;
    
    i = N_radial_intervals + 4;                 % aqueous NO2- conc.
        b(i) = C_old(i);
        b(i) = b(i) + (dt*qmax_1*C_old(i-2))/((C_old(i-2)^2/K_I1)+ KS_1 + C_old(i-2));
        b(i) = b(i) - (dt*qmax_2*C_old(i))/((C_old(i-2)^2/K_I2)+ KS_2 + C_old(i)); 
        
            
    i = N_radial_intervals + 5;                 % aqueous NO3- conc.
        b(i) = C_old(i);
        b(i) = b(i) +(dt*qmax_2*C_old(i-1))/((C_old(i-3)^2/K_I2)+ KS_2 + C_old(i-1));
            
    %%%%    Solve for C_pred
    C_pred = A \ b;
    
    %%%%    Linearize
    C_avg = 0.5*(C_old + C_pred);
    
    %%%%    CORRECTOR LOOP
    %%%%    The adsorption isotherm has non-linearity in it
    i = N_radial_intervals + 1;
        NH4_conc = C_avg(N_radial_intervals + 2);
        Na_conc = C_avg(N_radial_intervals + 3);
        b(i) = Q*K*NH4_conc / (Na_conc + K*NH4_conc);
    
    %%%%    The biology also has non-linearity    
    i = N_radial_intervals + 2;                 % aqueous NH4+ conc.
        b(i) = C_avg(i-2) * (3*M_over_V*g/(2*N_radial_intervals));
        b(i) = b(i) + C_avg(i-1) * (-3*M_over_V*g/(2*N_radial_intervals));
        b(i) = b(i) + C_avg(i) * (1 -(dt*qmax_1*C_avg(i))/((C_avg(i)^2/K_I1)+ KS_1 + C_avg(i)));
    
    i = N_radial_intervals + 4;                 % aqueous NO2- conc.
        b(i) = C_avg(i);
        b(i) = b(i) + (dt*qmax_1*C_avg(i-2))/((C_avg(i-2)^2/K_I1)+ KS_1 + C_avg(i-2));
        b(i) = b(i) - (dt*qmax_2*C_avg(i))/((C_avg(i-2)^2/K_I2)+ KS_2 + C_avg(i)); 
        
            
    i = N_radial_intervals + 5;                 % aqueous NO3- conc.
        b(i) = C_avg(i);
        b(i) = b(i) +(dt*qmax_2*C_avg(i-1))/((C_avg(i-3)^2/K_I2)+ KS_2 + C_avg(i-1));    
        
    %%%%    Solve for C_new
    C_new = A \ b;
    NH4_history(counter) = C_new(N_radial_intervals + 2);
    Na_history(counter) =  C_new(N_radial_intervals + 3);
    NO2_history(counter) = C_new(N_radial_intervals + 4);
    NO3_history(counter) =  C_new(N_radial_intervals + 5);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plot histories predicted by model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len = length(NH4_history);
t_vec = 0 : dt : (len-1)*dt;
figure(1);
NO2_plot = plot(t_vec, NO2_history, 'k-', 'LineWidth', 2.3);
xlabel('Time (hr)', 'FontSize', 16);
ylabel('Aqueous Concentration (meq/L)', 'FontSize', 16);
hold on;

NO3_plot = plot(t_vec, NO3_history, 'color',[1 0.5 0], 'LineWidth', 2.3);
xlabel('Time (hr)', 'FontSize', 16);
ylabel('Aqueous Concentration (meq/L)', 'FontSize', 16);
axis([0.0 200.0 0.0 15.0]);

legend([NO2_plot,NO3_plot], 'Model {NO_2}^- (meq/L)', 'Model {NO_3}^- (meq/L)', 'Location', 'east');
set(gca,'fontsize',16)
set(gca,'FontWeight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plot histories predicted by model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (2);
len = length(NH4_history);
t_vec = 0 : dt : (len-1)*dt;
NH4_plot = plot(t_vec, NH4_history, 'k-', 'LineWidth', 2);
xlabel('Time (hr)', 'FontSize', 16);
ylabel('Aqueous Concentration (meq/L)', 'FontSize', 16);
hold on;
Na_plot = plot(t_vec, Na_history, 'k--', 'LineWidth', 2);
legend([NH4_plot,Na_plot], 'Model {NH_4}^+ (meq/L)', 'Model {Na}^+ (meq/L)', 'Location', 'east');
set(gca,'fontsize',16)
set(gca,'FontWeight','bold')


end