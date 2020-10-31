% Script for plotting IX isotherm 
% and calculating NH4+ removal
% as a function of chabazite dosage

% Experimental ion concentrations
obs_ammconc = [45.93 38.33 33.07 20.05 10.53 7.55 0.0];  %meq/L
obs_adsamm  = [1.4 1.34 1.21 0.86 0.53 0.37 0.0];        %meq/g
obs_sodconc = [34.19 42.35 48.51 63.41 76.75 84.95];     %meq/L  

init_Q = 1.81;
init_K = 2.92;
err=[0.07 0.06 0.03 0.00 0.01 0.02 0.0];

subplot(2,2,[1 2]);
errorbar(obs_ammconc,obs_adsamm,err,'ko','Markersize',12,'Linewidth',2.0,'MarkerFaceColor', 'k');
title ('A');
ylabel('Adsorbed NH_4^+ Concentration (meq/g)','FontSize', 12);
xlabel('Aqueous NH_4^+ Concentration  (meq/L)','FontSize', 12);
hold on;

a=(0.0:0.001:50.0);
varx=sort(a,'descend');

b=(30.0:0.001:90.0);
varsod=b;

for i=1:size(varx,2)
   model(i)=(init_Q*init_K*varx(i))/(varsod(i)+(init_K*varx(i)));
end

plot(varx,model,'k-','Linewidth',2.3);
axis([0.0 50.0 0.0 2.0]);
hold on;
h = legend('Experimental','Model');
set(h, 'Box', 'on') ;
set(h,'Location','best');

hold on;
set(gca,'fontsize',12);
set(gca,'FontWeight','bold');
set(gca,'XTick',0:5:50);
set(gca,'YTick',0:0.2:2);
set(gca,'LineWidth',2.5);

%%%%%%%%%%%%%%%%%%%%
%Start of new plot %
%%%%%%%%%%%%%%%%%%%%

%Plot of NH4 Removal vs Chabazite dose

V = 1;    %volume of liquid
C_0 = 65.0;
Na_0 = 35.0;
Q = 1.81;
K = 2.92;

%Experimental removal data
dose = [0 13 19 25 50 100 150];  %g/L
y_data = [0 27.7 39.6 47.9 68.4 83.4 88.1]; 
err=[0.0 1.93 1.19 0.19 1.41 0.35 3.67];

subplot(2,2,[3 4]);
errorbar(dose,y_data,err,'ko','Markersize',12,'Linewidth',2.0,'MarkerFaceColor', 'w');
title ('B');

ylabel('NH_4^+ Removal (%)','FontSize', 12);
xlabel('Chabazite Dose (g chabazite/L)','FontSize', 12);
set(gca,'fontsize',12)
set(gca,'FontWeight','bold')
set(gca,'XTick',0:5:50);
set(gca,'YTick',0:0.2:2);
set(gca,'LineWidth',2.5);
axis([0 155 0 95]);

hold on;

M = linspace(0.0,155.0);
a = C_0*(K-1);
c = -(C_0+Na_0);
denom = 2*C_0*(K-1); 
term1 = 4*(C_0+Na_0)*C_0*(K-1);

for i=1:size(M,2)
    b(i)=(((M(i)/V)*Q*K)+(C_0+Na_0)-C_0*(K-1));
    model(i)=(-b(i)+sqrt(b(i)^2+term1))/(denom);
    removal(i)=100*(1-model(i));
end

plot(M,removal,'k-','Linewidth',2.3);
h2 = legend('Experimental','Model');
set(h2, 'Box', 'on') ;
set(gca,'fontsize',12)
set(gca,'XTick',0:25:150);
set(gca,'YTick',0:10:100);
set(gca,'FontWeight','bold')
set(h2,'Location','best');

