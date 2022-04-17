%04/08/2015 A linearized version of IX isotherm

raw_y=[1.40 1.34 1.21 0.86 0.53 0.37];  %raw adsorbed ammonium values (meq/g)
raw_sod=[34.19 42.35 48.51 63.41 76.75 84.95]; %aqueous sodium data (meq/L)
raw_amm=[45.93 38.33 33.07 20.05 10.53 7.55]; %aqueous ammonium data (meq/L)

for i=1:6
    ydata(i)=1/raw_y(i); %recipricol of adsorbed ammonium values
end    

for i=1:6
    xdata(i)=raw_sod(i)/raw_amm(i); %aqeuous sod/aqeous amm
end    
           
p = polyfit(xdata,ydata,1);  
%ydata=[];   
xlim([0 5]);
ylim([0 3.5]);
plot(xdata,ydata,'x');
hold;
%plot(raw_y)
yfit = polyval(p,xdata);
yresid = ydata - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(ydata)-1) * var(ydata);
rsq = 1 - SSresid/SStotal;

for i=1:6
    fit(i)=p(1)*xdata(i)+p(2);
end
plot(xdata,fit);

Q=1/p(2);   %exchange capacity
K=1/(p(1)*Q); %affinity coefficient