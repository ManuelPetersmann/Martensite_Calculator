a = 0.5^(1./3.)

eta = linspace(a,0.81);

dV = (eta.^3 .* 2 - 1.)*100;

plot(eta, dV);
xlabel('a_mart / a_aust')
ylabel('delta V [%]')

%% Marval
a_aust1 = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart1 = 2.8807346;
eta1 = a_mart1 / a_aust1;
dV2 = (eta1^3 * 2 - 1.)*100;
txt1 = ['\leftarrow ',num2str(dV2),' Marval: 12Cr-9Ni-2Mo-0.7Al-0.35Ti'];
text(eta1,dV2,txt1);

%% N404 Böhler
a_aust2 = 3.595; % 3.6017264; % for 200 Grad Celsius, 
a_mart2 = 2.885;  %2.8807346; % for 200 Grad Celsius, 
eta2 = a_mart2 / a_aust2;
dV2 = (eta2^3 * 2 - 1.)*100;
txt1 = ['\leftarrow ',num2str(dV2),' N404: 15Cr-5Ni-0.9Mo-0.4Si,Mn-0.04C'];
text(eta1,dV2,txt1);

%%