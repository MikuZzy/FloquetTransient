% This file generates synthetic experimental data by simulating Smolen
% (2004) model. A temperature change is simulated at t = 240hr by changing
% reaction parameters. Simulated data between t1 = 120hr to t2 = 360hr are
% collected to be comparable to experimental data.

% There are three conditions: 1) transient change of lambda = infinity; 2)
% fast change of lambda = 0.5; 3) slow change of lambda = 0.05.

% Model Parameters
lambda= 1.062;
V_per = lambda * 10.0;
V_vri = lambda * 72.0;
V_pdp = lambda * 324.0;
V_clk = lambda * 1.0;
R_pbas = lambda * 0.02;
R_vbas = lambda * 0.18;
R_cbas = lambda * 0.001;
R_pdbas = lambda * 0.36;
K_pv = 0.2;
K_pp = 0.24;
K_ppd = 0.1;
K_vc = 0.54;
K_pdc = 0.54;
K_cv = 0.083;
K_cp = 0.134;
K_cpd = 0.248;
v_pcyt = lambda * 1.6;
K_pcyt = 0.25;
v_pnuc = lambda * 0.3;
K_pnuc = 0.001;
v_trans = lambda * 1.6;
K_trans = 0.25;
v_dclk = lambda * 0.2;
v_dvri = lambda * 0.7;
v_dpdp = lambda * 0.65;
F_v = lambda * 1.0;
F_p = lambda * 1.0;
F_pd = lambda * 1.0;
N = 5;
k_vdeac = lambda * 0.2;
k_pdeac = lambda * 0.2;
k_pddeac = lambda * 0.2;
tau_pdp_op = 3.0 / lambda;
tau_per_op = 3.0 / lambda;
tau_vri_op = 3.0 / lambda;
v_degp = lambda * 5.0;
K_degp = 0.01;
k_d = lambda * 0.005;
k_light = 0.685;
tau_delay = 3.0 / lambda;
%%%%%
AC_vri = 1e-3;
AC_per = 1e-3;
AC_pdp = 1e-3;
OP_vri = 1e-3;
OP_per = 1e-3;
OP_pdp = 1e-3;

% Forward Euler Parameters
T_total = 480;
dt = 2e-4;
tspan = 0:dt:T_total;

% Variable Storage & Initial Condition
l = length(tspan)-1;
CLK = [1.0, zeros(1,l)];
VRI = [1.0, zeros(1,l)];
PDP = [1.8, zeros(1,l)];
P0_cyt = [1.0/3, zeros(1,l)];
P1_cyt = [1.0/3, zeros(1,l)];
P2_cyt = [1.0/3, zeros(1,l)];
P0_nuc = [1.0/3, zeros(1,l)];
P1_nuc = [1.0/3, zeros(1,l)];
P2_nuc = [1.0/3, zeros(1,l)];
PER_cyt = P0_cyt + P1_cyt + P2_cyt;
PER_nuc = P0_nuc + P1_nuc + P2_nuc;
PER_tot = PER_cyt + PER_nuc;
R_pdp = zeros(1,l);

% Forward Euler Iterations
for i = 1:l
    
    % Transient change at t = 240 hr.
    if (i*dt >= 240)
        Q10_1 = 2;
        Q10_2 = 0.83;
        lambda_1= 1.062 * Q10_1 * 1.0269;
        lambda_2= 1.062 * Q10_2 * 1.0269;
        V_per = lambda_1 * 10.0;
        V_vri = lambda_1 * 72.0;
        V_pdp = lambda_1 * 324.0;
        V_clk = lambda_1 * 1.0;
        R_pbas = lambda_1 * 0.02;
        R_vbas = lambda_1 * 0.18;
        R_cbas = lambda_1 * 0.001;
        R_pdbas = lambda_1 * 0.36;
        K_pv = 0.2;
        K_pp = 0.24;
        K_ppd = 0.1;
        K_vc = 0.54;
        K_pdc = 0.54;
        K_cv = 0.083;
        K_cp = 0.134;
        K_cpd = 0.248;
%         v_pcyt = lambda_2 * 1.6;
%         v_pcyt = lambda_2 * 1.6 - (lambda_2 * 1.6 - lambda * 1.6)*exp(-0.5*(i*dt-240));
        v_pcyt = lambda_2 * 1.6 - (lambda_2 * 1.6 - lambda * 1.6)*exp(-0.05*(i*dt-240));
        K_pcyt = 0.25;
        v_pnuc = lambda_1 * 0.3;
        K_pnuc = 0.001;
        v_trans = lambda_1 * 1.6;
        K_trans = 0.25;
        v_dclk = lambda_1 * 0.2;
        v_dvri = lambda_1 * 0.7;
        v_dpdp = lambda_1 * 0.65;
        F_v = lambda_1 * 1.0;
        F_p = lambda_1 * 1.0;
        F_pd = lambda_1 * 1.0;
        N = 5;
        k_vdeac = lambda_1 * 0.2;
        k_pdeac = lambda_1 * 0.2;
        k_pddeac = lambda_1 * 0.2;
        tau_pdp_op = 3.0 / lambda_1;
        tau_per_op = 3.0 / lambda_1;
        tau_vri_op = 3.0 / lambda_1;
        v_degp = lambda_1 * 5.0;
        K_degp = 0.01;
        k_d = lambda_1 * 0.005;
        k_light = 0.685;
        tau_delay = 3.0 / lambda_1;
    end
    
    delay = round(tau_delay/dt);
    delay_index = max(i-delay,1);
    % Eqn 1~3
    k_vacet = F_v * CLK(i) / (CLK(i) + K_cv) * K_pv / (PER_nuc(i) + K_pv);
    k_pacet = F_p * CLK(i) / (CLK(i) + K_cp) * K_pp / (PER_nuc(i) + K_pp);
    k_pdacet = F_pd * CLK(i) / (CLK(i) + K_cpd) * K_ppd / (PER_nuc(i) + K_ppd);
    % Eqn 4~6, increments
    d_AC_vri = dt * ( k_vacet*(1-AC_vri) - k_vdeac*AC_vri );
    d_AC_per = dt * ( k_pacet*(1-AC_per) - k_pdeac*AC_per );
    d_AC_pdp = dt * ( k_pdacet*(1-AC_pdp) - k_pddeac*AC_pdp );
    % Eqn 7~9
    OP_vri_ss = AC_vri ^ N;
    OP_per_ss = AC_per ^ N;
    OP_pdp_ss = AC_pdp ^ N;
    % Eqn 10~12, increments
    d_OP_vri = dt * ( OP_vri_ss - OP_vri ) / tau_vri_op;
    d_OP_per = dt * ( OP_per_ss - OP_per ) / tau_per_op;
    d_OP_pdp = dt * ( OP_pdp_ss - OP_pdp ) / tau_pdp_op;
    % Eqn 13~16
    R_per = V_per * OP_per + R_pbas;
    R_vri = V_vri * OP_vri + R_vbas;
    R_pdp(i) = V_pdp * OP_pdp + R_pdbas;
    R_clk = V_clk * PDP(i)^2 / (PDP(i)^2 + K_pdc^2) * K_vc^2 / (VRI(i)^2 + K_vc^2) + R_cbas;
    % Eqn 17~19, increments
    d_CLK = dt * ( R_clk - v_dclk * CLK(i) - k_d * CLK(i) );
    d_VRI = dt * ( R_vri - v_dvri * VRI(i) - k_d * VRI(i) );
    d_PDP = dt * ( R_pdp(delay_index) - v_dpdp * PDP(i) - k_d * PDP(i) );
    % Eqn 20~22, increments
    d_P0_cyt = dt * ( R_per - v_pcyt * P0_cyt(i) / (K_pcyt + P0_cyt(i)) - k_d * P0_cyt(i) );
    d_P1_cyt = dt * ( v_pcyt * P0_cyt(i) / (K_pcyt + P0_cyt(i)) - v_pcyt * P1_cyt(i) / (K_pcyt + P1_cyt(i)) - k_d * P1_cyt(i) );
    d_P2_cyt = dt * ( v_pcyt * P1_cyt(i) / (K_pcyt + P1_cyt(i)) - v_trans * P2_cyt(i) / (K_trans + P2_cyt(i)) - k_d * P2_cyt(i) );
    % Eqn 23~25, increments
    d_P0_nuc = dt * ( v_trans * P2_cyt(i) / (K_trans + P2_cyt(i)) - v_pnuc * P0_nuc(i) / (K_pnuc + P0_nuc(i)) - k_d * P0_nuc(i) );
    d_P1_nuc = dt * ( v_pnuc * P0_nuc(i) / (K_pnuc + P0_nuc(i)) - v_pnuc * P1_nuc(i) / (K_pnuc + P1_nuc(i)) - k_d * P1_nuc(i) );
    d_P2_nuc = dt * ( v_pnuc * P1_nuc(i) / (K_pnuc + P1_nuc(i)) - v_degp * P2_nuc(i) / (K_degp + P2_nuc(i)) - k_d * P2_nuc(i) );
    
    % Update increments from Eqn 4~6, 10~12, 17~19, 20~22, 23~25, 26~28
    AC_vri = AC_vri + d_AC_vri;
    AC_per = AC_per + d_AC_per;
    AC_pdp = AC_pdp + d_AC_pdp;
    OP_vri = OP_vri + d_OP_vri;
    OP_per = OP_per + d_OP_per;
    OP_pdp = OP_pdp + d_OP_pdp;
    CLK(i+1) = CLK(i) + d_CLK;
    VRI(i+1) = VRI(i) + d_VRI;
    PDP(i+1) = PDP(i) + d_PDP;
    P0_cyt(i+1) = P0_cyt(i) + d_P0_cyt;
    P1_cyt(i+1) = P1_cyt(i) + d_P1_cyt;
    P2_cyt(i+1) = P2_cyt(i) + d_P2_cyt;
    P0_nuc(i+1) = P0_nuc(i) + d_P0_nuc;
    P1_nuc(i+1) = P1_nuc(i) + d_P1_nuc;
    P2_nuc(i+1) = P2_nuc(i) + d_P2_nuc;
    PER_nuc(i+1) = P0_nuc(i+1) + P1_nuc(i+1) + P2_nuc(i+1);
    PER_cyt(i+1) = P0_cyt(i+1) + P1_cyt(i+1) + P2_cyt(i+1);
    PER_tot(i+1) = PER_cyt(i+1) + PER_nuc(i+1);

end

% Plot Figures
figure(1)
hold on
h1 = plot(tspan, CLK, 'b-', 'LineWidth', 2);
h2 = plot(tspan, PER_tot, 'r-', 'LineWidth', 2);
h3 = plot(tspan, VRI, 'g-', 'LineWidth', 2);
h4 = plot(tspan, PDP, 'y-', 'LineWidth', 2);
xlabel('Time (hr)')
ylabel('Concentration (nM)')
legend('CLK', 'PER', 'VRI', 'PDP-1')

[sharedval,idx] = intersect(tspan,120:1:360,'stable');
outdata = table((tspan(idx)-240)',CLK(idx)',PER_tot(idx)',VRI(idx)',PDP(idx)', ...
    'VariableNames',{'time','CLK','PER_total','VRI','PDP'});

figure(2)
hold on
h5 = plot(outdata.time, outdata.CLK, 'b-', 'LineWidth', 2);
h6 = plot(outdata.time, outdata.PER_total, 'r-', 'LineWidth', 2);
h7 = plot(outdata.time, outdata.VRI, 'g-', 'LineWidth', 2);
h8 = plot(outdata.time, outdata.PDP, 'y-', 'LineWidth', 2);
xlabel('Time (hr)')
ylabel('Concentration (nM)')
legend('CLK', 'PER', 'VRI', 'PDP-1')

% writetable(outdata,'Data_Synthetic_nolambda.csv','Delimiter',',','QuoteStrings',true)
% writetable(outdata,'Data_Synthetic_lambda0.50.csv','Delimiter',',','QuoteStrings',true)
% writetable(outdata,'Data_Synthetic_lambda0.05.csv','Delimiter',',','QuoteStrings',true)


