function sigma_yield = estimated_yield_strength( T_iso, wC, wSi, wCr, wMo, wN  )
% rough estimation of yield strength based on some element contents and 
% isothermal transformation temperature in [°C] due to: 
% Singh, Bhadeshia MSE-A 1998 - Estimation of bainite 
% plate-thickness in low-alloy steels

T_iso = T_iso - 25.;
    
sigma_yield = (1. -0.26e-2 *T_iso +0.47e-5* T_iso^2 - ...
              0.326e-8 * T_iso^3) * 15.4 * (4.4 + 23* wC + ...
              1.3*wSi +0.24*wCr + 0.94*wMo + 0.32*wN);


end

