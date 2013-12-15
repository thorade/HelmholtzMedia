clear;
syms D T;
syms dP_D dP_T dU_T;
syms d2P_D2 d2P_T2 d2P_TD d2U_T2;

dS_D = -1/D^2*dP_T
dS_T = 1/T*dU_T

d2S_D2 = -1/D^2*d2P_TD + 2/D^3*dP_T
d2S_T2 = 1/T*d2U_T2 - 1/T^2*dU_T
d2S_TD = -1/D^2*d2P_T2

d2P_D2_S = d2P_D2 -(dP_T*d2S_D2   + 2*d2P_TD*dS_D)/dS_T ...
                  +(d2P_T2*dS_D^2 + 2*dP_T*dS_D*d2S_TD)/dS_T^2 ...
                  -(dP_T  *dS_D^2*d2S_T2)/dS_T^3

simplify(d2P_D2_S)
expand(d2P_D2_S)
