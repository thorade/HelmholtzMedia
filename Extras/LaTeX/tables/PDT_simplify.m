clear;
syms dP_dD_T dP_dT_D
syms d2P_dD2_T d2P_dT2_D d2P_dTdD
syms d2T_dD2_P

d2T_dD2_P = -(d2P_dD2_T*dP_dT_D - dP_dD_T*d2P_dTdD)/dP_dT_D^2 ...
            +(d2P_dTdD*dP_dT_D - dP_dD_T*d2P_dT2_D)*dP_dD_T/dP_dT_D^3
simplify(d2T_dD2_P)
expand(d2T_dD2_P)

d2D_dT2_P = -(d2P_dT2_D*dP_dD_T - dP_dT_D*d2P_dTdD)/dP_dD_T^2 ...
            +(d2P_dTdD*dP_dD_T - dP_dT_D*d2P_dD2_T)*dP_dT_D/dP_dD_T^3
simplify(d2D_dT2_P)
expand(d2D_dT2_P)
