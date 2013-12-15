clear;

syms A B C D 
syms dA_dd_s dB_dd_s dC_dd_s dD_dd_s
syms AD AT BD BT CD CT DD DT

dA_dd_s = AD - AT*C/D
dB_dd_s = AT - BT*C/D
dC_dd_s = CD - CT*C/D
dD_dd_s = CT - DT*C/D

d2p_dd2_s = dA_dd_s - C/D*dB_dd_s -B/D*dC_dd_s + B*C/D/D*dD_dd_s

simplify(d2p_dd2_s)