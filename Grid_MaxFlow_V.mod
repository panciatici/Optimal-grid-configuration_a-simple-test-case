# Grid_MaxFlow_V.mod

param n >= 0,  integer;

set BUS := 1..n;

set BRANCH within {BUS,BUS};
set BREAKER within {BUS,BUS};

set N_X := {k in BUS, (k,m) in BRANCH} union {k in BUS, (m,k) in BRANCH};
set N_B := {k in BUS, (k,m) in BREAKER} union {k in BUS, (m,k) in BREAKER};

param Pgen {BUS};# Generation

param Pload {BUS}; # Load

param Zone {BUS}; #Zone Index

param BIGM = sum{k in BUS} Pload[k];

param Branch_X {BRANCH};	# branch reactance
param Branch_Fmax {BRANCH};	# branch maximun flow
param InterC {BRANCH}; 		# Interconnection index 

param B{(k,m) in N_X} := if ((k,m) in BRANCH) then 1000/Branch_X[k,m] else 1000/Branch_X[m,k];

# parameters for increasing the export

param TG1 := sum{k in BUS: Zone[k]==1} Pgen[k];
param TG2 := sum{k in BUS: Zone[k]==2} Pgen[k];

param TL1 := sum{k in BUS: Zone[k]==1} Pload[k];
param TL2 := sum{k in BUS: Zone[k]==2} Pload[k];

param alpha := TG1/TL2;
param beta := (TG2-TL1)/TL2;

var Breaker_F {N_B};

var Breaker_S {BREAKER} binary;

var Angle {BUS};

var Lambda;

maximize Export : Lambda;

subject to Balance_1 {k in BUS: Zone[k]==1}: 
Pload[k] + sum{(k,m) in N_X} B[k,m]*(Angle[k] - Angle[m]) + sum{(k,m) in N_B} Breaker_F[k,m] - Lambda*Pgen[k]=0;

subject to Balance_2 {k in BUS: Zone[k]==2}: 
(alpha*Lambda+beta)*Pload[k] + sum{(k,m) in N_X} B[k,m]*(Angle[k] - Angle[m]) + sum{(k,m) in N_B} Breaker_F[k,m] - Pgen[k]=0;

subject to Consistent_BF {(k,m) in BREAKER}:
   Breaker_F[k,m] = -Breaker_F[m,k];

subject to Status_breaker_p{(k,m) in BREAKER}:
   Breaker_F[k,m] <= BIGM *Breaker_S[k,m] ;
subject to Status_breaker_m{(k,m) in BREAKER}:
  -Breaker_F[k,m] <= BIGM * Breaker_S[k,m];

subject to Status_breaker_ap{(k,m) in BREAKER}:
  Angle[k]-Angle[m] <= BIGM *(1-Breaker_S[k,m]) ;
subject to Status_breaker_am{(k,m) in BREAKER}:
  Angle[m]-Angle[k] <= BIGM *(1-Breaker_S[k,m]) ;

subject to Lim_branch {(k,m) in BRANCH}:
  -Branch_Fmax[k,m] <=  B[k,m]*(Angle[k] - Angle[m]) <= Branch_Fmax[k,m] ;




