(* ::Package:: *)

BeginPackage["model`"]
Needs["contactmatrix`", FileNameJoin[{NotebookDirectory[], "contactmatrix.wl"}]]


modelGenerator::usage = "modelGenerator[beta, initial,alphaL,gammaa,gammas,gammah,ageDistribution,mu,ksi,asymp,pA,h]"
modelGenerator::usage = "modelGenerator[beta, initial,alphaL,gammaa,gammas,gammah,ageDistribution,mu,ksi,asymp,pA,h]"
modelSolver::usage = "modelSolver[beta, initial,alphaL,gammaa,gammas,gammah,ageDistribution,mu,ksi,asymp,pA,h]"
NGM::usage = "NGM[contact,r0, nAge, alphaL, gammaa,gammas,gammah,pA,asymp,ksi,h]"


latentPeriod::usage = "latentPeriod"
preSympInfPeriod::usage = "preSympInfPeriod"
infectPeriodA::usage = "gammaA"
infectPeriodS::usage = "gammaS"
infectPeriodH::usage = "gammaH"
infectPeriodC::usage = "gammaC"
infectPeriodCR::usage = "gammaC"
relInfectA::usage = "relInfectA"
mu::usage="mu"
pA::usage = "pA"
xi::usage = "xi"
h::usage = "h"


s;
varL={l1,l2}; varIp={ip}; nIp=1;
varIa={ia1,ia2,ia3};
varIs={is1,is2,is3};
varH={ih,ic,icr}; nH=Length[varH];
varOut={r,d,c}; nOut=Length[varOut];
nClassIa=3;
nClassIs=3;
nClassL=2;
t;


Begin["`Private`"]


(* ::Subsubsection::Closed:: *)
(*Parameters*)


assoc=Association[Import[StringJoin[ToString[NotebookDirectory[]],"/data/model_parameters.json"]]];
parameters=Association[Table[key->Association[assoc[key]]["value"],{key,Keys[assoc]}]];


(* Latent period *)
latentPeriod=1/parameters["alpha_l"]
(* pre-symptomatic infection period *)
preSympInfPeriod=1/parameters["alpha_p"];
(* asymptomatic infection period *)
infectPeriodA=1/parameters["gamma_a"];
(* symptomatic infection period *)
infectPeriodS= 1/parameters["gamma_s"];
(* hospitalization *)
infectPeriodH=1/parameters["gamma_h"];
(* intensive care until transition to Icr/D *)
infectPeriodC=1/parameters["gamma_c"];
(* intensive care recovery *)
infectPeriodCR=1/parameters["gamma_cr"];
(* relative infectiousness of Ia vs Is *)
relInfectA=parameters["inf_a"];
(* asymptomatic course *)
pA=parameters["p"];
(* probability of fatal outcome *)
mu=parameters["mu"];
(* hospitalization probability *)
h=parameters["h"];
(* probability of intensive care (given hospitalization) *)
xi=parameters["xi"];


(* ::Subsubsection::Closed:: *)
(*Variables*)


varGenerator[nAge_, patch_, nIa_, nIs_, nL_]:=Module[{vars, inits, var},
var=Join[{s}, Table[varL[[i]], {i,1,nL}], varIp, Table[varIa[[i]],{i,1,nIa}], Table[varIs[[i]],{i,1,nIs}], varH, varOut];
vars=Table[var/.Map[#->#[p,k][t]&,var],{p,1,patch},{k,1,nAge}];
inits=Table[var/.Map[#->#[p,k][0]&,var],{p,1,patch},{k,1,nAge}];
{vars,inits}
]


varGenerator[nAge_, patch_]:=Module[{nIa=nClassIa, nIs=nClassIs, nL=nClassL},
varGenerator[nAge,patch, nIa, nIs, nL]
]


varGenerator[nAge_, nIa_, nIs_, nL_]:=Module[{vars, inits,var},
var=Join[{s},Table[varL[[i]],{i,1,nL}],varIp,Table[varIa[[i]],{i,1,nIa}],Table[varIs[[i]],{i,1,nIs}],varH,varOut];
vars=Table[var/.Map[#->#[k][t]&,var],{k,1,nAge}];
inits=Table[var/.Map[#->#[k][0]&,var],{k,1,nAge}];
{{vars},{inits}}
]


(* ::Subsubsection:: *)
(*Model*)


modelGenerator[vars_, patch_, travel_, beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_, nIa_, nIs_, nL_]:=
Module[
{nAge, var, id, initv, initvalues, incidence, 
eqns, eqnl1, eqnl2, eqnip, eqnIaOthers, eqnIaFirst, eqnIa, eqnIsOthers, eqnIsFirst, eqnIs, eqnIh, eqnIc, eqnIcr, eqnR, eqnD, eqnC,
eqn, model},

(*This data structure allows us to refer to the index of a variable via its name*)
id=Association[Table[x[[1]]->x[[2]],{x,Transpose[{
(*compartments*)
Join[{s},varL[[Range[nL]]],varIp,varIa[[Range[nIa]]],varIs[[Range[nIs]]],varH,varOut], 
(*compartment indices*)
Range[1+nL+nIp+nIa+nIs+nH+nOut]
}]}]];

nAge=Length[ageDistribution[[1]]];
var=vars[[1]];
initv=vars[[2]];
initvalues=Flatten[Table[initv[[p,k,j]]==initial[[p,k,j]],{p,1,patch},{k,1,nAge},{j,1,Dimensions[initv][[3]]}]];

(*Incidence rate*)
incidence=Table[
Sum[
(*infection from pre-symptomatic class*)
beta[t][[p,j,k]]*var[[p,j,id[ip]]] 
(*infection from class with mild symptoms*)
+asymp*Sum[beta[t][[p,j,k]]*var[[p,j,id[o]]],{o,varIa[[Range[nIa]]]}] 
(*infection from symptomatic class*) 
+Sum[beta[t][[p,j,k]]*var[[p,j,id[o]]],{o,varIs[[Range[nIs]]]}],
{j,1,nAge}],
{p,1,patch},{k,1,nAge}];

(*S'[t]*)
eqns=Table[D[var[[p,k,id[s]]],t]==
(*infection*)
-var[[p,k,id[s]]]/ageDistribution[[p,k]]*incidence[[p,k]]
(*travel-out*)
-var[[p,k,id[s]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[s]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}];

(*L1'[t]*)
eqnl1=Table[D[var[[p,k,id[l1]]],t]==
(*infection*)
var[[p,k,id[s]]]/ageDistribution[[p,k]]incidence[[p,k]]
(*latency*)
-nL*alphaL[[k]]var[[p,k,id[l1]]]
(*travel-out*)
-var[[p,k,id[l1]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[l1]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}];

(*L2'[t]*)
eqnl2=Table[D[var[[p,k,id[l2]]],t]==
(*latency*)
nL*alphaL[[k]]var[[p,k,id[l1]]]
(*latency*)
-nL*alphaL[[k]]var[[p,k,id[l2]]]
(*travel-out*)
-var[[p,k,id[l2]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[l2]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}];

(*Ip'[t]*)
eqnip=Table[D[var[[p,k,id[ip]]],t]==
(*latency*)
nL*alphaL[[k]]var[[p,k,id[l2]]]
(*pre-symptomatic infection period*)
-alphaP[[k]]var[[p,k,id[ip]]]
(*travel-out*)
-var[[p,k,id[ip]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[ip]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}];

(*first Ia eqn]*)
eqnIaFirst={Table[D[var[[p,k,id[ia1]]],t]==
(*pre-symptomatic infection period*)
pA[[k]]*alphaP[[k]]var[[p,k,id[ip]]]
(*asymptomatic infection period*)
-nIa*gammaa[[k]]var[[p,k,id[ia1]]]
(*travel-out*)
-var[[p,k,id[ia1]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[ia1]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}]};

(*further Ia eqns*)
eqnIaOthers={};
If[nIa>1,
eqnIaOthers=Table[
Table[D[var[[p,k,id[ia1]+o]],t]==
(*asymptomatic infection period*)
nIa*gammaa[[k]]var[[p,k,id[ia1]+o-1]]
(*asymptomatic infection period contd*)
-nIa*gammaa[[k]]var[[p,k,id[ia1]+o]]
(*travel-out*)
-var[[p,k,id[ia1]+o]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[ia1]+o]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}],
{o,1,nIa-1}]];

(*All Ia eqns*)
eqnIa=Join[eqnIaFirst,eqnIaOthers];

(*first Is eqn*)
eqnIsFirst={Table[D[var[[p,k,id[is1]]],t]==
(*pre-symptomatic infection period*)
(1-pA[[k]])*alphaP[[k]]var[[p,k,id[ip]]]
(*symptomatic infection period*)
-nIs*gammas[[k]]var[[p,k,id[is1]]]
(*travel-out*)
-var[[p,k,id[is1]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[is1]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}]};

(*further Is eqn*)
eqnIsOthers={};
If[nIs>1,
eqnIsOthers=Table[
Table[D[var[[p,k,id[is1]+o]],t]==
(*symptomatic infection period*)
nIs*gammas[[k]]var[[p,k,id[is1]+o-1]]
(*symptomatic infection period contd*)
-nIs*gammas[[k]]var[[p,k,id[is1]+o]]
(*travel-out*)
-var[[p,k,id[is1]+o]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[is1]+o]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}],
{o,1,nIs-1}]];

(*All Is eqns*)
eqnIs=Join[eqnIsFirst,eqnIsOthers];

(*Ih'[t]*)
eqnIh=Table[D[var[[p,k,id[ih]]],t]==
(*symptomatic infection period contd*)
h[[k]] (1-ksi[[k]])*nIs*gammas[[k]]var[[p,k,id[is1]+nIs-1]]
(*hospitalization*)
-gammah[[k]]*var[[p,k,id[ih]]],
{p,1,patch},{k,1,nAge}];

(*Ic'[t])*)
eqnIc=Table[D[var[[p,k,id[ic]]],t]==
(*symptomatic infection period contd*)
h[[k]] ksi[[k]]*nIs*gammas[[k]]var[[p,k,id[is1]+nIs-1]]
(*intensive care until transition to D/Icr*)
-gammac[[k]]*var[[p,k,id[ic]]],{p,1,patch},{k,1,nAge}];

(*Icr'[t])*)
eqnIcr=Table[D[var[[p,k,id[icr]]],t]==
(*intensive care until transition to D/Icr*)
(1-mu[[k]])gammac[[k]]*var[[p,k,id[ic]]]
(*intensive care recovery*)
-gammacr[[k]]*var[[p,k,id[icr]]],{p,1,patch},{k,1,nAge}];

(*R'[t]*)
eqnR=Table[D[var[[p,k,id[r]]],t]==
(*asymptomatic infection period contd*)
nIa*gammaa[[k]]var[[p,k,id[ia1]+nIa-1]]
(*asymptomatic infection period contd*)
+(1-h[[k]])*nIs*gammas[[k]]*var[[p,k,id[is1]+nIs-1]]
(*hospitalization*)
+gammah[[k]]*var[[p,k,id[ih]]]
(*intensive care recovery*)
+gammacr[[k]]*var[[p,k,id[icr]]]
(*travel-out*)
-var[[p,k,id[r]]]*Sum[Which[(p1==p),0,True,travel[[p,p1,k]]],{p1,1,patch}]
(*travel-in*)
+Sum[Which[p2==p,0,True,var[[p2,k,id[r]]]*travel[[p2,p,k]]],{p2,1,patch}],
{p,1,patch},{k,1,nAge}];

(*D'[t]*)
eqnD=Table[D[var[[p,k,id[d]]],t]==
(*intensive care until transition to D/Icr*)
mu[[k]]* gammac[[k]]*var[[p,k,id[ic]]],{p,1,patch},{k,1,nAge}];

(*C'[t]*)
eqnC=Table[D[var[[p,k,id[c]]],t]==
(*pre-symptomatic infection period*)
var[[p,k,id[s]]]/ageDistribution[[p,k]]*incidence[[p,k]],{p,1,patch},{k,1,nAge}];

(*All eqns*)
eqn=Transpose[Join[{eqns, eqnl1, eqnl2, eqnip}, eqnIa, eqnIs, {eqnIh, eqnIc, eqnIcr, eqnR, eqnD, eqnC}]];
model=Flatten[Join[eqn,initvalues]];
{var, model}]


modelGenerator[vars_, beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_, nIa_, nIs_, nL_]:=Module[{patch=1, travel},
travel=Table[0,{i,1,patch},{j,1,patch},{k,1,Length[vars]}];
modelGenerator[vars, patch, travel, beta, initial, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h, nIa, nIs, nL]]


modelGenerator[beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_]:=Module[
{nIa=nClassIa,nIs=nClassIs,nL=nClassL,vars, nAge},
nAge=Length[ageDistribution];
vars=varGenerator[nAge,nIa,nIs, nL];
modelGenerator[vars, beta, {initial}, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, {ageDistribution}, mu, ksi, asymp, pA, h, nIa, nIs, nL]]


modelGenerator[patch_, travel_, beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_]:=Module[
{nIa=nClassIa, nIs=nClassIs, nL=nClassL, vars, nAge},
nAge=Length[ageDistribution[[1]]];
vars=varGenerator[nAge, patch, nIa, nIs, nL];
modelGenerator[vars, patch, travel, beta, initial, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h, nIa, nIs, nL]]


modelSolver[beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_]:=
Module[{T,modeloutput,vars,model,MO},
T=5000;
modeloutput = modelGenerator[beta, initial, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];
vars=modeloutput[[1]];
model=modeloutput[[2]];
MO=NDSolve[model,Flatten[vars],{t,0,T}]]


modelSolver[patch_, travel_, beta_, initial_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, ageDistribution_, mu_, ksi_, asymp_, pA_, h_]:=
Module[{T,modeloutput,vars,model,MO},
T=5000;
modeloutput = modelGenerator[patch, travel, beta, initial, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];
vars=modeloutput[[1]];
model=modeloutput[[2]];
MO=NDSolve[model,Flatten[vars],{t,0,T}]]


(* ::Subsubsection::Closed:: *)
(*NGM*)


NGM[patch_, travel_, contact_, r0_, nAge_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, pA_, asymp_, nIa_, nIs_, nL_]:=
Module[{id, nComp, V, F, NGMLargeDomain, EE, NextGen, betaBase, nIncubation, V2, travelMatrix, F2},
If[Length[contact[[1]]]>nAge, Message["Kontaktm\[AAcute]trix m\[EAcute]rete nagyobb, mint a nAgecsoportok sz\[AAcute]ma!"; Return[$Failed]]];

nIncubation=nIp+nL;
nComp=nIncubation+nIs+nIa;

id=Association[Table[x[[1]]->x[[2]],{x,Transpose[{
(*compartments*)
Join[varL,varIp,varIa[[Range[nIa]]],varIs[[Range[nIs]]]], 
(*indices of compartments*)
Range[nIncubation+nIa+nIs]
}]}]];

(*Transition matrix*)
V=Table[0,{i,1,nComp*nAge},{j,1,nComp*nAge}];
For[k=0,k<nAge,
{
(*Latency*)
V[[(k-1)*nComp+id[l1],(k-1)*nComp+id[l1]]]=nL*alphaL[[k]],
For[o=1,o<nL,o++,
(*Latency*)
{V[[(k-1)*nComp+(id[l1]+o),(k-1)*nComp+(id[l1]+o)-1]]=-nL*alphaL[[k]],
(*Latency contd*)
V[[(k-1)*nComp+(id[l1]+o),(k-1)*nComp+(id[l1]+o)]]=nL*alphaL[[k]]}
],

(*Latency*)
V[[(k-1)*nComp+id[ip],(k-1)*nComp+id[l1]+(nL-1)]]=-nL*alphaL[[k]],
(*Pre-symptomatic infection period*)
V[[(k-1)*nComp+id[ip],(k-1)*nComp+id[ip]]]=alphaP[[k]],

(*Pre-symptomatic infection period*)
V[[(k-1)*nComp+id[ia1],(k-1)*nComp+id[ip]]]=-alphaP[[k]]pA[[k]],
(*Asymptomatic infection period*)
V[[(k-1)*nComp+id[ia1],(k-1)*nComp+id[ia1]]]=nIa*gammaa[[k]],

If[nIa>1,
For[o=0,o<nIa-1,
(*Asymptomatic infection period*)
{V[[(k-1)*nComp+(id[ia1]+o),(k-1)*nComp+(id[ia1]+o)-1]]=-nIa gammaa[[k]],
(*Asymptomatic infection period contd*)
V[[(k-1)*nComp+(id[ia1]+o),(k-1)*nComp+(id[ia1]+o)]]=nIa gammaa[[k]]},
o++]
],

(*Pre-symptomatic infection period*)
V[[(k-1)*nComp+id[is1],(k-1)*nComp+id[ip]]]=-alphaP[[k]](1-pA[[k]]),
(*Symptomatic infection period*)
V[[(k-1)*nComp+id[is1],(k-1)*nComp+id[is1]]]=nIs*gammas[[k]],

If[nIs>1,
For[o=0,o<nIs-1,
(*Symptomatic infection period*)
{V[[(k-1)*nComp+(id[is1]+o),(k-1)*nComp+(id[is1]+o)-1]]=-nIs gammas[[k]],
(*Symptomatic infection period contd*)
V[[(k-1)*nComp+(id[is1]+o),(k-1)*nComp+(id[is1]+o)]]=nIs gammas[[k]]},
o++]
]},
k++];

(*NGM for multiple patches*)
If[patch>1,
(*block-diagonal matrix*)
V2=Table[Which[i==j,V,True,0*V],{i,1,patch},{j,1,patch}];
(*travel matrix*)
travelMatrix=Table[
Which[
(*travel matrix: travel-out appears in the diagonal*)
i==j,DiagonalMatrix[Flatten[
Table[
Sum[Which[(p==i),0,True,travel[[j,p,k]]],{p,1,patch}],
{k,1,nAge},{n,1,nComp}]]],
(*travel matrix: travel-in appears in the off-diagonal*)
True,DiagonalMatrix[Flatten[
Table[
-travel[[j,i,k]],
{k,1,nAge},{n,1,nComp}]]]],
{i,1,patch},{j,1,patch}];
(*final matrix (4D)*)
V2=V2+travelMatrix;
(*final matrix (2D)*)
V=Table[
V2[[ Quotient[i,nComp*nAge]+1, Quotient[j,nComp*nAge]+1, Mod[i,nComp*nAge]+1, Mod[j,nComp*nAge]+1 ]],
{i,0,nComp*nAge*patch-1},{j,0,nComp*nAge*patch-1}];
];

(*Transmission matrix*)
F=Table[0,{i,1,nComp*nAge},{j,1,nComp*nAge}];
For[k=0,k<nAge,
{For[j=0,j<nAge,
{
(*infection from pre-symptomatic class*)
F[[(k-1)*nComp+1,(j-1)*nComp+id[ip]]]=contact[[j,k]],
(*infection from asymptomatic class*)
For[o=0,o<nIa, F[[(k-1)*nComp+1,(j-1)*nComp+id[ia1]-1+o]]=asymp*contact[[j,k]], o++],
(*infection from symptomatic class*)
For[o=0,o<nIs, F[[(k-1)*nComp+1,(j-1)*nComp+id[is1]-1+o]]=contact[[j,k]], o++]
},
j++]},
k++];

(*NGM for multiple patches*)
If[patch>1,
(*block-diagonal: F matrix appears in the diagonal (4D)*)
F2=Table[Which[i==j,F,True,0*F],{i,1,patch},{j,1,patch}];
(*final matrix (2D)*)
F=Table[
F2[[ Quotient[i,nComp*nAge]+1, Quotient[j,nComp*nAge]+1, Mod[i,nComp*nAge]+1, Mod[j,nComp*nAge]+1 ]],
{i,0,nComp*nAge*patch-1},{j,0,nComp*nAge*patch-1}]
];

(*NGM calculation for large domain*)
NGMLargeDomain=F.Inverse[V];
(*transformation from large domain*)
EE=Table[0,{ps,1,nAge*patch},{pf,1,nComp*nAge*patch}];
For[ps=1,ps<(nAge*patch)+1,ps++,EE[[ps,nComp(ps-1)+1]]=1];
(*NGM*)
NextGen=EE.NGMLargeDomain.Transpose[EE];
(*calculation of beta based on R0 formula*)
betaBase=x/.FindRoot[x*Max[Abs[Eigenvalues[NextGen]]]==r0,{x,0.02}]
]


NGM[contact_, r0_, nAge_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_,  pA_, asymp_, nIa_, nIs_, nL_]:=Module[{patch=1,travel},
travel=Table[0,{i,1,patch},{j,1,patch},{k,1,1+nL+nIp+nIa+nIs+nH+nOut}];
NGM[patch, travel, contact, r0, nAge, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, pA, asymp, nIa, nIs, nL]]


NGM[contact_, r0_, nAge_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_,  pA_, asymp_]:=Module[{nIa=nClassIa, nIs=nClassIs, nL=nClassL},
NGM[contact, r0, nAge, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, pA, asymp, nIa, nIs, nL]
]


NGM[patch_, travel_, contact_, r0_, nAge_, alphaL_, alphaP_, gammaa_, gammas_, gammah_, gammac_, gammacr_, pA_, asymp_]:=Module[{nIa=nClassIa, nIs=nClassIs, nL=nClassL},
NGM[patch, travel, contact, r0, nAge, alphaL, alphaP, gammaa, gammas, gammah, gammac, gammacr, pA, asymp, nIa, nIs, nL]
]


(* ::Subsubsection::Closed:: *)
(*End*)


End[ ]


EndPackage[ ]
