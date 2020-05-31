#!/usr/bin/env wolframscript -cloud
(* ::Package:: *)

Import[StringJoin[NotebookDirectory[],"/modules/contactmatrix.wl"]]
Import[StringJoin[NotebookDirectory[],"/modules/model.wl"]]


generated=StringJoin[NotebookDirectory[],"/generated/"];


(* ::Section::Closed:: *)
(*Preliminary calculations*)


ifr={0.002, 0.006, 0.03, 0.08, 0.15, 0.6, 2.2, 5.1, 9.3}; (*Ferguson*)
hh={0.03, 0.03, 0.39, 1.45, 2.55, 4.95, 7.75, 17.88, 32.97}; (*McCaw*)
icu={0.01, 0.01, 0.11, 0.43, 0.75, 1.45, 2.27, 5.25, 9.68}; (*McCaw*)
hPs={0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3}; (*Ferguson*)
iPh={5.0, 5.0, 5.0, 5.0, 6.3, 12.2, 27.4, 43.2, 70.9}; (*Ferguson*)
asympPinf={0.95, 0.95,0.90,0.8,0.7,0.6,0.4,0.2,0.2}; (**)
icuPhosp=icu/hh; (*McCaw - derived*)
fPicu=(ifr/100)/((hPs/100)*(iPh/100)) (*Ferguson - derived*)


binAll={Sum[population1Y[[i]],{i,1,5}], (*0-4*)
Sum[population1Y[[i]],{i,6,10}], (*5-9*)
Sum[population1Y[[i]],{i,11,15}], (*10-14*)
Sum[population1Y[[i]],{i,16,20}], (*15-19*)
Sum[population1Y[[i]],{i,21,30}], (*20-29*)
Sum[population1Y[[i]],{i,31,40}], (*30-39*)
Sum[population1Y[[i]],{i,41,50}], (*40-49*)
Sum[population1Y[[i]],{i,51,60}], (*50-59*)
Sum[population1Y[[i]],{i,61,70}], (*60-69*)
Sum[population1Y[[i]],{i,71,80}], (*70-79*)
Sum[population1Y[[i]],{i,81,Length[population1Y]}]} (*80+ *)
sumAll=Sum[binAll[[i]],{i,1,Length[binAll]}];


(*ICU | Hospitalization*)
(*Based on McCaw parametrization*)
{
icuPhosp[[1]], (*0-4*)

icuPhosp[[1]] * binAll[[2]]/(binAll[[2]]+binAll[[3]]) + icuPhosp[[2]] * binAll[[3]]/(binAll[[2]]+binAll[[3]]), (*5-14*)

icuPhosp[[2]] * binAll[[4]]/(binAll[[4]]+binAll[[5]]) + icuPhosp[[3]] * binAll[[5]]/(binAll[[4]]+binAll[[5]]),(*15-29*)

icuPhosp[[4]] * binAll[[6]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + icuPhosp[[5]] * binAll[[7]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + icuPhosp[[6]] * binAll[[8]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]), (*30-59*)

icuPhosp[[7]], (*60-69*)

icuPhosp[[8]], (*70-79*)

icuPhosp[[9]] (*80+*)
}


(*Hospitalization | Symptoms*)
(*Based on McCaw parametrization*)
{
hh[[1]], (*0-4*)

hh[[1]] * binAll[[2]]/(binAll[[2]]+binAll[[3]]) + hh[[2]] * binAll[[3]]/(binAll[[2]]+binAll[[3]]), (*5-14*)

hh[[2]] * binAll[[4]]/(binAll[[4]]+binAll[[5]]) + hh[[3]] * binAll[[5]]/(binAll[[4]]+binAll[[5]]),(*15-29*)

hh[[4]] * binAll[[6]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + hh[[5]] * binAll[[7]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + hh[[6]] * binAll[[8]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]), (*30-59*)

hh[[7]], (*60-69*)

hh[[8]], (*70-79*)

hh[[9]] (*80+*)
}/100


(*Fatal | ICU - temporal, scaled*)
(*Based on Ferguson parametrization*)
fPiScaled={
fPicu[[1]], (*0-4*)

fPicu[[1]] * binAll[[2]]/(binAll[[2]]+binAll[[3]]) + fPicu[[2]] * binAll[[3]]/(binAll[[2]]+binAll[[3]]), (*5-14*)

fPicu[[2]] * binAll[[4]]/(binAll[[4]]+binAll[[5]]) + fPicu[[3]] * binAll[[5]]/(binAll[[4]]+binAll[[5]]),(*15-29*)

fPicu[[4]] * binAll[[6]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + fPicu[[5]] * binAll[[7]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + fPicu[[6]] * binAll[[8]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]), (*30-59*)

fPicu[[7]], (*60-69*)

fPicu[[8]], (*70-79*)

fPicu[[9]] (*80+*)
}


wsum=fPiScaled[[1]]*binAll[[1]]/sumAll+
fPiScaled[[2]]*(binAll[[2]]+binAll[[3]])/sumAll+
fPiScaled[[3]]*(binAll[[4]]+binAll[[5]])/sumAll+
fPiScaled[[4]]*(binAll[[6]]+binAll[[7]]+binAll[[8]])/sumAll+
fPiScaled[[5]]*(binAll[[9]])/sumAll+
fPiScaled[[6]]*binAll[[10]]/sumAll+
fPiScaled[[7]]*binAll[[11]]/sumAll


fPiScaled*0.5/wsum


(*Mild symptoms | Infection*)
(*Based on  parametrization*)
mildScaled={
asympPinf[[1]], (*0-4*)

asympPinf[[1]] * binAll[[2]]/(binAll[[2]]+binAll[[3]]) + asympPinf[[2]] * binAll[[3]]/(binAll[[2]]+binAll[[3]]), (*5-14*)

asympPinf[[2]] * binAll[[4]]/(binAll[[4]]+binAll[[5]]) + asympPinf[[3]] * binAll[[5]]/(binAll[[4]]+binAll[[5]]),(*15-29*)

asympPinf[[4]] * binAll[[6]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + asympPinf[[5]] * binAll[[7]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]) + asympPinf[[6]] * binAll[[8]]/(binAll[[6]]+binAll[[7]]+binAll[[8]]), (*30-59*)

asympPinf[[7]], (*60-69*)

asympPinf[[8]], (*70-79*)

asympPinf[[9]] (*80+*)
}


wsum2=
asympPinf[[1]]*binAll[[1]]/sumAll+
asympPinf[[2]]*(binAll[[2]]+binAll[[3]])/sumAll+
asympPinf[[3]]*(binAll[[4]]+binAll[[5]])/sumAll+
asympPinf[[4]]*(binAll[[6]]+binAll[[7]]+binAll[[8]])/sumAll+
asympPinf[[5]]*(binAll[[9]])/sumAll+
asympPinf[[6]]*binAll[[10]]/sumAll+
asympPinf[[7]]*binAll[[11]]/sumAll


mildScaled*0.5/wsum2


(* ::Section::Closed:: *)
(*Our model*)


(* ::Subsection::Closed:: *)
(*Initialization*)


(* ::Subsubsection::Closed:: *)
(*Parameters*)


nAge=Length[ageDistribution];
contactMatrix=aggregatePREM[contactPREM];


alpha=Table[1/latentPeriod,{i,1,7}];
alpha2=Table[1/preSympInfPeriod,{i,1,7}];
gammaa=Table[1/infectPeriodA,{i,1,7}];
gammas=Table[1/infectPeriodS,{i,1,7}];


gammah=Table[1/infectPeriodH,{i,1,7}];
gammac=Table[1/infectPeriodC,{i,1,7}];
gammacr=Table[1/infectPeriodCR,{i,1,7}];


ksi=xi;
asymp=relInfectA;


(* ::Subsubsection::Closed:: *)
(*Initial values: \[CapitalSigma] I = 5000*)


R0=1.35;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


l10={0,0,1,1,1,0,0};
l20={0,0,0,0,0,0,0};
ip0={0,0,0,0,0,0,0};
ia10={0,0,0,0,0,0,0};
ia20={0,0,0,0,0,0,0};
ia30={0,0,0,0,0,0,0};
is10={0,0,0,0,0,0,0};
is20={0,0,0,0,0,0,0};
is30={0,0,0,0,0,0,0};
ih0={0,0,0,0,0,0,0};
ic0={0,0,0,0,0,0,0};
icr0={0,0,0,0,0,0,0};
r0={0,0,0,0,0,0,0};
d0={0,0,0,0,0,0,0};
c0={0,0,0,0,0,0,0};
s0=ageDistribution-l10;
init=Transpose[{s0,l10,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];


sol=modelSolver[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


sumOfI[t_]:=Evaluate[Sum[ip[i][t]+ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.sol];


Plot[sumOfI[t],{t,0,200}]


t0=t/.FindRoot[sumOfI[t]==5000,{t,150}][[1]]


initialValuesI=Flatten[modelGenerator[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h][[1]]/.sol/.t->t0,2]


initialValuesI[[7,-5]]+=10;
initialValuesI[[6,-5]]+=7;


initialValuesI=Round[initialValuesI]


(* ::Subsubsection::Closed:: *)
(*Functions*)


mortalityCumulative=Table[Sum[d[j][t],{j,1,i}],{i,1,nAge}];
cumulative=Table[Sum[c[j][t],{j,1,i}],{i,1,nAge}];


(* ::Section::Closed:: *)
(*Scenarios post-lockdown*)


ticklabel1={{0,"1.May"},{31,"1.June"},{31+30,"1.July"},{31+30+31,"1.Aug"},{31+30+31+31,"1.Sept"},{31+30+31+31+30,"1.Oct"},{31+30+31+31+30+31,"1.Nov"},{31+30+31+31+30+31+30,"1.Dec"},{31+30+31+31+30+31+30+31,"1.Jan"}};
ticklabel2={{0,"1.May"},{31+30,"1.July"},{31+30+31+31,"1.Sept"},{31+30+31+31+30+31,"1.Nov"},{31+30+31+31+30+31+30+31,"1.Jan"},{31+30+31+31+30+31+30+31+31+28,"1.Feb"}};
ticklabel3={{0,"1.May"},{31+30+31+31,"1.Sept"},{31+30+31+31+30+31+30+31,"1.Jan"},{365,"1.May"},{365+31+30+30+31,"1.Sept"},{365+31+30+30+31+30+31+30+31,"1.Jan"}};


twoPlotsStyle={
{Thickness[0.01],ColorData["Crayola"]["MidnightBlue"]},
{Thickness[0.01],ColorData["Crayola"]["Blue"]}};


ageStructureLegend=LineLegend[{"0-4","5-14","15-29","30-59","60-69","70-79","80+"}, LegendLayout->"ReversedColumn", LegendMarkers->Graphics[{Thickness[0.4], Line[{{-10, 0}, {10, 0}}]}]];


(* ::Subsection::Closed:: *)
(*Unmitigated scenario*)


R0=2.2;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


(*ICU needs*)
plt=Plot[
{Sum[ih[i][t],{i,1,nAge}]/.sol, 
Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,150},
PlotRange->{0,All},
LabelStyle->Directive[18],
PlotLabel->Style["ICU Needs"],
PlotStyle->twoPlotsStyle,
PlotRange->All,
Ticks->{ticklabel1,Automatic}]


Export[StringJoin[generated,"unmitigated_hosp.pdf"],plt]


(*Epidemic curve*)
plt=Plot[{
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[is1[i][t]+is2[i][t]+is3[i][t],{i,1,nAge}]/.sol]},{t,0,150},
PlotRange->{0,All},
LabelStyle->Directive[18],
PlotLabel->Style["Epidemic curve"],
PlotStyle->twoPlotsStyle,
PlotRange->All,
Ticks->{ticklabel1,Automatic}]


Export[StringJoin[generated,"unmitigated_epi.pdf"],plt]


(*Mortality*)
plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,150},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{ticklabel1,Automatic}]


Export[StringJoin[generated,"unmitigated_mort.pdf"],plt]


(*Cumulative*)
plt=Plot[Evaluate[cumulative/.sol],{t,0,150},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, Black, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{ticklabel1,
{{1000000,"1 Million"},{2000000,"2 Million"},{3000000,"3 Million"},{4000000,"4 Million"},{5000000,"5 Million"},{6000000,"6 Million"},{7000000,"7 Million"}}}]


Export[StringJoin[generated,"unmitigated_cumu.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Weak control (25% transmission reduction)*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Sum[ih[i][t],{i,1,nAge}]/.sol, 
Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,200},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
AxesStyle->Directive[Thickness[0.005]],
PlotLabel->Style["ICU Needs"],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel1,Automatic}, ImageSize->650]


Export[StringJoin[generated,"weak_hosp.pdf"],plt]


(*Epidemic curve*)
plt=Plot[
{
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[is1[i][t]+is2[i][t]+is3[i][t],{i,1,nAge}]/.sol]},{t,0,200},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
PlotLabel->Style["Epidemic curve"],
AxesStyle->Directive[Thickness[0.005]],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel1,Automatic}, ImageSize->650]


Export[StringJoin[generated,"weak_epi.pdf"],plt]


(*Mortality*)
plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,200},
PlotLabel->Style["Mortality"],
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black,30],
AxesStyle->Directive[Thickness[0.0075]],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},Background->White,
Ticks->{ticklabel1,Automatic}, ImageSize->800]


Export[StringJoin[generated,"weak_mort.pdf"],plt]


(*Cumulative*)
plt=Plot[Evaluate[cumulative/.sol],{t,0,200},
PlotRange->{0,All},
PlotLabel->Style["Recovery"],
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black,30],
AxesStyle->Directive[Thickness[0.005]],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},Background->White,
Ticks->{ticklabel1,
{{1000000,"1 Million"},{2000000,"2 Million"},{3000000,"3 Million"},{4000000,"4 Million"},{5000000,"5 Million"},{6000000,"6 Million"},{7000000,"7 Million"}}}, ImageSize->800]


Export[StringJoin[generated,"weak_cumu.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Moderate control (40% transmission reduction)*)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Sum[ih[i][t],{i,1,nAge}]/.sol, 
Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,250},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
AxesStyle->Directive[Thickness[0.005]],
PlotLabel->Style["ICU Needs"],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel1,Automatic}, ImageSize->800]


Export[StringJoin[generated,"moderate_hosp.pdf"],plt]


(*Epidemic curve*)
plt=Plot[
{
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[is1[i][t]+is2[i][t]+is3[i][t],{i,1,nAge}]/.sol]},{t,0,250},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
AxesStyle->Directive[Thickness[0.005]],
PlotLabel->Style["Epidemic curve"],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel1,Automatic}, ImageSize->800]


Export[StringJoin[generated,"moderate_epi.pdf"],plt]


(*Mortality*)
plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,350},
PlotRange->{0,All},
PlotLabel->Style["Mortality"],
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black,30],
AxesStyle->Directive[Thickness[0.005]],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},Background->White,
Ticks->{ticklabel2,Automatic},ImageSize->800]


Export[StringJoin[generated,"moderate_mort.pdf"],plt]


(*Cumulative*)
plt=Plot[Evaluate[cumulative/.sol],{t,0,350},
PlotRange->{0,All},
PlotLabel->Style["Recovery"],
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black, 30],
AxesStyle->Directive[Thickness[0.005]],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{ticklabel2,
{{1000000,"1 Million"},{2000000,"2 Million"},{3000000,"3 Million"},{4000000,"4 Million"},{5000000,"5 Million"},{6000000,"6 Million"},{7000000,"7 Million"}}},Background->White,ImageSize->800]


Export[StringJoin[generated,"moderate_cumu.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Strong control (50% transmission reduction)*)


R0=1.1;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Sum[ih[i][t],{i,1,nAge}]/.sol, 
Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
PlotLabel->Style["ICU Needs"],
AxesStyle->Directive[Thickness[0.005]],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel2,Automatic}, ImageSize->800]


Export[StringJoin[generated,"strong_hosp.pdf"],plt]


(*Epidemic curve*)
plt=Plot[
{
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[is1[i][t]+is2[i][t]+is3[i][t],{i,1,nAge}]/.sol]},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black,25],
PlotLabel->Style["Epidemic curve"],
AxesStyle->Directive[Thickness[0.005]],
PlotStyle->twoPlotsStyle,
PlotRange->All,Background->White,
Ticks->{ticklabel2,Automatic}, ImageSize->800]


Export[StringJoin[generated,"strong_epi.pdf"],plt]


(*Mortality*)
plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,630},
PlotRange->{0,All},
PlotLabel->Style["Mortality"],
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black,30],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},Background->White,
Ticks->{ticklabel3,Automatic},ImageSize->800]


Export[StringJoin[generated,"strong_mort.pdf"],plt]


(*Cumulative*)
plt=Plot[Evaluate[cumulative/.sol],{t,0,630},
PlotRange->{0,All},
PlotLabel->Style["Recovery"],
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->Placed[ageStructureLegend,{Right, Right}], AxesStyle->Directive[Thickness[0.005],Black,FontSize->25], LabelStyle->Directive[Bold, Black,30],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{ticklabel3,{500000,{1000000,"1 Million"},{1500000,"1.5 Million"},{2000000,"2 Million"},{2500000,"2.5 Million"}}},Background->White,ImageSize->800]


Export[StringJoin[generated,"strong_cumu.pdf"],plt]


(* ::Section::Closed:: *)
(*Age-specific analysis*)


(* ::Subsection::Closed:: *)
(*School closures*)


(* ::Subsubsection::Closed:: *)
(*Contact matrix*)


contactMatrixPREMSchool=Table[
Which[(4>=i>=2 || 4>=j>=2),contactPREM[[i,j]]-contactPREMschool[[i,j]]-0.5*contactPREMother[[i,j]],
True,contactPREM[[i,j]]-contactPREMschool[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixSchool=aggregatePREM[contactMatrixPREMSchool];


(* ::Subsubsection::Closed:: *)
(*Model*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrixSchool};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


R0*baseBeta/NGM[contactMatrixSchool, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]


(* ::Subsubsection::Closed:: *)
(*Plots*)


plt=Plot[
{Sum[ih[i][t],{i,1,nAge}]/.sol, 
Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, 20],
PlotLabel->Style["ICU Needs"],
PlotStyle->twoPlotsStyle,
PlotRange->All,
Ticks->{{{0,"1.May"},{31+30+31+31,"1.September"},{31+30+31+30+31+30+31+31,"1.February"}},Automatic}]


Export[StringJoin[generated,"school_hosp.pdf"],plt]


(*Epidemic curve*)
plt=Plot[
{
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[is1[i][t]+is2[i][t]+is3[i][t],{i,1,nAge}]/.sol]},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, 20],
PlotLabel->Style["Epidemic curve"],
PlotStyle->twoPlotsStyle,
PlotRange->All,
Ticks->{{{0,"1.May"},{31+30+31+31,"1.September"},{31+30+31+30+31+30+31+31,"1.February"}},Automatic}]


Export[StringJoin[generated,"school_epi.pdf"],plt]


(*Mortality*)
plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,300},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{{{0,"1.May"},{31+30+31+31,"1.September"},{31+30+31+30+31+30+31+31,"1.February"}},Automatic}]


Export[StringJoin[generated,"school_mort.pdf"],plt]


(*Cumulative*)
plt=Plot[Evaluate[cumulative/.sol],{t,0,300},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{{{0,"1.May"},{31+30+31+31,"1.September"},{31+30+31+30+31+30+31+31,"1.February"}},
{{1000000,"1 Million"},{2000000,"2 Million"},{3000000,"3 Million"},{4000000,"4 Million"},{5000000,"5 Million"},{6000000,"6 Million"},{7000000,"7 Million"}}}]


Export[StringJoin[generated,"school_cumu.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Elderly protection*)


(* ::Subsubsection::Closed:: *)
(*Contact matrix*)


contactMatrixPREMElderlyFull=Table[
Which[(i>=13 || j>=13),contactPREM[[i,j]]-contactPREMother[[i,j]]-contactPREMschool[[i,j]]-contactPREMwork[[i,j]],
True,contactPREM[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixElderlyFull=aggregatePREM[contactMatrixPREMElderlyFull];


contactMatrixPREMElderlyModerate=Table[
Which[(i>=13 || j>=13),contactPREM[[i,j]]-contactPREMschool[[i,j]]-0.5*contactPREMother[[i,j]]-0.5*contactPREMwork[[i,j]],
True,contactPREM[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixElderlyModerate=aggregatePREM[contactMatrixPREMElderlyModerate];


(* ::Subsubsection::Closed:: *)
(*Model*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*contactMatrixElderlyModerate};
solElderlyModerate=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*contactMatrixElderlyFull};
solElderlyFull=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


(* ::Subsubsection::Closed:: *)
(*Plots: ICU*)


plt=Plot[
{Evaluate[Sum[ic[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyModerate],
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyFull]
},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, 20],
PlotLabel->Style["ICU needs"],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["OuterSpace"]},
{Thickness[0.01],ColorData["Crayola"]["Silver"]},
{Thickness[0.01],ColorData["Crayola"]["Periwinkle"]}},
PlotRange->All,
Ticks->{
{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"}},
All}]


Export[StringJoin[generated,"elderly_hosp.pdf"],plt]


(* ::Subsubsection::Closed:: *)
(*Plots: mortality*)


plt=Plot[Evaluate[mortalityCumulative/.sol],{t,0,250},
PlotRange->{0,20000},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLabel->"Mortality",
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"}}, 
Automatic}]


Export[StringJoin[generated,"elderly_hosp_1.pdf"],plt]


plt=Plot[Evaluate[mortalityCumulative/.solElderlyModerate],{t,0,250},
PlotRange->{0,20000},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLabel->"Mortality",
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"}}, 
Automatic}]


Export[StringJoin[generated,"elderly_hosp_2.pdf"],plt]


plt=Plot[Evaluate[mortalityCumulative/.solElderlyFull],{t,0,250},
PlotRange->{0,20000},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLabel->"Mortality",
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
Ticks->{{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"}}, 
Automatic}]


Export[StringJoin[generated,"elderly_hosp_3.pdf"],plt]


(* ::Section::Closed:: *)
(*Sensitivity analysis*)


runThisSection=False;


R0=1.65;
stepSize=40;


hmultInit=0.7;
hmultEnd=2.3;
R0init=1.0;
R0end=2.5;
maximumList=Table[{0,0,0},{hh,1,1+stepSize*(hmultEnd-hmultInit)},{rr0,1,1+stepSize*(R0end-R0init)}];


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];
attack=(Flatten[(Table[c[i][t],{i,1,nAge}]/.sol)/.t->300]/ageDistribution);
If[runThisSection,
	rr0=1;
	For[Rt=R0init,Rt<=R0end,Rt+=1/stepSize,
		hh=1;
		baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
		beta[t_]={baseBeta*(Rt/R0)*contactMatrix};
		For[hmult=hmultInit,hmult<=hmultEnd,hmult+=1/stepSize,
			sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h*hmult];
			target[t_]:=Sum[ic[i][t],{i,1,nAge}]/.sol;
			ifrr=Total[ageDistribution*attack*(1-pA)*(h*hmult)*xi*mu]/Total[ageDistribution*attack];
			maximumList[[hh,rr0]]={ifrr, Rt, Quiet[FindMaximum[target[t],{t,80}]][[1]]};
			hh++;
		];
		rr0++;
	]
]//Timing


IFRdistinct=0; maximumDistinct=0;
distinctVector=maximumList[[1+Round[(1.0-hmultInit)*stepSize],1+Round[(1.65-R0init)*stepSize]]]


legendMax=40000;
plot=Show[ListContourPlot[Flatten[maximumList,1],
PlotRange->{{0.00289,0.0089},{1.1,2.4}},
Contours->({#, If[Round[#] == 8000,
				Directive[Thickness[0.01], White],
				Null
		  ]}& /@ Range[ 0, legendMax, 1000]),
PlotLegends->BarLegend[{Automatic, {0, legendMax}}, Table[i,{i,0,40000,1000}]], 
ColorFunction->ColorData[{"RedBlueTones", {legendMax,0}}],
ColorFunctionScaling -> False,
FrameLabel->{"IFR",Style["\[ScriptCapitalR]",FontSize->15]},LabelStyle->Directive[Bold, 12, Black]],
Graphics[{White,PointSize[0.025],Point[{0.00393415,1.695}]}]]


Export[StringJoin[generated,"contourPlot.pdf"],plot]


(* ::Section::Closed:: *)
(*Seasonality*)


(* ::Subsection::Closed:: *)
(*Seasonality*)


maxChange=0.3;
period[change_,t_]:=0.5*change*Cos[(t*2*Pi)/366]+(1-0.5*change);
inset=Plot[{1,period[0.1,t],period[0.2,t],period[0.3,t]},{t,0,365},PlotRange->{1-1.5*maxChange,1.05},LabelStyle->Directive[Bold, Black, 20],
PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
AxesStyle-> Directive[Thickness[0.005]],
Ticks->{
{{0,"1. December"},{31+31+28+31+30+31,"1. June"},{364,"30. November"}},
{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}}, ImageSize->420]


Export[StringJoin[generated,"seas_inset.pdf"],inset]


(* ::Subsection::Closed:: *)
(*Analysis*)


(* ::Subsubsection::Closed:: *)
(*Weak control*)


R0=1.65;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]
},{t,0,300},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black, 25],
AxesLabel->{"time","number of infected individuals"},
AxesStyle->Directive[Thickness[0.005]],
PlotLabel->Style["R=1.65"],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotRange->All,
Ticks->{
{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"}},
All},
Epilog->Inset[inset,{240,480000}],ImageSize->1200, Background->White]


Export[StringJoin[generated,"seas_weak.pdf"],plt]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP03],{t,150}]


N[321565/293400]


(* ::Subsubsection::Closed:: *)
(*Moderate control*)


R0=1.32;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]
},{t,0,500},
PlotRange->{{0,500},{0,270000}},
LabelStyle->Directive[Bold, Black, 25],
PlotLabel->Style["R=1.32"],
AxesStyle->Directive[Thickness[0.005]],
AxesLabel->{"time","number of infected individuals"},
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotRange->All,
Ticks->{
{{0,"1. May"},{31+30+31,"1. August"},{31+30+31+31+30+31+30,"1. December"},{365,"1. May"}},
All},
Epilog->Inset[inset,{390,200000}],ImageSize->1200,Background->White]


Export[StringJoin[generated,"seas_moderate.pdf"],plt]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP01],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP03],{t,150}]


N[148787/128675]


N[160612/128675]


N[160612/148787]


(* ::Subsubsection::Closed:: *)
(*Strong control*)


R0=1.1;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesI, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]
},{t,0,700},
PlotRange->{0,All},
LabelStyle->Directive[Bold, Black, 25],
AxesStyle->Directive[Thickness[0.005]],
PlotLabel->Style["R=1.1"],
AxesLabel->{"time","number of infected individuals"},
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotRange->All,
Ticks->{
{{0,"1. May"},{31+30+31+31+30+31+30,"1. December"},{365,"1. May"},{365+31+30+31+31+30+31+30,"1. December"}},
Automatic},
Epilog->Inset[inset,{520,75000}],ImageSize->1200, Background->White]


Export[StringJoin[generated,"seas_strong.pdf"],plt]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP01],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP03],{t,150}]


N[17807/71697]


(* ::Section::Closed:: *)
(*Spatial heterogeneity*)


(* ::Subsection::Closed:: *)
(*Bp+agglomeration*)


(* ::Subsubsection::Closed:: *)
(*Initialization*)


patch=2;


ageDistributioncountryside={346119,744349,1339657,3190135,858195,560583,284205};
ageDistributionBp={129438,228716,478981,1138123,315001,190753,112523};
nAge=Length[ageDistributionBp];


ageDistributionpatch={ageDistributionBp,ageDistributioncountryside}; (*special for 2-patch model*)
popPatch=Table[Sum[ageDistributionpatch[[p1,k]],{k,1,nAge}],{p1,1,patch}];


contactMatrix=aggregatePREM[contactPREM];


alpha=Table[1/latentPeriod,{i,1,7}];
alpha2=Table[1/preSympInfPeriod,{i,1,7}];
gammaa=Table[1/infectPeriodA,{i,1,7}];
gammas=Table[1/infectPeriodS,{i,1,7}];


gammah=Table[1/infectPeriodH,{i,1,7}];
gammac=Table[1/infectPeriodC,{i,1,7}];
gammacr=Table[1/infectPeriodCR,{i,1,7}];


ksi=xi;
asymp=relInfectA;


R0Bp=1.7;
R0Countryside=1.5;


baseBetaBp=NGM[contactMatrix, R0Bp, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
baseBetaCountryside=NGM[contactMatrix, R0Countryside, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];


beta[t_]={baseBetaBp*contactMatrix, baseBetaCountryside*contactMatrix}; (*special for 2-patch model*)


l10Bp={0,0,1,1,1,0,0};
l10Countryside={0,0,0,0,0,0,0};
l20={0,0,0,0,0,0,0};
ip0={0,0,0,0,0,0,0};
ia10={0,0,0,0,0,0,0};
ia20={0,0,0,0,0,0,0};
ia30={0,0,0,0,0,0,0};
is10={0,0,0,0,0,0,0};
is20={0,0,0,0,0,0,0};
is30={0,0,0,0,0,0,0};
ih0={0,0,0,0,0,0,0};
ic0={0,0,0,0,0,0,0};
icr0={0,0,0,0,0,0,0};
r0={0,0,0,0,0,0,0};
d0={0,0,0,0,0,0,0};
c0={0,0,0,0,0,0,0};
s0Bp=ageDistributionBp-l10Bp;
s0Countryside=ageDistributioncountryside-l10Countryside;
initBp=Transpose[{s0Bp,l10Bp,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];
initCountryside=Transpose[{s0Countryside,l10Countryside,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];
initPatch={initBp,initCountryside}; (*special for 2-patch model*)


(* ::Subsubsection::Closed:: *)
(*Model: 400K travel/day*)


travelBpCountryside=4*101411;
travel=Table[Which[p1==p2,0,True,travelBpCountryside/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


solWeitz100k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfI[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.solWeitz100k];


weitz = Plot[sumOfI[t],{t,0,300},PlotRange->All,PlotStyle-> {Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],AxesLabel->{"time (days)","number of infected individuals"}]


plt=Plot[{Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.solWeitz100k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,2,2},{i,1,nAge}]/.solWeitz100k]},
{t,0,300},PlotRange->All,LabelStyle->Directive[Bold, 20],
AxesLabel->{"time (days)","number of infected individuals"},PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotLegends->Placed[{"Bp", "Countryside"},{Right, Right}],ImageSize->1500]


Export[StringJoin[generated,"hetero_Bp_400k.pdf"],plt]


(* ::Subsubsection::Closed:: *)
(*Model: 40K travel/day*)


travelBpCountryside=4*101411/(4*10);
travel=Table[Which[p1==p2,0,True,travelBpCountryside/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


solWeitz10k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfI[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.solWeitz10k];


Plot[sumOfI[t],{t,0,300},PlotRange->All]


Plot[{Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.solWeitz10k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,2,2},{i,1,nAge}]/.solWeitz10k]},
{t,0,300},PlotRange->All,LabelStyle->Directive[Bold, 20],
AxesLabel->{"time (days)","number of infected individuals"},PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotLegends->Placed[{"Bp", "Countryside"},{Right, Right}],ImageSize->1500]


plt=Plot[{Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.solWeitz10k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.solWeitz10k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,2,2},{i,1,nAge}]/.solWeitz10k]},{t,0,300},PlotRange->{0,All},
LabelStyle->Directive[Bold, 20],
AxesLabel->{"time (days)","number of infected individuals"},
PlotLegends->Placed[{"homogeneous model", "2-patch model with 100,000 daily travel","2-patch model with 10,000 daily travel"},{Right,Right}],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotRange->All,ImageSize->1500]


Export[StringJoin[generated,"hetero_Bp_40k.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Pest county*)


(* ::Subsubsection::Closed:: *)
(*Initialization*)


patch=2;


ageDistributionBp={152071,273328,549435,1316738,360795,216668,125293};
ageDistributioncountryside={323486,699737,1269203,3011520,812401,534668,271435};
nAge=Length[ageDistributionBp];


ageDistributionpatch={ageDistributionBp,ageDistributioncountryside}; (*special for 2-patch model*)
popPatch=Table[Sum[ageDistributionpatch[[p1,k]],{k,1,nAge}],{p1,1,patch}];


contactMatrix=aggregatePREM[contactPREM];


alpha=Table[1/latentPeriod,{i,1,7}];
alpha2=Table[1/preSympInfPeriod,{i,1,7}];
gammaa=Table[1/infectPeriodA,{i,1,7}];
gammas=Table[1/infectPeriodS,{i,1,7}];


gammah=Table[1/infectPeriodH,{i,1,7}];
gammac=Table[1/infectPeriodC,{i,1,7}];
gammacr=Table[1/infectPeriodCR,{i,1,7}];


ksi=xi;
asymp=relInfectA;


R0Bp=1.7;
R0Countryside=1.5;


baseBetaBp=NGM[contactMatrix, R0Bp, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
baseBetaCountryside=NGM[contactMatrix, R0Countryside, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];


beta[t_]={baseBetaBp*contactMatrix, baseBetaCountryside*contactMatrix}; (*special for 2-patch model*)


l10Bp={0,0,1,1,1,0,0};
l10Countryside={0,0,0,0,0,0,0};
l20={0,0,0,0,0,0,0};
ip0={0,0,0,0,0,0,0};
ia10={0,0,0,0,0,0,0};
ia20={0,0,0,0,0,0,0};
ia30={0,0,0,0,0,0,0};
is10={0,0,0,0,0,0,0};
is20={0,0,0,0,0,0,0};
is30={0,0,0,0,0,0,0};
ih0={0,0,0,0,0,0,0};
ic0={0,0,0,0,0,0,0};
icr0={0,0,0,0,0,0,0};
r0={0,0,0,0,0,0,0};
d0={0,0,0,0,0,0,0};
c0={0,0,0,0,0,0,0};
s0Bp=ageDistributionBp-l10Bp;
s0Countryside=ageDistributioncountryside-l10Countryside;
initBp=Transpose[{s0Bp,l10Bp,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];
initCountryside=Transpose[{s0Countryside,l10Countryside,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];
initPatch={initBp,initCountryside}; (*special for 2-patch model*)


(* ::Subsubsection::Closed:: *)
(*Model: 400K travel/day*)


travelBpCountryside=4*84789;
travel=Table[Which[p1==p2,0,True,travelBpCountryside/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


sol400k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfI400k[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol400k];


Plot[sumOfI400k[t],{t,0,300},PlotRange->All]


plt=Plot[{Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol400k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.sol400k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,2,2},{i,1,nAge}]/.sol400k]},
{t,0,350},PlotRange->{{0, 350},{0, 370000}},
LabelStyle->Directive[Bold, Black, 30], AxesStyle->Directive[Thickness[0.005]],
AxesLabel->{"time (days)","number of infected individuals"},PlotStyle->{
{CapForm["Butt"], Dashed, Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotLegends->Placed[{"cumulated epidemic curve","Pest county", "other parts"},{0.8, 0.8}],ImageSize->1500, Background->White]


Export[StringJoin[generated,"hetero_Pest_400k.pdf"],plt]


(* ::Subsubsection:: *)
(*Model: 10K travel/day*)


travelBpCountryside=4*84789/(4*10);
travel=Table[Which[p1==p2,0,True,travelBpCountryside/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


sol10k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfI10k[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol10k];


Plot[sumOfI10k[t],{t,0,300},PlotRange->All]


plt=Plot[{Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol10k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.sol10k],
Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,2,2},{i,1,nAge}]/.sol10k]},
{t,0,300},PlotRange->{{0, 350},{0, All}},
LabelStyle->Directive[Bold, Black, 30], AxesStyle->Directive[Thickness[0.005]], 
AxesLabel->{"time (days)","number of infected individuals"},PlotStyle->{
{CapForm["Butt"], Dashed, Thickness[0.01],ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}},
PlotLegends->Placed[{"cumulated epidemic curve","Pest county", "other parts"},{0.8, 0.8}],ImageSize->1500, Background->White]


Export[StringJoin[generated,"hetero_Pest_40k.pdf"],plt]


(* ::Section::Closed:: *)
(*What if ... impact of implemented measures*)


(* ::Subsection::Closed:: *)
(*Initial values*)


l10={0,0,0,1,0,0,0};
l20={0,0,0,0,0,0,0};
ip0={0,0,0,0,0,0,0};
ia10={0,0,0,0,0,0,0};
ia20={0,0,0,0,0,0,0};
ia30={0,0,0,0,0,0,0};
is10={0,0,0,0,0,0,0};
is20={0,0,0,0,0,0,0};
is30={0,0,0,0,0,0,0};
ih0={0,0,0,0,0,0,0};
ic0={0,0,0,0,0,0,0};
icr0={0,0,0,0,0,0,0};
r0={0,0,0,0,0,0,0};
d0={0,0,0,0,0,0,0};
c0={0,0,0,0,0,0,0};
s0=ageDistribution-(l10+l20+ip0+ia10);
init=Transpose[{s0,l10,l20,ip0,ia10,ia20,ia30,is10,is20,is30,ih0,ic0,icr0,r0,d0,c0}];


R0=2.2;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


(*ICU needs*)
Plot[Sum[(ic[i][t]),{i,1,nAge}]/.sol,{t,0,50},PlotRange->{0,All}]


(Table[ia1[i][t]+ia2[i][t]+ia3[i][t]+is1[i][t]+is2[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol)/.t->35


(*Epidemic curve*)
Plot[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],{t,0,50},PlotRange->{0,All}]


init1=Flatten[modelGenerator[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h][[1]]/.sol/.t->50,2]


(* ::Subsection::Closed:: *)
(*Model*)


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, init1, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


(*ICU needs*)
Plot[
Sum[ic[i][t],{i,1,nAge}]/.sol,{t,0,46},
PlotRange->{0,All},
LabelStyle->Directive[18,Black,Bold],
PlotLabel->Style["ICU Needs"],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["MidnightBlue"]}},
PlotRange->All,
Ticks->{{{0,"15.March"}, {16,"1.April"}, {31, "15.April"}, {46, "1.May"}},Automatic}]


(*Cumulative curve*)
plot=Show[Plot[
Evaluate[Sum[c[i][t],{i,1,nAge}]/.sol],{t,0,46},
PlotRange->{0,All},
LabelStyle->Directive[18,Black,Bold],
PlotLabel->Style["Number of cumulative cases"],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["MidnightBlue"]}},
PlotRange->All,
Ticks->{{{0,"15.March"}, {16,"1.April"}, {31, "15.April"}, {46, "1.May"}},Automatic}]]


hunData=20*{32,39,50,58,73,85,103,131,167,187,226,261,
300,343,408,447,492,525,585,623,678,733,744,817,895,
980,1190,1310,1410,1458,1512,1579,1652,1763,1834,1916,
1984,2098,2168,2284,2383,2443,2500,2583,2649,2727,2775};


func=Interpolation[Table[{i-1,hunData[[i]]},{i,1,47}]]


(*Cumulative curve*)
plot=Show[LogPlot[
{Evaluate[Sum[c[i][t],{i,1,nAge}]/.sol],func[t]},{t,0,46},
LabelStyle->Directive[35,Black,Bold], PlotRange->{{0,All},{20,All}}, AxesOrigin->{0.,20.},
PlotLabel->Style["Number of cumulative cases"],
PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["MidnightBlue"]},
Opacity[0]},
Filling->{1->{2}},
FillingStyle->Directive[ColorData["Crayola"]["BlizzardBlue"], Opacity[0.25]],
Ticks->{{{0,"15.March"}, {16,"1.April"}, {31, "15.April"}, {46, "1.May"}},{1,10,200,1000,{10000,10000},{100000,100000},{600000,600000}}}],
ListLogPlot[Table[{i-1,hunData[[i]]},{i,1,47}],PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["Blue"]}}],
ListLogPlot[Table[{i-1,1/20*hunData[[i]]},{i,1,47}],PlotStyle->{
{Thickness[0.01],ColorData["Crayola"]["Red"]}}],Background->White, ImageSize->{1000,500}]


Export[StringJoin[generated,"impactOfMeasures.pdf"],plot]


(*Mortality*)
Plot[Evaluate[mortalityCumulative/.sol],{t,0,46},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
Ticks->{{{0,"15.March"}, {16,"1.April"}, {31, "15.April"}, {46, "1.May"}},Automatic},
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}}]


(*Cumulative*)
Plot[Evaluate[cumulative/.sol],{t,0,45},
PlotRange->{0,All},
PlotStyle->{
{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
PlotLegends->ageStructureLegend, AxesStyle->Directive[Gray,FontSize->20], LabelStyle->Directive[Bold, 20],
Filling->{
{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}}]


(* ::Section::Closed:: *)
(*Technical appendix*)


(* ::Subsection::Closed:: *)
(*Contact matrix - heatmap*)


a = ColorData[{"ValentineTones", "Reverse"}];


b=	Function[
	x, (
		x; 
		x / (Max[contactMatrix] - Min[contactMatrix]) + 0.1 - Min[contactMatrix] / (Max[contactMatrix] - Min[contactMatrix])
		)
	];


c = Function[x, a[b[x]]];


contactMatrix = {{2.85265757236992`,1.0168035823185497`,0.9685513913908559`,3.6394710640997423`,0.7636165091831196`,0.26465811449632815`,0.16478042597673156`},{0.4920476243619219`,5.5591531856025815`,1.3504814369626075`,4.594274005138205`,0.7201168699258064`,0.5079625876327662`,0.2712341242084005`},{0.27685605410253394`,0.7977189385389726`,6.219897755428148`,6.522055041065205`,0.6110526839596531`,0.3066735514563918`,0.14958519657799538`},{0.4176563153332957`,1.0895009251108805`,2.618387301006593`,5.4218448320596275`,0.9159669127390052`,0.43426172617935965`,0.19979273798286082`},{0.26756270229416934`,0.5214145882147075`,0.7490271587498547`,2.7967218055259924`,1.6643410756692467`,0.6027761432448209`,0.17003271688765734`},{0.13848711461847765`,0.5492699095199908`,0.5613957491247491`,1.9801355696025458`,0.9001812106572092`,0.7311654280330745`,0.1455885636741889`},{0.16475717166363135`,0.5604191442260885`,0.5232340130971493`,1.7407548882408745`,0.4852001613787113`,0.2781901927944575`,0.5505304214784429`}}


Transpose[contactMatrix];


labelfuncH[s_] = Rotate[Style[s, 15], Pi/2];
labelfuncV[s_] = Style[s, 15];


AgelabelsH=Thread[{Range[7], {labelfuncH["0-4 "], labelfuncH["5-14 "], labelfuncH["15-29 "], labelfuncH["30-59 "], labelfuncH["60-69 "], labelfuncH["70-79 "], labelfuncH["80-     "]}}];
AgelabelsV=Thread[{Range[7], {labelfuncV["0-4"], labelfuncV["5-14"], labelfuncV["15-29"], labelfuncV["30-59"], labelfuncV["60-69"], labelfuncV["70-79"], labelfuncV["80-    "]}}];


contactMatrix
Max[contactMatrix]
Min[contactMatrix]
cntxmatrixplot = Labeled[MatrixPlot[Transpose[contactMatrix], ColorFunction->
	Function[
	x, (
		x; 
		a[x / (Max[contactMatrix] - Min[contactMatrix]) + 0.1 - Min[contactMatrix] / (Max[contactMatrix] - Min[contactMatrix])]
		)
	], ColorFunctionScaling->False, PlotLegends -> Automatic, DataReversed -> True, FrameTicks->{AgelabelsV, AgelabelsH}, ImageSize->Scaled[.3]], {Rotate[Style["Age of contact", 20, Bold], Pi/2], Style["Age of individual", 20, Bold]}, {Left, Bottom}, Spacings->{2, 1}]
(* Export["contactMatrixHeatMap.png",cntxmatrixplot] *)


Export[StringJoin[generated,"contactMatrix_heatMap.pdf"],cntxmatrixplot]
