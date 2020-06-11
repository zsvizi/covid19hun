#!/usr/bin/env wolframscript -cloud
(* ::Package:: *)

Import[StringJoin[NotebookDirectory[],"modules/contactmatrix.wl"]]
Import[StringJoin[NotebookDirectory[],"modules/model.wl"]]


If[!DirectoryQ[StringJoin[NotebookDirectory[],"generated"]], CreateDirectory[StringJoin[NotebookDirectory[],"generated"]]];


generated=StringJoin[NotebookDirectory[],"generated/"];


(* ::Chapter:: *)
(*Initialization*)


(* ::Section::Closed:: *)
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


(* ::Section::Closed:: *)
(*Initial values*)


R0=1.35;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


l10={0,0,10,10,10,0,0};
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


(* ::Section::Closed:: *)
(*Functions*)


sumOfI[t_]:=Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.sol];


fdie[t_]:=Evaluate[Sum[d[i][t],{i,1,nAge}]/.sol];


fhosp[t_]:=Evaluate[Sum[ih[i][t],{i,1,nAge}]/.sol];


fincidence[s_]:=Evaluate[(Sum[c[i][t],{i,1,nAge}]/.sol)/.t->(s+1)]-Evaluate[(Sum[c[i][t],{i,1,nAge}]/.sol)/.t->s];


mortalityCumulative=Table[Sum[d[j][t],{j,1,i}],{i,1,nAge}];
cumulative=Table[Sum[c[j][t],{j,1,i}],{i,1,nAge}];


t0=t/.FindRoot[sumOfI[t]==5000,{t,150}][[1]]


die0=t/.FindRoot[fdie[t]==38,{t,150}][[1]]


(* ::Section::Closed:: *)
(*Finalizing the initial values*)


initialValuesImod=Flatten[modelGenerator[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h][[1]]/.sol/.t->die0,2]


R0=1.05;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


initialValuesImod=Round[modelGenerator[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h][[1]]/.sol/.t->24][[1,1]]


(* ::Chapter:: *)
(*Scenarios post-lockdown*)


(* ::Section::Closed:: *)
(*Styling*)


twoPlotsStyle={
{Thickness[0.01],ColorData["Crayola"]["Orange"]},
{Thickness[0.01],RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}};


twoPlotsStylehosp={
{Thickness[0.01],RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]},
{Thickness[0.01],ColorData["Crayola"]["Red"]}};


ageStructureLegend=LineLegend[{"0-4","5-14","15-29","30-59","60-69","70-79","80+"}, LegendLayout->"ReversedColumn", LegendMarkerSize->{40,12},LegendMarkers -> 
    Graphics[{AbsoluteThickness[3000], Line[{{-100, 100}, {0, 0}}]}]];


ticklabel1={{0,"1.May"},{31,"1.June"},{31+30,"1.July"},{31+30+31,"1.Aug"},{31+30+31+31,"1.Sept"},{31+30+31+31+30,"1.Oct"},{31+30+31+31+30+31,"1.Nov"},{31+30+31+31+30+31+30,"1.Dec"},{31+30+31+31+30+31+30+31,"1.Jan"}};
ticklabel2={{0,"1.May"},{31+30,"1.July"},{31+30+31+31,"1.Sept"},{31+30+31+31+30+31,"1.Nov"},{31+30+31+31+30+31+30+31,"1.Jan"},{31+30+31+31+30+31+30+31+31+28,"1.Feb"}};
ticklabel3={{0,"1.May"},{31+30+31+31,"1.Sept"},{31+30+31+31+30+31+30+31,"1.Jan"},{365,"1.May"},{365+31+30+30+31,"1.Sept"},{365+31+30+30+31+30+31+30+31,"1.Jan"}};


(* ::Section:: *)
(*Plots for Figure 6*)


(* ::Subsection:: *)
(*Plots for hospitalization*)


(* ::Subsubsection::Closed:: *)
(*R = 1.65*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Sum[ih[i][t] + icr[i][t],{i,1,nAge}]/.sol, 
	Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,200},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black,35],
		AxesStyle->Directive[Thickness[0.005]],
		PlotLabel->Style["Weak control, R = 1.65\nTotal hospital bed and ICU need"],
		PlotStyle->twoPlotsStylehosp,
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		Background->White,
		FrameTicks->{{0, 30, 60, 90, 120, 150,180},
					{{5000, "5,000 "},{10000, "10,000 "},{15000,"15,000 "},{20000,"20,000 "}}}, 
		ImageSize->900]


Export[StringJoin[generated,"weak_hosp.pdf"],plt]


(* ::Text:: *)
(*Calculations for Table 2*)


FindMaximum[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol,{t,70}]


FindMaximum[Sum[ic[i][t],{i,1,nAge}]/.sol,{t,70}]


FindRoot[Sum[ic[i][t],{i,1,nAge}]-1000/.sol,{t,30}]


(* ::Subsubsection::Closed:: *)
(*R = 1.32*)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol, 
	Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,300},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black,35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotLabel->Style["Moderate control, R = 1.32\nTotal hospital bed and ICU need"],
		PlotStyle->twoPlotsStylehosp,
		Background->White,
		FrameTicks->{{0, 60, 120,180,240,300},
					{{1500," 1,500 "},{3000," 3,000 "},{4500," 4,500 "},{6000," 6,000 "},{7500," 7,500 "}}}, 
		ImageSize->900]


Export[StringJoin[generated,"moderate_hosp.pdf"],plt]


(* ::Text:: *)
(*Calculations for Table 2*)


FindMaximum[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol,{t,100}]


FindMaximum[Sum[ic[i][t],{i,1,nAge}]/.sol,{t,100}]//Quiet


FindRoot[Sum[ic[i][t],{i,1,nAge}]-1000/.sol,{t,60}]


(* ::Subsubsection::Closed:: *)
(*R = 1.1*)


R0=1.1;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol, 
	Sum[ic[i][t],{i,1,nAge}]/.sol},{t,0,650},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[ Black,35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotLabel->Style["Strong control, R = 1.1\nTotal hospital bed and ICU need"],
		PlotStyle->twoPlotsStylehosp,
		PlotLegends->Placed[
					LineLegend[{"hospitalized\nindividuals","individuals\non ICU"},
					LegendMarkerSize->{40,12},
					LegendMarkers->Graphics[{AbsoluteThickness[3000], Line[{{-100, 100}, {0, 0}}]}]],{0.75, 0.6}],
		Background->White,
		FrameTicks->{{0, 120, 240, 360, 480, 600},
					{{250, "250 "},{500, "500 "},{750, "750 "},{1000, "1,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"strong_hosp.pdf"],plt]


(* ::Text:: *)
(*Calculations for Table 2*)


FindMaximum[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol,{t,100}]


FindMaximum[Sum[ic[i][t],{i,1,nAge}]/.sol,{t,100}]


(* ::Subsection:: *)
(*Plots for incidence*)


(* ::Subsubsection::Closed:: *)
(*R = 1.65*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	fincidence[s],{s,1,200},
		PlotRange->{0, All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black, 35],
		PlotLabel->Style["Incidence"],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->twoPlotsStyle,
		Background->White,
		FrameTicks->{{0, 30, 60, 90, 120, 150,180},
					{{20000,"20,000 "},{40000,"40,000 "},{60000,"60,000 "},{80000,"80,000 "},{100000,"100,000 "},{120000,"120,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"weak_inc.pdf"],plt]


(* ::Text:: *)
(*Calculations for Chapter 2.4.2*)


FindRoot[fincidence[s]-100000,{s,50}]


FindRoot[fincidence[s]-100000,{s,75}]


Sum[fincidence[s],{s,54,73,1}]


(* ::Subsubsection::Closed:: *)
(*R = 1.32*)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	fincidence[s],{s,1,300},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black,35],
		PlotLabel->Style["Incidence"],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->twoPlotsStyle,
		Background->White,
		FrameTicks->{{0, 60, 120, 180,240,300},
					{{10000,"10,000 "},{20000,"20,000 "},{30000,"30,000 "},{40000,"40,000 "},{50000,"50,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"moderate_inc.pdf"],plt]


(* ::Subsubsection::Closed:: *)
(*R = 1.1*)


R0=1.1;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	fincidence[s],{s,1,650},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black,35],
		PlotLabel->Style["Incidence"],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->twoPlotsStyle,
		Background->White,
		FrameTicks->{{0,120,240,360,480,600},
					{{1000,"1,000 "},{2000,"2,000 "},{3000,"3,000 "},{4000,"4,000 "},{5000,"5,000 "},{6000,"6,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"strong_inc.pdf"],plt]


(* ::Section:: *)
(*Plots for Figure 7*)


(* ::Subsection:: *)
(*Mortality*)


(* ::Subsubsection::Closed:: *)
(*R = 1.65*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[mortalityCumulative/.sol],{t,0,200},
		PlotLabel->Style["Mortality\nR=1.65"],
		PlotRange->{0,All},
		PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
					{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
					{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
					{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
					{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
					{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
					{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
		AxesStyle->Directive[Thickness[0.005], Black, FontSize->35], 
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
				{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
				{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
				{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
				{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
				{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
				{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
		Background->White,
		FrameTicks->{{0, 30, 60, 90, 120, 150,180},
					{{5000,"5,000 "},{10000,"10,000 "},{15000,"15,000 "},{20000,"20,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"weak_mort.pdf"],plt]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[mortalityCumulative[[7]]/.sol,{t,70}](*calculation for Table 2*)


(* ::Subsubsection::Closed:: *)
(*R = 1.32*)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[mortalityCumulative/.sol],{t,0,360},
	PlotRange->{0,All},
	PlotLabel->Style["Mortality\nR=1.32"],
	PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
	PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
				{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
				{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
				{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
				{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
				{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
				{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
	AxesStyle->Directive[Thickness[0.005],Black,FontSize->35],
	LabelStyle->Directive[ Black,35],
	Frame->{{True, False},{True, False}},
	FrameStyle-> Directive[Thickness[0.005], 35],
	FrameLabel->{"days passed",""},
	Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
			{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
			{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
			{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
			{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
			{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
			{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
	Background->White,
	FrameTicks->{{0, 60, 120, 180,240,300, 360},
				{{2000, "2,000 "},{4000, "4,000 "},{6000, "6,000 "},{8000, "8,000 "},{10000,"10,000 "},{12000, "12,000 "}}},
	ImageSize->900]


Export[StringJoin[generated,"moderate_mort.pdf"],plt]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[mortalityCumulative[[7]]/.sol,{t,70}](*calculation for Table 2*)


(* ::Subsubsection::Closed:: *)
(*R = 1.1*)


R0=1.1;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[mortalityCumulative/.sol],{t,0,630},
		PlotRange->{0,All},
		PlotLabel->Style["Mortality\nR=1.1"],
		PlotRangePadding->0,
		PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
					{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
					{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
					{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
					{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
					{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
					{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		LabelStyle->Directive[ Black,35],
		Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
				{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
				{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
				{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
				{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
				{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
				{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
		Background->White,
		FrameTicks->{{0,120,240,360,480,600},
					{{1000,"1,000 "},{2000,"2,000 "},{3000,"3,000 "},{4000,"4,000 "}}},
		ImageSize->900]


Export[StringJoin[generated,"strong_mort.pdf"],plt]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[mortalityCumulative[[7]]/.sol,{t,70}]//Quiet


(* ::Subsection:: *)
(*Recovered*)


(* ::Subsubsection::Closed:: *)
(*R = 1.65*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[cumulative/.sol],{t,0,200},
		PlotRange->{0,All},
		PlotLabel->Style["Recovered\nR=1.65"],
		PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
					{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
					{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
					{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
					{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
					{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
					{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
		AxesStyle->Directive[Thickness[0.005],Black,FontSize->35],
		LabelStyle->Directive[Black,35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
				{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
				{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
				{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
				{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
				{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
				{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
		Background->White,
		FrameTicks->{{0, 30, 60, 90, 120, 150,180},
					{{1000000,"1 Mio. "},{2000000,"2 Mio. "},{3000000,"3 Mio. "},{4000000,"4 Mio. "},{5000000,"5 Mio. "}}},
		ImageSize->900]


Export[StringJoin[generated,"weak_cumu.pdf"],plt]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[cumulative[[7]]/.sol,{t,70}](*calculation for Table 2*)


(* ::Subsubsection::Closed:: *)
(*R = 1.32*)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[cumulative/.sol],{t,0,360},
		PlotRange->{0,All},
		PlotLabel->Style["Recovered\nR=1.32"],
		PlotRangePadding->0,
		PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
					{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
					{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
					{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
					{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
					{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
					{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
		AxesStyle->Directive[Thickness[0.005],Black,FontSize->35],
		LabelStyle->Directive[ Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
				{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
				{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
				{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
				{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
				{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
				{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
		FrameTicks->{{0, 60, 120, 180,240,300, 360},
					{{1000000,"1 Mio. "},{2000000,"2 Mio. "},{3000000,"3 Mio. "}}},
		Background->White,
		ImageSize->900]


Export[StringJoin[generated,"moderate_cumu.pdf"],plt]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[cumulative[[7]]/.sol,{t,70}](*calculation for Table 2*)


(* ::Subsubsection::Closed:: *)
(*R = 1.1*)


R0=1.1;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	Evaluate[cumulative/.sol],{t,0,630},
		PlotRange->{0,All},
		PlotLabel->Style["Recovered\nR=1.1"],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		PlotStyle->{{RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`], Thickness[0.0075]},
					{RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`], Thickness[0.0075]},
					{RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`], Thickness[0.0075]},
					{RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`], Thickness[0.0075]},
					{RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`], Thickness[0.0075]},
					{RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`], Thickness[0.0075]},
					{RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`], Thickness[0.0075]}},
		AxesStyle->Directive[Thickness[0.005],Black,FontSize->35],
		LabelStyle->Directive[Black,35],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		Filling->{{7-> {{6}, RGBColor[0.06666666666666667`,0.06666666666666667`,0.09803921568627451`]}},
				{6-> {{5}, RGBColor[0.1764705882352941`,0.2823529411764706`,0.39215686274509803`]}},
				{5-> {{4}, RGBColor[0.8823529411764706`,0.2549019607843137`,0.07058823529411765`]}},
				{4-> {{3}, RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]}},
				{3-> {{2}, RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
				{2-> {{1}, RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]}},
				{1-> {Axis, RGBColor[0.615686274509804`,0.8196078431372549`,0.7333333333333333`]}}},
		FrameTicks->{{0,120,240,360,480,600},
					{{300000,"300,000 "},{600000,"600,000 "},{900000,"900,000 "},{1200000,"1.2 Mio. "}}},
		Background->White,
		ImageSize->900]


(* ::Text:: *)
(*Calculation for Table 2*)


FindMaximum[cumulative[[7]]/.sol,{t,70}]//Quiet


Export[StringJoin[generated,"strong_cumu.pdf"],plt]


(* ::Chapter:: *)
(*Age-specific analysis*)


(* ::Section:: *)
(*Plots for Figure 8*)


(* ::Subsection::Closed:: *)
(*Contact matrix*)


contactMatrixPREMSchool=
Table[
	Which[
	(4>=i>=2 || 4>=j>=2),
			contactPREM[[i,j]]-contactPREMschool[[i,j]]-0.5*contactPREMother[[i,j]],
	True,contactPREM[[i,j]]-contactPREMschool[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixSchool=aggregatePREM[contactMatrixPREMSchool];


(* ::Subsection::Closed:: *)
(*Hospital bed needs*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrixSchool};
solschool=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol, 
	Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solschool},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->twoPlotsStyle,
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}}, 
		FrameStyle -> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		PlotLabel->Style["Hospital bed needs"],
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		Background->White,
		FrameTicks->{{0,60,120,180,240, 300, 360, 420},
					{{5000, "5,000 "},{10000, "10,000 "},{15000, "15,000 "},{20000, "20,000 "}}}]


Export[StringJoin[generated,"school_hosp_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact *)


FindMaximum[
Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol],{t,100}]


FindMaximum[
Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solschool],{t,100}]


(* ::Subsection::Closed:: *)
(*ICU needs*)


plt=Plot[
	{Sum[ic[i][t],{i,1,nAge}]/.sol, 
	Sum[ic[i][t],{i,1,nAge}]/.solschool},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->twoPlotsStyle,
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle -> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		PlotLabel->Style["ICU needs"],
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		Background->White
		FrameTicks->{{0,60,120,180,240, 300, 360, 420},
					{{2000, "2000 "},{4000, "4,000 "},{6000, "6,000 "},{8000, "8,000 "}}}]


Export[StringJoin[generated,"school_icu_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact *)


FindMaximum[
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.sol],{t,60}](*calculation for impact of school closures*)


FindMaximum[
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solschool],{t,100}](*calculation for impact of school closures*)


(* ::Subsection::Closed:: *)
(*Mortality*)


plt=Plot[
	{Evaluate[Total[mortalityCumulative[[7]]]/.sol],
	Evaluate[Total[mortalityCumulative[[7]]]/.solschool]},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->twoPlotsStyle,
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle -> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		PlotLabel->Style["Mortality"],
		Background->White,
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		FrameTicks->{{0,60,120,180,240, 300, 360, 420, 480},
					{{10000, "10,000 "},{20000, "20,000 "},{30000, "30,000 "},{40000, "40,000 "}}}]


Export[StringJoin[generated,"school_mort_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact *)


FindMaximum[mortalityCumulative[[7]]/.sol,{t,70}](*calculation for impact of measure*)


FindMaximum[mortalityCumulative[[7]]/.solschool,{t,70}](*calculation for impact of measure*)


(* ::Section:: *)
(*Plots for Figure 9*)


(* ::Subsection::Closed:: *)
(*Contact matrix*)


contactMatrixPREMElderlyFull=
Table[
	Which[
	(i>=13 || j>=13),
		contactPREM[[i,j]]-contactPREMother[[i,j]]-contactPREMschool[[i,j]]-contactPREMwork[[i,j]],
	True,
		contactPREM[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixElderlyFull=aggregatePREM[contactMatrixPREMElderlyFull];


contactMatrixPREMElderlyModerate=
Table[
	Which[(i>=13 || j>=13),
		contactPREM[[i,j]]-contactPREMschool[[i,j]]-0.5*contactPREMother[[i,j]]-0.5*contactPREMwork[[i,j]],
	True,
		contactPREM[[i,j]]],
{i,1,Length[contactPREM]},{j,1,Length[contactPREM]}];
contactMatrixElderlyModerate=aggregatePREM[contactMatrixPREMElderlyModerate];


(* ::Subsection::Closed:: *)
(*Hospital bed needs*)


R0=1.65;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*contactMatrixElderlyModerate};
solElderlyModerate=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*contactMatrixElderlyFull};
solElderlyFull=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol],
	Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solElderlyModerate],
	Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solElderlyFull]},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["OuterSpace"]},
					{Thickness[0.01],ColorData["Crayola"]["Silver"]},
					{Thickness[0.01],ColorData["Crayola"]["Periwinkle"]}},
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle -> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		PlotLabel->Style["Hospital bed needs"],
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		Background->White,
		FrameTicks->{{0,60,120,180,240},
					{{5000, "5,000 "},{10000, "10,000 "},{15000, "15,000 "},{20000, "20,000 "}}}]


Export[StringJoin[generated,"elderly_total_hosp_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact*)


FindMaximum[
Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.sol],{t,60}]//Quiet


FindMaximum[
Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solElderlyModerate],{t,60}]//Quiet


FindMaximum[
Evaluate[Sum[ih[i][t]+icr[i][t],{i,1,nAge}]/.solElderlyFull],{t,60}]//Quiet


(* ::Subsection::Closed:: *)
(*ICU needs*)


plt=Plot[
	{Evaluate[Sum[ic[i][t],{i,1,nAge}]/.sol],
	Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyModerate],
	Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyFull]},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["OuterSpace"]},
					{Thickness[0.01],ColorData["Crayola"]["Silver"]},
					{Thickness[0.01],ColorData["Crayola"]["Periwinkle"]}},
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle -> Directive[Thickness[0.005], 35], 
		FrameLabel->{"days passed",""},
		PlotLabel->Style["ICU needs"],
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		Background->White,
		FrameTicks->{{0,60,120,180,240},
					{{1500, "1,500 "},{3000, "3,000 "},{4500, "4,500 "},{6000, " 6,000 "}}}]


Export[StringJoin[generated,"elderly_hosp_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact*)


FindMaximum[
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.sol],{t,60}]//Quiet


FindMaximum[
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyModerate],{t,60}]//Quiet


FindMaximum[
Evaluate[Sum[ic[i][t],{i,1,nAge}]/.solElderlyFull],{t,60}]//Quiet


(* ::Subsection::Closed:: *)
(*Mortality*)


plt=Plot[
	{Evaluate[Total[mortalityCumulative[[7]]]/.sol],
	Evaluate[Total[mortalityCumulative[[7]]]/.solElderlyModerate],
	Evaluate[Total[mortalityCumulative[[7]]]/.solElderlyFull]},{t,0,250},
		PlotRange->{0,All},
		PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["OuterSpace"]},
					{Thickness[0.01],ColorData["Crayola"]["Silver"]},
					{Thickness[0.01],ColorData["Crayola"]["Periwinkle"]}},
		AxesStyle->Directive[Thickness[0.005], Black, FontSize -> 35],
		LabelStyle->Directive[Black, 35],
		Frame->{{True, False},{True, False}},
		FrameStyle -> Directive[Thickness[0.005], 35],
		FrameLabel->{"days passed",""},
		PlotLabel->Style["Mortality"],
		Background->White,
		FrameLabel->{"days since release",""},
		ImageSize -> 900,
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		FrameTicks->{{0,60,120,180,240},
					{{10000, "10,000 "},{20000, "20,000 "},{30000, "30,000 "},{40000, "40,000 "},{50000, "50,000 "}}}]


Export[StringJoin[generated,"elderly_mort_1.65.pdf"],plt]


(* ::Text:: *)
(*Calculation for determining the impact*)


FindMaximum[mortalityCumulative[[7]]/.sol,{t,70}](*calculation for impact of measure*)


FindMaximum[mortalityCumulative[[7]]/.solElderlyModerate,{t,70}](*calculation for impact of measure*)


FindMaximum[mortalityCumulative[[7]]/.solElderlyFull,{t,70}](*calculation for impact of measure*)


(* ::Chapter:: *)
(*Seasonality - Plots for Figure 10*)


(* ::Section::Closed:: *)
(*Styling*)


fourplotsStyle = {{Thickness[0.01],
ColorData["Crayola"]["VividViolet"]},
{Thickness[0.01],ColorData["Crayola"]["VioletBlue"]},
{Thickness[0.01],ColorData["Crayola"]["BlueGreen"]},
{Thickness[0.01],ColorData["Crayola"]["JungleGreen"]}};


(* ::Section::Closed:: *)
(*Seasonality function*)


maxChange=0.3;
period[change_,t_]:=0.5*change*Cos[(t*2*Pi)/366]+(1-0.5*change);
inset=Plot[
		{1,period[0.1,t],
		period[0.2,t],
		period[0.3,t]},{t,0,366},
			PlotRange->{1-1.16*maxChange,1.05},
			LabelStyle->Directive[Black, 30],
			PlotStyle->fourplotsStyle,
			PlotLabel->Style["Seasonality function", 35],
			Frame->{{True, False},{True, False}},
			FrameStyle-> Directive[Thickness[0.005], 35],
			FrameLabel->{"","rel. transmissibility"},
			FrameTicks->{{{0,"December 1"},{31+31+28+31+30+31,"June 1"},{366,"December 1"}},
						{0.0,0.1,0.2,0.3,0.4,0.5,0.65,{0.7, "0.7 "},{0.8, "0.8 "},{0.9, "    0.9 "},{1.0, "1 "}}},
			ImageSize->900,
			Background->White]


Export[StringJoin[generated,"seas_inset.pdf"],inset]


(* ::Section::Closed:: *)
(*R = 1.65*)


R0=1.65;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]},{t,0,250},
		PlotRange->{0,All},
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black, 30],
		AxesLabel->{"time","number of infected individuals"},
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"","no. of infected"},
		PlotLabel->Style["Epidemic curve with seasonality, R = 1.65"],
		PlotStyle->fourplotsStyle,
		PlotRange->All,
		FrameTicks->{{{0,"May 1"},{61,"July 1"},{123, "Sept. 1"},{184, "Nov. 1"},{245,"Jan. 1"},300},
					{{100000,"100,000 "},{200000,"200,000 "},{300000,"300,000 "},{400000,"400,000 "},{500000,"500,000 "}}},
		ImageSize->900,
		Background->White]


Export[StringJoin[generated,"seas_weak.pdf"],plt]


(* ::Text:: *)
(*Calculation for impact of seasonality*)


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP03],{t,150}]


(* ::Section::Closed:: *)
(*R = 1.32*)


R0=1.32;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]},{t,0,400},
		LabelStyle->Directive[Black, 30],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		PlotLabel->Style["Epidemic curve with seasonality, R = 1.32"],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"","no. of infected"},
		PlotStyle->fourplotsStyle,
		PlotRange->{0,All},
		FrameTicks->{{{0,"May 1"}, {92,"Aug. 1"},{184,"Nov. 1"},{276, "Feb. 1"},{366, "May 1"},450},
					{{50000,"50,000 "},{100000,"100,000 "},{150000,"150,000 "},{200000,"200,000 "}}}, 
		ImageSize->900, 
		Background->White]


Export[StringJoin[generated,"seas_moderate.pdf"],plt]


(* ::Text:: *)
(*Calculation for impact of seasonality*)


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP01],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP03],{t,150}]


(* ::Section::Closed:: *)
(*R = 1.1*)


R0=1.1;
t0501=31+31+28+31+30;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.1,t+t0501]*contactMatrix};
solP01=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.2,t+t0501]*contactMatrix};
solP02=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


beta[t_]={baseBeta*period[0.3,t+t0501]*contactMatrix};
solP03=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


plt=Plot[
	{Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.sol],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP01],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP02],
	Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t],{i,1,nAge}]/.solP03]},{t,0,700},
		PlotRange->{0,All},
		LabelStyle->Directive[Black, 30],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"","no. of infected"},
		PlotLabel->Style["Epidemic curve with seasonality, R = 1.1"],
		AxesLabel->{"time","number of infected individuals"},
		PlotStyle->fourplotsStyle,
		PlotRange->All,
		FrameTicks->{{{0,"May 1"},{123, "Sept. 1"},{245,"Jan. 1"},{366, "May 1"},{489, "Sept. 1"},{612, "Jan. 1"}},
					{{5000,"5,000 "},{10000,"10,000 "},{15000, "15,000 "},{20000, "20,000 "},{25000, " 25,000 "}}},
		ImageSize->900, 
		Background->White]


Export[StringJoin[generated,"seas_strong.pdf"],plt]


(* ::Text:: *)
(*Calculation for impact of seasonality*)


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP01],{t,150}]//Quiet


FindMaximum[Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.solP02],{t,150}]//Quiet


(* ::Chapter:: *)
(*Spatial heterogeneity*)


(* ::Section:: *)
(*Plots for Figure 11*)


(* ::Subsection::Closed:: *)
(*Initialization *)


R0=1.32;


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


ageDistributionotherO={308442,673820,1107068,2748640,896129,616662,310281};
ageDistributionBp = ageDistribution - ageDistributionotherO;
ageDistributionother = ageDistributionotherO;
sol1patch=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr,ageDistribution, mu, ksi, asymp, pA, h];


sumOfI1patch[t_]:=Evaluate[Sum[ia1[i][t]+is1[i][t]+ia2[i][t]+is2[i][t]+ia3[i][t]+is3[i][t]+ip[i][t]+icr[i][t],{i,1,nAge}]/.sol1patch];


patch=2;


nAge=Length[ageDistributionBp];


ageDistributionpatch={ageDistributionBp,ageDistributionother}; (*special for 2-patch model*)
popPatch=Table[Total[ageDistributionpatch[[p1]]],{p1,1,patch}];


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


(* ::Subsection::Closed:: *)
(*Initial values*)


R0Bp=1.32;
R0other=1.32;


baseBetaBp=NGM[contactMatrix, R0Bp, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
baseBetaother=NGM[contactMatrix, R0other, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];


beta[t_]={baseBetaBp*contactMatrix, baseBetaother*contactMatrix}; (*special for 2-patch model*)


initBp = initialValuesImod*0.55;
transposeBp = Transpose[initBp];
initother = initialValuesImod*0.45;
transposeother = Transpose[initother];


{s0Bp, l10Bp, l20Bp,ip0Bp, ia10Bp, ia20Bp, ia30Bp, is10Bp, is20Bp, is30Bp, ih0Bp, ic0Bp, icr0Bp, r0Bp, d0Bp, c0Bp} = Table[initValue, {initValue,transposeBp}];


{s0other, l10other, l20other,ip0other, ia10other, ia20other, ia30other, is10other, is20other, is30other, ih0other, ic0other, icr0other, r0other, d0other, c0other} = Table[initValue, {initValue,transposeother}];


s0Bp=ageDistributionBp-Sum[transposeBp[[i]],{i,2,15}];
s0other=ageDistributionother-Sum[transposeother[[i]],{i,2,15}];
initBudapest=Transpose[{s0Bp,l10Bp,l20Bp,ip0Bp,ia10Bp,ia20Bp,ia30Bp,is10Bp,is20Bp,is30Bp,ih0Bp,ic0Bp,icr0Bp,r0Bp,d0Bp,c0Bp}];
initother=Transpose[{s0other,l10other,l20other,ip0other,ia10other,ia20other,ia30other,is10other,is20other,is30other,ih0other,ic0other,icr0other,r0other,d0other,c0other}];
initPatch={initBudapest,initother};


travelBpother=400000;
travel=Table[Which[p1==p2,0,True,travelBpother/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


sol400k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfI400k[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol400k];


(* ::Subsection::Closed:: *)
(*First scenario*)


plt=Plot[
	{sumOfI1patch[t],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.sol400k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.sol400k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,patch,patch},{i,1,nAge}]/.sol400k]},{t,0,250},
		PlotRange->{0, All},
		PlotLabel->Style[""],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black, 35],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->{{CapForm["Butt"], Dashed, Thickness[0.01],RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]},
					{CapForm["Butt"], Dashed, Thickness[0.01],ColorData["Crayola"]["Brown"]},
					{Thickness[0.01],RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]},
					{Thickness[0.01],RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
		ImageSize->900,
		FrameTicks->{{0, 30, 60, 90, 120, 150, 180, 210, 240},
					{{50000,"50,000 "},{100000,"100,000 "},{150000,"150,000 "},{200000,"200,000 "}}},
		Background->White]


Export[StringJoin[generated,"hetero_400k_BP1.32Country1.32Uniform1.32.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Second scenario*)


travelBpother=10000;
travel=Table[Which[p1==p2,0,True,travelBpother/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


solequalR10k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfIequalR10k[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.solequalR10k];


plt=Plot[
	{sumOfI1patch[t],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.solequalR10k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.solequalR10k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,patch,patch},{i,1,nAge}]/.solequalR10k]},{t,0,250},
		PlotRange->{0, All},
		PlotLabel->Style[""],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black, 35],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->{{CapForm["Butt"], Dashed, Thickness[0.01],RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]},
					{CapForm["Butt"], Dashed, Thickness[0.01],ColorData["Crayola"]["Brown"]},
					{Thickness[0.01],RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]},
					{Thickness[0.01],RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
		ImageSize->900,
		FrameTicks->{{0, 30, 60, 90, 120, 150, 180, 210, 240},
					{{50000,"50,000 "},{100000,"100,000 "},{150000,"150,000 "},{200000,"200,000 "}}},
		Background->White]


Export[StringJoin[generated,"hetero_10k_BP1.32Country1.32Uniform1.32.pdf"],plt]


(* ::Subsection::Closed:: *)
(*Third scenario*)


R0Bp=1.47;
R0other=1.25;


baseBetaBp=NGM[contactMatrix, R0Bp, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
baseBetaother=NGM[contactMatrix, R0other, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];


beta[t_]={baseBetaBp*contactMatrix, baseBetaother*contactMatrix}; (*special for 2-patch model*)


travelBpother=10000;
travel=Table[Which[p1==p2,0,True,travelBpother/popPatch[[p1]]],{p1,1,patch},{p2,1,patch},{k,1,nAge}]; (*special for 2-patch model*)


soldifferentR10k=modelSolver[patch, travel, beta, initPatch, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistributionpatch, mu, ksi, asymp, pA, h];


sumOfIdifferentR10k[t_]:=Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.soldifferentR10k];


plt=Plot[
	{sumOfI1patch[t],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,patch},{i,1,nAge}]/.soldifferentR10k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,1,1},{i,1,nAge}]/.soldifferentR10k],
	Evaluate[Sum[ia1[p,i][t]+is1[p,i][t]+ia2[p,i][t]+is2[p,i][t]+ia3[p,i][t]+is3[p,i][t]+ip[p,i][t],{p,patch,patch},{i,1,nAge}]/.soldifferentR10k]},{t,0,250},
		PlotRange->{0, All},
		PlotLabel->Style[""],
		PlotRangePadding->{{0, 0}, {0, Scaled[0.02]}},
		LabelStyle->Directive[Black, 35],
		AxesStyle->Directive[Thickness[0.005]],
		Frame->{{True, False},{True, False}},
		FrameStyle-> Directive[Thickness[0.005], 35],
		FrameLabel->{"days since release",""},
		PlotStyle->{{CapForm["Butt"], Dashed, Thickness[0.01],RGBColor[0.2549019607843137`,0.7176470588235294`,0.7058823529411764`]},
					{CapForm["Butt"], Dashed, Thickness[0.01],ColorData["Crayola"]["Brown"]},
					{Thickness[0.01],RGBColor[0.9372549019607843`,0.4627450980392157`,0.19215686274509802`]},
					{Thickness[0.01],RGBColor[0.9568627450980391`,0.7607843137254902`,0.3254901960784314`]}},
		ImageSize->900,
		FrameTicks->{{0, 30, 60, 90, 120, 150, 180, 210, 240},
					{{50000,"50,000 "},{100000,"100,000 "},{150000,"150,000 "},{200000,"200,000 "}}},
		Background->White]


Export[StringJoin[generated,"hetero_10k_BP1.47Other1.25Uniform1.32.pdf"],plt]


(* ::Chapter:: *)
(*Impact of implemented measures*)


(* ::Section::Closed:: *)
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


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


init1=Flatten[modelGenerator[beta, init, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h][[1]]/.sol/.t->45,2];
(*get appropriate time to fit to the initial phase of the epidemic in Hungary*)


(* ::Section::Closed:: *)
(*Plot - Figure 12 (left)*)


baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp]
beta[t_]={baseBeta*contactMatrix};


sol=modelSolver[beta, init1, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];


underreporting=17;


hunData=underreporting*{32,39,50,58,73,85,103,131,167,187,226,261,
300,343,408,447,492,525,585,623,678,733,744,817,895,
980,1190,1310,1410,1458,1512,1579,1652,1763,1834,1916,
1984,2098,2168,2284,2383,2443,2500,2583,2649,2727,2775,2863,2942,2998,3035,3065,3111,3150,3178,3213,3263,3284};


func=Interpolation[Table[{i-1,hunData[[i]]},{i,1,57}]]


(*Cumulative curve*)
plot=Show[LogPlot[
			{Evaluate[Sum[c[i][t],{i,1,nAge}]/.sol],func[t]},{t,0,56},
				LabelStyle->Directive[35,Black,Bold], 
				PlotRange->{{0,All},{20,All}}, 
				AxesOrigin->{0.,20.},
				PlotLabel->Style["Cumulative number of cases"],
				PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["MidnightBlue"]}, Opacity[0]},
				Filling->{1->{2}},
				FillingStyle->Directive[ColorData["Crayola"]["BlizzardBlue"], Opacity[0.25]],
				Ticks->{{{0,"March 15"}, {16,"April 1"}, {31, "April 15"}, {46, "May 1"},{56, "May 10"}},
						{1,10,100,1000,{10000,10000},{100000,100000},{1000000,1000000},{10000000,10000000}}}],
		  ListLogPlot[
			Table[{i-1,hunData[[i]]},{i,1,57}],
				PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["Blue"]}}],
		  ListLogPlot[
			Table[{i-1,1/underreporting*hunData[[i]]},{i,1,57}],
				PlotStyle->{{Thickness[0.01],ColorData["Crayola"]["Red"]}}],
		Background->White, 
		ImageSize->{1000,500}]


Export[StringJoin[generated,"impactOfMeasures.pdf"],plot]


(* ::Chapter:: *)
(*Sensitivity analysis (for generate this plot, set Boolean to True!)*)


(* ::Section::Closed:: *)
(*Data generation for plotting*)


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
sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h];
attack=(Flatten[(Table[c[i][t],{i,1,nAge}]/.sol)/.t->300]/ageDistribution);
If[runThisSection,
	rr0=1;
	For[Rt=R0init,Rt<=R0end,Rt+=1/stepSize,
		hh=1;
		baseBeta=NGM[contactMatrix, R0, nAge, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, pA, asymp];
		beta[t_]={baseBeta*(Rt/R0)*contactMatrix};
		For[hmult=hmultInit,hmult<=hmultEnd,hmult+=1/stepSize,
			sol=modelSolver[beta, initialValuesImod, alpha, alpha2, gammaa, gammas, gammah, gammac, gammacr, ageDistribution, mu, ksi, asymp, pA, h*hmult];
			target[t_]:=Sum[ic[i][t],{i,1,nAge}]/.sol;
			ifrr=Total[ageDistribution*attack*(1-pA)*(h*hmult)*xi*mu]/Total[ageDistribution*attack];
			maximumList[[hh,rr0]]={ifrr, Rt, Quiet[FindMaximum[target[t],{t,80}]][[1]]};
			hh++;
		];
		rr0++;
	]
]//Timing


IFRdistinct=0; maximumDistinct=0;
distinctVector=maximumList[[1+Round[(1.0-hmultInit)*stepSize],1+Round[(1.65-R0init)*stepSize]]];


(* ::Section::Closed:: *)
(*Plot - Figure 12 (right)*)


legendMax=40000;
plot=Show[
		ListContourPlot[
			Flatten[maximumList,1],
				PlotRange->{{0.00289,0.0089},{1.1,2.4}},
				Contours->({#, If[Round[#] == 8000,
								Directive[Thickness[0.01], White],
								Null
								]}& /@ Range[ 0, legendMax, 1000]),
				PlotLegends->BarLegend[{Automatic, {0, legendMax}}, 
										Table[i,{i,0,40000,1000}]], 
				ColorFunction->ColorData[{"RedBlueTones", {legendMax,0}}],
				ColorFunctionScaling -> False,
				FrameLabel->{"IFR",
							 Style["\[ScriptCapitalR]",FontSize->15]},
							 LabelStyle->Directive[Bold, 12, Black]],
		Graphics[
			{White, PointSize[0.025], Point[{0.00393415,1.695}]}]]


Export[StringJoin[generated,"contourPlot.pdf"],plot]


(* ::Chapter:: *)
(*Technical appendix*)


(* ::Section:: *)
(*Contact matrix - heatmap*)


(* ::Subsection::Closed:: *)
(*Data and functions*)


a=ColorData[{"ValentineTones", "Reverse"}];


b=Function[x, (x; x / (Max[contactMatrix] - Min[contactMatrix]) + 0.1 - Min[contactMatrix] / (Max[contactMatrix] - Min[contactMatrix]))];


c=Function[x, a[b[x]]];


contactMatrix = 
{{2.85265757236992`,1.0168035823185497`,0.9685513913908559`,3.6394710640997423`,0.7636165091831196`,0.26465811449632815`,0.16478042597673156`},
{0.4920476243619219`,5.5591531856025815`,1.3504814369626075`,4.594274005138205`,0.7201168699258064`,0.5079625876327662`,0.2712341242084005`},
{0.27685605410253394`,0.7977189385389726`,6.219897755428148`,6.522055041065205`,0.6110526839596531`,0.3066735514563918`,0.14958519657799538`},
{0.4176563153332957`,1.0895009251108805`,2.618387301006593`,5.4218448320596275`,0.9159669127390052`,0.43426172617935965`,0.19979273798286082`},
{0.26756270229416934`,0.5214145882147075`,0.7490271587498547`,2.7967218055259924`,1.6643410756692467`,0.6027761432448209`,0.17003271688765734`},
{0.13848711461847765`,0.5492699095199908`,0.5613957491247491`,1.9801355696025458`,0.9001812106572092`,0.7311654280330745`,0.1455885636741889`},
{0.16475717166363135`,0.5604191442260885`,0.5232340130971493`,1.7407548882408745`,0.4852001613787113`,0.2781901927944575`,0.5505304214784429`}};


Transpose[contactMatrix];


labelfuncH[s_] = Rotate[Style[s, 15], Pi/2];
labelfuncV[s_] = Style[s, 15];


AgelabelsH=Thread[{Range[7], {labelfuncH["0-4 "], labelfuncH["5-14 "], labelfuncH["15-29 "], labelfuncH["30-59 "], labelfuncH["60-69 "], labelfuncH["70-79 "], labelfuncH["80-     "]}}];
AgelabelsV=Thread[{Range[7], {labelfuncV["0-4"], labelfuncV["5-14"], labelfuncV["15-29"], labelfuncV["30-59"], labelfuncV["60-69"], labelfuncV["70-79"], labelfuncV["80-    "]}}];


(* ::Subsection::Closed:: *)
(*Plot - Figure C1.*)


cntxmatrixplot = 
Labeled[
	MatrixPlot[
		Transpose[contactMatrix], 
			ColorFunction->Function[x, (x; a[x / (Max[contactMatrix] - Min[contactMatrix]) + 0.1 - Min[contactMatrix] / (Max[contactMatrix] - Min[contactMatrix])];)], 
			ColorFunctionScaling->False, 
			PlotLegends -> Automatic, 
			DataReversed -> True, 
			FrameTicks->{AgelabelsV, AgelabelsH}, 
			ImageSize->Scaled[.3]], 
			{Rotate[Style["Age of contact", 20, Bold], Pi/2], Style["Age of individual", 20, Bold]}, 
			{Left, Bottom}, 
			Spacings->{2, 1}
	]


Export[StringJoin[generated,"contactMatrix_heatMap.pdf"],cntxmatrixplot];
