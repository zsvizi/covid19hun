(* ::Package:: *)

BeginPackage["contactmatrix`"]


(*Age-specific population data*)
population1Y::usage = "population1Y"
ageDistribution::usage = "ageDistribution"


(*Contact matrices*)
ageDistribution::usage = "ageDistribution"
contactPREM::usage = "contactPREM"
contactPREMhome::usage = "contactPREMhome"
contactPREMschool::usage = "contactPREMschool"
contactPREMother::usage = "contactPREMother"
contactPREMwork::usage = "contactPREMwork"


(*Contact matrix transformators*)
aggregatePREM::usage = "aggregatePREM[matrix]"
transformMatrix::usage = "transformMatrix[matrix]"


Begin["`Private`"]


dataFolder=StringJoin[ToString[NotebookDirectory[]],"/data/"]


(* ::Subsubsection::Closed:: *)
(*Age distribution*)


(*Hungarian Statistical Office (KSH) data*)
population1Y=Import[StringJoin[dataFolder,"age_distribution.xls"]][[1]][[All,2]];
ageDistribution={
Sum[population1Y[[i]],{i,1,5}],
Sum[population1Y[[i]],{i,6,15}],
Sum[population1Y[[i]],{i,16,30}],
Sum[population1Y[[i]],{i,31,60}],
Sum[population1Y[[i]],{i,61,70}],
Sum[population1Y[[i]],{i,71,80}],
Sum[population1Y[[i]],{i,81,Length[population1Y]}]};


(*Age distribution with 5y bins+*)
population5Y = Flatten[{Table[Sum[population1Y[[i]],{i,5*j+1,5*j+5}],{j,0,14}],Sum[population1Y[[i]],{i,76,Length[population1Y]}]}];


(* ::Subsubsection::Closed:: *)
(*Contact matrices*)


(* ::Text:: *)
(*Prem et. al*)


sheets=Import[StringJoin[dataFolder,"contact_matrices.xls"], "Sheets"];
contactData=Import[StringJoin[dataFolder,"contact_matrices.xls"]];
sid=Association[Table[sheets[[i]]->i,{i,Range[Length[sheets]]}]];


contactPREM = contactData[[sid["all"]]];


contactPREMhome =contactData[[sid["home"]]];


contactPREMwork=contactData[[sid["work"]]];


contactPREMschool=contactData[[sid["school"]]];


contactPREMother=contactData[[sid["other"]]];


(* ::Subsubsection::Closed:: *)
(*Transformators*)


aggregatePREM[matrix_]:=Module[
{ageGroup80, ageGroup75, ratio75, numberOfContacts, bins, newMatrix},
ageGroup80=Sum[population1Y[[i]],{i,81,Length[population1Y]}];
ageGroup75=Sum[population1Y[[i]],{i,76,80}];
ratio75 = N[ageGroup75/(ageGroup75+ageGroup80)];
numberOfContacts=Table[(matrix[[i,j]]*population5Y[[i]]+matrix[[j,i]]*population5Y[[j]])/(2),{i,1,Length[matrix]},{j,1,Length[matrix]}];
bins={1,2,4,7,13,15,16};
newMatrix=Table[
Which[
(i!=j&&i<Length[bins]-1&&j<Length[bins]-1),Sum[numberOfContacts[[m,n]],{m,bins[[i]],bins[[i+1]]-1},{n,bins[[j]],bins[[j+1]]-1}],
(i==j&&i<Length[bins]-1&&j<Length[bins]-1),Sum[numberOfContacts[[m,n]],{m,bins[[i]],bins[[i+1]]-1},{n,m,bins[[j+1]]-1}],
(j==Length[bins]-1&&i<Length[bins]),Sum[numberOfContacts[[m,bins[[j]]]]+ratio75*numberOfContacts[[m,bins[[j+1]]]],{m,bins[[i]],bins[[i+1]]-1}],
(i==Length[bins]-1&&j<Length[bins]-1),Sum[numberOfContacts[[bins[[i]],n]]+ratio75*numberOfContacts[[bins[[i+1]],n]],{n,bins[[j]],bins[[j+1]]-1}],
(j==Length[bins]&&i<Length[bins]),Sum[(1-ratio75)*numberOfContacts[[m,bins[[j]]]],{m,bins[[i]],bins[[i+1]]-1}],
(i==Length[bins]&&j<Length[bins]),Sum[(1-ratio75)*numberOfContacts[[bins[[i]],n]],{n,bins[[j]],bins[[j+1]]-1}],
(j==Length[bins]&&i==Length[bins]),(1-ratio75)*numberOfContacts[[bins[[i]],bins[[j]]]],
True,0],
{i,1,Length[bins]},{j,1,Length[bins]}];
Table[(newMatrix[[i,j]]/ageDistribution[[i]]),{i,1,Length[newMatrix]},{j,1,Length[newMatrix]}]
]


transformMatrix[matrix_]:=Module[{numberOfContacts},
numberOfContacts=Table[(matrix[[i,j]]*population5Y[[i]]+matrix[[j,i]]*population5Y[[j]])/(2),{i,1,Length[matrix]},{j,1,Length[matrix]}];
Table[(numberOfContacts[[i,j]]/ageDistribution[[i]]),{i,1,Length[numberOfContacts]},{j,1,Length[numberOfContacts]}]
]


(* ::Subsubsection::Closed:: *)
(*End*)


End[ ]


EndPackage[ ]
