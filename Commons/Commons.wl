(* ::Package:: *)

(* ::Section::Closed:: *)
(*Additive Groups*)


BeginPackage["Commons`"];



Commons::usage = "Methods used often but is not provided by mathematica. Aliases pointing to mathematica functions. Think it is intended to be that way.";


(* ::Section::Closed:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*Common tools*)


Subset::usage = "The ordinary proper-subset operator.";
SubsetEqual::usage = "The ordinary subset operator."; 
Congruent::usage = "Shorthand for the maping operator that maps two levels into a list instead of at level one.";
Pos::usage = " The fist position of the occurrence of an element in a list.";
NthLargest::usage = " The n:th next to largest in an unordered list.";

PrintIf::usage = "Write string to terminal depending on first boolean argument.";
TIMING::usage = "Boolean constant wether to enable time prognosis or not. Groups computation times grows as a power two polynomial with group order.";
EstimatedTime::usage = " Int --> Int  |  Group with subgroups's computation time grows as a power two polynomial with group order. This polynomial.";
EstimatedTimeCosets::usage = " Int --> Int  |  Group with subgroups and quotient groups's computation time including grows as a power two polynomial with group order. This polynomial.";


(* ::Section::Closed:: *)
(*Code*)


(* ::Subsection::Closed:: *)
(*Common tools - abstract *)


Subset[X_,Y_]:=ContainsAll[Y,X]\[And](Sort[X]!=Sort[Y]);

SubsetEqual[X_,Y_]:=ContainsAll[Y,X];

Congruent:=OperatorApplied[Commons`Private`MapLevel2,2];

Pos[x_,xs_]:=FirstPosition[xs,x][[1]];

NthLargest[n_,M_]:=(Sort@*DeleteDuplicates@*Flatten)[M][[-n]]


(* ::Subsection:: *)
(*Common tools - run-time*)


(* ::Subsection::Closed:: *)
(*Common tools - io*)


TIMING=True;

PrintIf[debug_,s_]:=If[debug,Print[s],None];

EstimatedTime[x_]:= 8.64*^-6 x^2;
EstimatedTimeCosets[x_]:= 1.02*^-8 x^2.75;


(* ::Subsection::Closed:: *)
(*Non-exports*)


Begin["`Private`"];

Commons`Private`MapLevel2[f_,L_]:=Map[f,L,{2}]

End[];


(* ::Section::Closed:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


Remove[f,s,x,L,X,Y,n,M,debug]
EndPackage[];
