(* ::Package:: *)

(* ::Section:: *)
(*Additive Group Minimal Package*)


PrependTo[$ContextPath,"Commons`"];
BeginPackage["AdditiveGroupMinimal`"];
<< Commons`


AdditiveGroupMinimalPackage::usage = "This is the first module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suite. The AdditiveGroup package suite. " <>
   									"It contains the minimal functionality needed for investigating \!\(\*SubscriptBox[\(Z\), \(n\)]\) and it's subgroups.";

Print["AdditivegroupMinimal`: See Docs[\"Minimal\"] for documentation."];


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Constants, operators, constructors  and primitive mappings*)


N0::usage = " Int  \n  Modulus of the operators and size of \!\(\*SubscriptBox[\(Z\), \(n\)]\) that is n."; 
N1::usage = " Int  \n  Number of subgroups in \!\(\*SubscriptBox[\(Z\), \(n\)]\). Modulus of the total quotient group operator."; 

Zn::usage = " {Int}  \n  List of the elements in \!\(\*SubscriptBox[\(Z\), \(n\)]\) the additive group modulo N0.";
Sns::usage = " {{Int}}  \n  Subgroups of the group currently in use.";

CirclePlus::usage = " Int,Int --> Int  \n  Group operator (addition).";
SuperMinus::usage = " Int --> Int  \n  The inverse of an element modulo N0.";

MakeMinimalGroup::usage = " Int --> Global`  \n  Computes an additive group of order n and also precomputes it's subgroups. All variables it computes is global to current context.";


(* ::Subsection:: *)
(*Independent instances*)


MakeMinimalGroupInstance::usage = " Int --> {{Int},Int,Int,Int,{Int},{Int}}  \n  Computes an additive group instance containing also it's subgroups, though not it's quotient groups." <>
      																			   " It temporarily alters but restores the current context. It returns a 'context' given as a six-tuple" <>
      																			   " {Zn,N0,N1,N2,Sns,Css} where N2 and Css are null in this module. See MakeGroupInstance in Quotients package." <>
      																			   "See N0. See N1. See Sns.";
InstanceSubgroups::usage = " Int -->  {{Int}}  \n  Returns the subgroups of the additive group of a given order. It temporarily alters but restores the current context group.";


(* ::Subsection:: *)
(*Elementwise*)


ElementOrder::usage = " Int--> Int  \n The order of an element in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
ElementOrders::usage = " {Int}  \n  All the orders of the elements in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";


(* ::Subsection:: *)
(*Subsets and related elements*)


Subgroups::usage = " {{Int}}  \n  The subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
SubgroupsAndGenerator::usage = " <|Int -> {Int}|>   \n  Mapping from the generators to the subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) .";
SubgroupGenerators::usage = " {Int}  \n  The generators of all subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) given in subgroup size order.";

Zeros::usage = " The zeros \!\(\*SubscriptBox[\(ofZ\), \(n\)]\) given as pairs.";


(* ::Subsection:: *)
(*Structure and graphical overview*)


CayleyTable::usage = " []|{Int},(Int,Int -> Int) --> Grid[Int]  \n  Multiplication table of \!\(\*SubscriptBox[\(Z\), \(n\)]\) or of some given group and operator.";

Docs::usage = "Documentation of a package as an association between method names and their usage descriptions. Argument one of \"Minimal\",\"Basic\",\"\", \"Quotients\",\"Theorems\".";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Constants, operators, constructors  and primitive mappings*)


CirclePlus:= AdditiveGroupMinimal`Private`ModularAddition

SuperMinus[g_]:= Mod[-g,N0]
Remove[g];

MakeMinimalGroup[n_]:= Module[{t},
							Zn={};N0=0;Sns={};N1=0;Css={};N2=0;
							t = Commons`EstimatedTime[n];
							If[t>60, 
								If[t>3600, 
									Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/3600]<>" hours."];
									,
									Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/60]<>" minutes."];
								]
								,
								Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t]<>" seconds."]
							];	
							t=TimeUsed[]; 
							N0=n; 
							Zn=Range[0,N0-1]; 
							Sns=Subgroups[];
							N1=Length[Sns];
							t=TimeUsed[]-t;
							If[t > 60, 
								If[t > 3600, 
									Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t/3600]<>" hours."];
									,
									Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t/60]<>" minutes."];
									];
									,
								Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t]<>" seconds."];
							];
							runtime={N[N0],N[t]};
							Save[$HomeDirectory<>"\\makegrouplog.ma",runtime];
						]
Remove[t,n];


(* ::Subsection:: *)
(*Independent instances*)


MakeMinimalGroupInstance[n_]:= Module[{C0,C1},
									C0 = {Zn,N0,N1,Null,Sns,Null};
									MakeMinimalGroup[n];
									C1 = {Zn,N0,N1,Null,Sns,Null};
									{Zn,N0,N1,x1,Sns,x2}=C0;
									Return[C1];
								]
								
InstanceSubgroups[n_]:= MakeMinimalGroupInstance[n][[5]];
Remove[C0,C1,n,x1,x2];


(* ::Subsection:: *)
(*Elementwise*)


ElementOrder[g_]:= Module[{order=1,h=g},
							If[g==0,Return[1],None];
							While[h!=0, h=h\[CirclePlus]g;order++];
							Return[order];]
Remove[order,h,g,Css];

ElementOrders[]:= Association[#->ElementOrder[#]& /@ Zn]


(* ::Subsection:: *)
(*Subsets and related elements*)


Subgroups[]:=Subgroups[N0]
Subgroups[k_]:= ReverseSort[DeleteDuplicates[Union[{0},#1]& /@ AdditiveGroupMinimal`Private`Subcycles[k]]];
Remove[k];

SubgroupsAndGenerator[]:= Union[{0},#]& /@ AdditiveGroupMinimal`Private`SubcycleAndGenerator[]
SubgroupGenerators:= Sort@*Keys@*SubgroupsAndGenerator

Zeros[]:= Table[{Zn[[i]],SuperMinus[Zn[[i]]]},{i,1,N0}]
Remove[i];


(* ::Subsection:: *)
(*Structure and graphical overview*)


CayleyTable[]:= CayleyTable[Zn,CirclePlus]

CayleyTable[G_,op_]:= Grid[Table[op[G[[i]],G[[j]]],
									{i,1,Length[G]},
									{j,1,Length[G]}]/.(0->Item[0,Frame->True])]
Remove[G,op,i,j];							


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];					

Docs[str_]:= With[{package="AdditiveGroup"<>str<>"`*"}, Return[AssociationThread[Names[package]->(Information[#,LongForm->False]& /@ Names[package])]]];

ModularAddition[x1_,x2_]:=If[Dimensions[x1]=={}\[And]Dimensions[x2]=={},
							Return[Mod[x1+x2,N0]],
							If[Dimensions[x1[[1]]]=={}\[And]Dimensions[x2[[1]]]=={},
								Return[Sort@DeleteDuplicates@Flatten@Outer[Mod[#1+#2,N0]&,x1,x2]],
								Return[-1];];]
Remove[x1,x2]

Subcycles::usage = "Int|Int,Int --> {{Int},{Int}}  |  Repeated addition by the start element gives subcycles in cyclic \!\(\*SubscriptBox[\(Z\), \(n\)]\). "<>
													  "Like previous method but does not compute subcycles larger than the second argument given.";
Subcycles[]:= Subcycles[N0]

Subcycles[k_]:= With[{G=Zn},Module[{zero,g,g0,cyclic={},cyclics={},i},
										zero=G[[1]]; 
										For[i=2,i<=N0,i++,
											cyclic={};
											g0=G[[i]];
											g=g0;
											While[g!=zero\[And]Length[cyclic]<=k,
												AppendTo[cyclic,g];
												g=g\[CirclePlus]g0;];
											cyclic=Sort[cyclic];
											If[Length[cyclic]>=k\[Or]MemberQ[cyclics,cyclic],
												None,
												AppendTo[cyclics,cyclic]];
										];
										AppendTo[cyclics,{zero}];
										Return[cyclics];
									]]
Remove[G,zero,g,g0,cyclic,cyclics,i,k]
																
SubcycleAndGenerator::usage = "Int --> <|Int -> {Int}|>   |  Association between generators and their subcycles.";		
SubcycleAndGenerator[]:= With[{Cs=Subcycles[]},
								Return[Association[MapThread[(#1->#2)&,{Min/@Cs,Cs}]]];]
Remove[k,Cs]
								
End[];


Remove[runtime,Css,N2]
EndPackage[];


(* ::Author:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)
