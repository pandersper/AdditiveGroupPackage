(* ::Package:: *)

(* ::Section:: *)
(*Additive Groups Minimal Package*)


PrependTo[$ContextPath,"Commons`"];
BeginPackage["AdditiveGroupMinimal`"];
<< Commons`


AdditiveGroupMinimalPackage::usage = "This is the first module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit. The AdditiveGroup package suite.\
									It contains the minimal functionality needed for investigating \!\(\*SubscriptBox[\(Z\), \(n\)]\) and it's subgroups.";


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Operators, constants  and primitive mappings*)


N0::usage = " Int  |  Modulus of the operators and size of Zn."; 
N1::usage = " Int  |  Number of subgroups in Zn. Modulus of the total quotient group operator."; 

Zn::usage = " {Int}  |  List of the elements in the additive group modulo N0.";
Sns::usage = " {{Int}}  |  Subgroups of the group currently in use.";

CirclePlus::usage = " Int,Int --> Int  |  Group operator (addition).";
SuperMinus::usage = " Int --> Int  |  The inverse of g modulo N0.";

MakeMinimalGroup::usage = " Int --> Null  |  Computes an additive group of order n and also precomputes it's subgroups. All variables it computes is global to current context.";
MakeMinimalGroupInstance::usage = " Int --> {{Int},Int,Int,Int,{Int},{Int}}  |  Computes an additive group instance containing also it's subgroups, though not it's quotient groups.\
																			    It temporarily alters but restores the current context.";
InstanceSubgroups::usage = " Int -->  {{Int}}  |  Returns the subgroups of the additive group of a given order. It temporarily alters but restores the current context group.";


(* ::Subsection:: *)
(*Properties of single elements in the group*)


ElementOrder::usage = " Int,Int--> Int  | The order of an element in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
ElementOrders::usage = " Int -->  Int  |  All the orders of the elements in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";


(* ::Subsection:: *)
(*Subsets and their properties*)


Subgroups::usage = " {{Int}}  |  The subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
SubgroupsAndGenerator::usage = " <|Int -> {Int}|>   |  The subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) mapped to by their generator as key.";
SubgroupGenerators::usage = " {Int}  |  The generators of all subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) given in subgroup size order.";


(* ::Subsection:: *)
(*Structure and graphical overview*)


CayleyTable::usage = " Int --> SquareMatrix[Int]  |  Multiplication table of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
CayleyTable::usage = " {Int},(Int,Int -> Int) --> Grid[Int]  |  Multiplication table of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Operators, constants  and primitive sets and mappings*)


CirclePlus:= AdditiveGroupMinimal`Private`ModularAddition
SuperMinus[g_]:= Mod[-g,N0]

Remove[g];


(* ::Subsection:: *)
(*Undependent theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


MakeMinimalGroup[n_]:= Module[{t,runtime},
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
Remove[t,runtime,n];
							
ElementOrder[g_]:= Module[{order=1,h=g},
							If[g==0,Return[1],None];
							While[h!=0, h=h\[CirclePlus]g;order++];
							Return[order];
						]
Remove[order,h,g];


(* ::Subsubsection:: *)
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
(*Theory*)


(* ::Subsubsection::Closed:: *)
(*Elementwise*)


ElementOrders[]:= Association[#->ElementOrder[#]& /@ Zn]


(* ::Subsubsection:: *)
(*Subsets*)


(* * * subgroups * * *)
SubgroupsAndGenerator[]:= Union[{0},#]& /@ AdditiveGroupMinimal`Private`SubcycleAndGenerator[]

SubgroupGenerators:= Sort@*Keys@*SubgroupsAndGenerator

Subgroups[]:=Subgroups[N0]

(* * * lossy faster versions * * *)
Subgroups[k_]:= ReverseSort[DeleteDuplicates[Union[{0},#1]& /@ AdditiveGroupMinimal`Private`Subcycles[k]]];
Remove[k];


(* ::Subsection::Closed:: *)
(*Structure and graphical overview*)


CayleyTable[]:= Module[{G=Zn,M={},row,i},
							For[i=1,i<=N0,i++,
								row=(G[[i]]\[CirclePlus]#1)& /@ G;
								AppendTo[M,row]
							];
							Return[Grid[M/.(0->Item[0,Frame->True])]];
						]
Remove[G,M,row,i];

CayleyTable[G_,op_]:= Module[{M={},row,i},
								For[i=1,i<=Length[G],i++,
									row=op[G[[i]],#1]&/@G;
									If[!AtomQ[row[[1]]],row=Sort/@row,None];
									AppendTo[M,row]
								];
								Return[Grid[M/.G[[1]]->Item[G[[1]],Frame->True]]];
							]
Remove[G,op,M,row,i];							


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];					
					
ModularAddition[x1_,x2_]:=If[Dimensions[x1]=={}\[And]Dimensions[x2]=={},
							Return[Mod[x1+x2,N0]],
							If[Dimensions[x1[[1]]]=={}\[And]Dimensions[x2[[1]]]=={},
								Return[Sort@DeleteDuplicates@Flatten@Outer[Mod[#1+#2,N0]&,x1,x2]],
								Return[-1];];]
Remove[x1,x2]

Subcycles::usage = "Int|Int,Int --> {{Int},{Int}}  |  Repeated addition by the start element gives subcycles in cyclic \!\(\*SubscriptBox[\(Z\), \(n\)]\). \
													  Like previous method but does not compute subcycles larger than the second argument given.";
Subcycles[]:= Subcycles[N0]

SubcycleAndGenerator::usage = "Int --> <|Int -> {Int}|>   |  Association between generators and their subcycles.";		
Subcycles[k_]:= With[{G=Zn},Module[{zero,g,g0,cyclic={},cyclics,i},
										zero=G[[1]]; 
										cyclics ={{zero}};
										For[i=2,i<=N0,i++,
											cyclic={};
											g0=G[[i]];
											g=g0;
											While[g!=zero\[And]Length[cyclic]<=k,AppendTo[cyclic,g];g=g\[CirclePlus]g0;];
											cyclic=Sort[cyclic];
											If[Length[cyclic]>=k\[Or]MemberQ[cyclics,cyclic],None,AppendTo[cyclics,cyclic]];
										];
										Return[cyclics];
									]]
Remove[G,zero,g,g0,cyclic,cyclics,i,k]
																
SubcycleAndGenerator[]:= With[{Cs=Subcycles[]},
								Return[Association[MapThread[(#1->#2)&,{Min/@Cs,Cs}]]];]
Remove[k,Cs]
								
End[];


(* ::Section:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
