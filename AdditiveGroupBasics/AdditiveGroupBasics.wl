(* ::Package:: *)

(* ::Section:: *)
(*Additive Groups Basics Package*)


PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
BeginPackage["AdditiveGroupBasics`"];
<<AdditiveGroupMinimal`


AdditiveGroupBasicsPackage::usage = "This is the second module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit. The AdditiveGroup package suite."<>
									" It contains the basic functionality. Not quotientgroups and fancy subgroup functionality.";

Print["AdditivegroupBasics`: See Docs[\"Basics\"] for documentation."]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Operators, constants  and primitive mappings*)


(* ::Subsection:: *)
(*Subsets and their properties*)


Zeros::usage = " The zeros of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
SubgroupOrders::usage = " {Int}  |  All orders of the subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
GeneratorsAndOrder::usage = "  <|Int -> Int|>  |  Returns an association of generators and their respective order.";
NonUnitMinimalSubgroupGenerators::usage = " Int --> {Int}  |  The subgroup generators except the trivial ones, 0 and 1.";
ZeroMeetingSubgroups::usage = " {{Int},{Int}}}  |  Pair of subgroups that only have the zero element in common.";


(* ::Subsection:: *)
(*Structure and graphical overview*)


Containment::usage = " Int --> {{Int}}  |  All subgroups that contains a subgroup. The subgroup argument is given as a size-order index. NOTE: the groups in a containment"<>
										  " does not contain each other consecutively. Compare with SubsettingPaths in main package.";
ContainmentNonTrivial::usage = " Int --> {{Int}}  |  Like Containment but with trivial subgroups removed. Se Containment.";

ContainmentHierarchy::usage = " {{{Int}}}  |  Table of all containments of every subgroup.";
ContainmentHierarchyNonTrivial::usage = " {{Int}}  |  Like ContainmentHierarchy but with trivial subgroups removed.";

ContainmentIndexes::usage = " Int --> {Int}  |  Like Containment but returns the size-order indexes of the subgroups only.";
ContainmentIndexesNonTrivial::usage = " Int --> {Int}  |  Like ContainmentIndexes but trivial subgroups are removed.";
ContainmentIndexesHierarchy::usage = " {{Int}} |  Like ContainmentHierarchy but returns the size-order indexes instead of the subgroups themselves.";
ContainmentIndexesHierarchyNonTrivial::usage = " {{Int}}  |  Like ContainmentHieracy but trivial subgroups are removed.";

ContainmentSizes::usage = " {Int}  | The sizes of the subgroups in all containments.";
ContainmentGenerators::usage = " {Int}  |  The generators of subgroups in the containment of a subgroup.";

ContainmentMatrix::usage = " Int -> {{Int}}  |  Matrix of common elements of subgroups in a containment.";

HasseDiagram::usage = " Int|[] --> {{1|0},{{Int}}} |  "<>
								" Adjancy matrix describing the HasseDiagram of the current group's subgroups together with a table of the subgroups"<>
								" concerned. The output is indented to be given to HasseGraph for displaying the HasseDiagram graph of the groups subgroup."<>
								" If an argument integer is given only subgroups of order less than the argument is computed (incomplete to inspect"<>
								" the first parst of the Hassediagram for large groups).";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Operators, constants  and primitive sets and mappings*)


Diamond:=AdditiveGroupBasics`Private`ModularAdditionDirect

Remove[gs1,gs2]


(* ::Subsection:: *)
(*Undependent theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


(* ::Subsubsection:: *)
(*Subsets*)


Zeros[]:= With[{G=Zn},Module[{M={},i},
								For[i=1,i<=N0,i++,
									AppendTo[M,{G[[i]],SuperMinus[G[[i]]]}]];
								Return[M];]]
Remove[G,i]


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


(* ::Subsubsection:: *)
(*Subsets*)


SubgroupOrders[]:= SubgroupOrders[N0];
SubgroupOrders[k_]:= Reverse[Length[#]& /@ Subgroups[k]]

GeneratorsAndOrder[]:= <|KeyValueMap[(#1->Length[#2])&,SubgroupsAndGenerator[]]|>

ZeroMeetingSubgroups[]:= Module[{nonzeros,pairs={},unique,i,j}, 
										nonzeros=Delete[#,1]& /@ Sns;
										nonzeros=nonzeros[[1;;-2]];
										Table[
											If[nonzeros[[i]]\[Intersection]nonzeros[[j]]=={},
												AppendTo[pairs,{nonzeros[[i]],nonzeros[[j]]}]
												,
												None],
											{i,1,N1-1},{j,1,N1-1}
										];	
										unique=DeleteDuplicates[pairs,#1[[1]]==#2[[2]]\[And]#1[[2]]==#2[[1]]&];
										Return[Prepend[#,0]& \[Congruent] unique]]
Remove[nonzeros,pairs,unique,i,j]


(* ::Subsection:: *)
(*Structure and graphical overview*)


(* * * subset sequences * * *)
Containment[k_]:= Sns[[ContainmentIndexes[k]]]
ContainmentNonTrivial[k_]:= Sns[[ContainmentIndexesNonTrivial[k]]];

ContainmentIndexes[k_]:= With[{HS=Sns},
								Module[{H,contained,containing={}},
										contained=HS[[k]];
										Table[H=HS[[i]];
											If[ContainsAll[H,contained],
												AppendTo[containing,i],
												None];
											,{i,1,Length[HS]}];
										Return[containing];]]
Remove[k,HS,H,contained,containing]
									
ContainmentIndexesNonTrivial[k_]:= With[{HS=Sns[[2;;-2]]},
									Module[{H,contained,containing={}},
												contained=HS[[k]];
												Table[H=HS[[i]];
													If[ContainsAll[H,contained],
														AppendTo[containing,i+1],
														None];
													,{i,1,Length[HS]}];
												Return[containing];]]
Remove[k,HS,H,contained,containing]

ContainmentHierarchyNonTrivial[]:=Table[ContainmentNonTrivial[k],{k,1,N1-2,1}]
Remove[k];

ContainmentHierarchy[]:=Table[Containment[i],{i,1,N1}];
ContainmentIndexesHierarchy[]:=Table[ContainmentIndexes[i],{i,1,N1}];
ContainmentIndexesHierarchyNonTrivial[]:=Table[ContainmentIndexesNonTrivial[k],{k,1,N1-2,1}]
Remove[i,k]

ContainmentSizes[]:= Length /@ Containment[#]& /@ Range[N1]
ContainmentGenerators[k_]:= SubgroupGenerators[][[ContainmentIndexes[k]]]
Remove[k]

NonUnitMinimalSubgroupGenerators[]:= Sort[#[[2]]& /@ (Last /@ ContainmentHierarchyNonTrivial[])]
						
HasseDiagram[m_] := With[{HS=Reverse[Subgroups[m]]}, 
							Module[{M,K,i,j,l,connected},
									l=Length[HS];
									M=Table[0,l,l];
									M[[1]]=Table[1,{l}];
									M[[1]][[1]]=0;
									M=PadLeft[#1,l]&/@Table[AdditiveGroupBasics`Private`Weight[HS[[i]],HS[[j]],1&],{i,1,l},{j,i,l}]+M;
									M=AdditiveGroupBasics`Private`OnlyLeadingOne/@M;
									Table[If[Total[M[[All,i]]]==0,M[[1,i]]=1,None],{i,2,l}];
									Return[{M,HS}];
							]]
Remove[m,HS,M,K,i,j,l,connected]
							
HasseDiagram[] := HasseDiagram[N0]


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];

ModularAdditionDirect[x1_,x2_]:=If[Dimensions[x1]=={}\[And]Dimensions[x2]=={},
											Return[Mod[x1+x2,N0]]
											,
											If[Dimensions[x1[[1]]]=={}\[And]Dimensions[x2[[1]]]=={},
												Return[Outer[Mod[#1+#2,N0]&,x1,x2]]
												,
												Return[-1];];]
Remove[x1,x2]
																								
(* naming - edge weight - also more generic *)
Weight[H1_,H2_,w_] := With[{},
							If[ContainsAll[H2,H1]\[And]H1!=H2,
								Return[w[H1,H2]]
								,
								Return[0]
							];]
Remove[H1,H2,w]

(* to complicated as well *) 
OnlyLeadingOne[xs_] := Module[{all, leading, trailing},
									all = SequenceSplit[xs,{0, 1}:>{0, 1}];
									If[Length[all] >= 3,
										leading = Flatten[all[[1 ;; 2]]];
										trailing=ConstantArray[0,Length[xs]-Length[leading]];
										,
										None];
									If[Length[all]==2,
										If[all[[1]]=={0, 1},
											leading=all[[1]];trailing=ConstantArray[0,Length[all[[2]]]];
											,
											leading=ConstantArray[0,Length[all[[1]]]];trailing=all[[2]];]
										,
										None];
										If[Length[all]==1,
											leading={};trailing=ConstantArray[0,Length[all[[1]]]];
											,
											None;];
										Return[Join[leading,trailing]];
								]
Remove[k,all,leading,trailing,xs]
End[];


(* ::Subsection:: *)
(*Seldom used and helpers*)


ContainmentMatrix[k_]:= With[{all=Containment[k]},Module[{n=Length[all]},
							Return[Table[If[ContainsAll[all[[i]],all[[j]]]\[Or]ContainsAll[all[[j]],all[[i]]],i,0],{i,1,n},{j,1,n}]];]]
Remove[k,n,all,i,j]							


(* ::Section::Closed:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];

