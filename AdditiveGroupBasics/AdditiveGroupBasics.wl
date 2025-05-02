(* ::Package:: *)

(* ::Section::Closed:: *)
(*Additive Groups Basics Package*)


PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
BeginPackage["AdditiveGroupBasics`"];
<<AdditiveGroupMinimal`


AdditiveGroupBasicsPackage::usage = "This is the second module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit. The AdditiveGroup package suite."<>
									" It contains the basic functionality. Not quotientgroups and fancy subgroup functionality.";

Print["AdditivegroupBasics`: See Docs[\"Basics\"] for documentation."]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*Operators, constants  and primitive sets and mappings*)


Diamond::usage = " Modular addition on cosets with cosets preserved as elements, not representatives.";


(* ::Subsection::Closed:: *)
(*Subsets*)


SubgroupOrders::usage = " {Int}  |  All orders of the subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

GeneratorsAndOrder::usage = "  <|Int -> Int|>  |  Returns an association of generators and their respective order.";

ZeroMeetingSubgroups::usage = " {{Int},{Int}}}  |  Pair of subgroups that only have the zero element in common.";

SubgroupProducts::usage = " {{Int}}  |  Returns the matrix i.e the multiplication table of all products of two subgroups.";

SubgroupIntersections::usage = " {{Int}}  |  Returns the matrix i.e the multiplication table of all intersections of two subgroups.";


(* ::Subsection::Closed:: *)
(*Structure and graphical overview*)


Containment::usage = " Int --> {{Int}}  |  All subgroups that contains a subgroup. The subgroup argument is given as a size-order index. NOTE: the groups in a containment"<>
										  " does not contain each other consecutively. Compare with SubsettingPaths in main package.";

ContainmentIndexes::usage = " Int --> {Int}  |  Like Containment but returns the size-order indexes of the subgroups only.";

ContainmentGenerators::usage = " {Int}  |  The generators of subgroups in the containment of a subgroup.";

HasseDiagram::usage = " Int|[] --> {{1|0},{{Int}}} |  "<>
								" Adjancy matrix describing the HasseDiagram of the current group's subgroups together with a table of the subgroups"<>
								" concerned. The output is indented to be given to HasseGraph for displaying the HasseDiagram graph of the groups subgroup."<>
								" If an argument integer is given only subgroups of order less than the argument is computed (incomplete to inspect"<>
								" the first parst of the Hassediagram for large groups).";
								
ShortestHassePath::usage = " {<|Int->{Int}|>}  |  The shortest path in the Hasse diagram of \!\(\*SubscriptBox[\(Z\), \(n\)]\). Association of rules of group indexes to groups. See HasseGraphEdges.";
LongestHassePath::usage = " {Int}  |  The longest path in the Hasse diagram of \!\(\*SubscriptBox[\(Z\), \(n\)]\). Association of rules of group indexes to groups. See HasseGraphEdges.";

ShortestHasseGeneratorPath::usage = " {<|Int->Int|>}  |  The shortest paths of generator extensions in the Hasse diagram of \!\(\*SubscriptBox[\(Z\), \(n\)]\). List of associations between group index and generator. See HasseGraphEdges.";
LongestHasseGeneratorPath::usage = " {<|Int->Int|>} |  The longest paths of generator extensions in the Hasse diagram of \!\(\*SubscriptBox[\(Z\), \(n\)]\).  List of associations between group index and generator. See HasseGraphEdges.";


(* ::Section:: *)
(*Code*)


(* ::Subsection::Closed:: *)
(*Operators, constants  and primitive sets and mappings*)


Diamond:=AdditiveGroupBasics`Private`ModularAdditionDirect
Remove[gs1,gs2]


(* ::Subsection::Closed:: *)
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

SubgroupProducts[]:=Table[Sns[[i]]\[CirclePlus]Sns[[j]],{i,1,N1},{j,1,N1}]
Remove[i,j]

SubgroupIntersections[]:=Table[Sns[[i]]\[Intersection]Sns[[j]],{i,1,N1},{j,1,N1}]
Remove[i,j]


(* ::Subsection::Closed:: *)
(*Structure and graphical overview*)


(* * * subset sequences * * *)
Containment[k_]:= Sns[[ContainmentIndexes[k]]]

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

ContainmentGenerators[k_]:= SubgroupGenerators[][[ContainmentIndexes[k]]]
Remove[k]
															
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

ShortestHassePath[]:= With[{MHS=HasseDiagram[]},Module[{M,HS,A,G,VS,ES,lmin,p,ps},
			
									M=MHS[[1]]; HS=MHS[[2]];
									
									A = AdjacencyGraph[M,DirectedEdges->True];
									ES = EdgeList[A]; VS = VertexList[A];
																		
									G=Graph[ES, EdgeWeight -> {_ -> 1}];									
									
									p=FindShortestPath[G,VS[[1]],VS[[-1]]];
									lmin = Length[p];
									ps = FindPath[G,VS[[1]],VS[[-1]],{lmin-1},Infinity];
									
									Return[AssociationMap[HS[[#]]&,##]&/@ps];
								]]

ShortestHasseGeneratorPath[]:= ShortestHassePath[]/.KeyValueMap[(#2->#1)&,SubgroupsAndGenerator[]]

LongestHassePath[]:= With[{MHS=HasseDiagram[]}, Module[{M,HS,A,G,VS,ES,lmin,p,ps},

									M=MHS[[1]]; HS=MHS[[2]];
								
									A = AdjacencyGraph[M,DirectedEdges->True];
									ES = EdgeList[A]; VS = VertexList[A];
									
									G = Graph[ES, EdgeWeight -> {_ -> 1}];									
								
									p=FindShortestPath[G,VS[[1]],VS[[-1]]];
									lmin = Length[p];
									ps = FindPath[G,VS[[1]],VS[[-1]],{lmin,Length[VS]},Infinity];
									ps = Select[ps,(Length[#]== Max[Length /@ ps])&];
									
									If[ps=={},
										Return[ShortestFactorisations[n]],
										Return[AssociationMap[HS[[#]]&,##]&/@ps]];
								]]

LongestHasseGeneratorPath[]:= LongestHassePath[]/.KeyValueMap[(#2->#1)&,SubgroupsAndGenerator[]]
Remove[M,HS,A,G,VS,ES,lmin,p,ps];



(* ::Subsection::Closed:: *)
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


(* ::Section::Closed:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
