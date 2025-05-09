(* ::Package::"Tags"-><|"UppercaseParameter" -> <|Enabled -> False|>|>:: *)

(* ::Section:: *)
(*The Additive Group  Main Package*)


PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
PrependTo[$ContextPath,"AdditiveGroupBasics`"];
BeginPackage["AdditiveGroup`"];
<< AdditiveGroupBasics`


AdditiveGroupPackage::usage = "This is is the third and main module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit. The AdditiveGroup package suite."<>
							  " It contains the the more fancy functionality related to investiganting \!\(\*SubscriptBox[\(Z\), \(n\)]\) and it's subgroups"<>
							  " without transcending to quotient groups.";
							  
Print["Additivegroup`: See Docs[\"\"] for documentation."]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Operators, constants  and primitive mappings*)


Permutation::usage = " Int --> {Int}  |  Returns the Mathematica permutation representation corresponding to the permutation "<>
										"of \!\(\*SubscriptBox[\(Z\), \(n\)]\)that an element gives rise to.";
InnerAut::usage = " Int,Int --> Int  |  The inner g-automorpism of the element x in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

CyclesList::usage = " {Cycles}  |  List of permutations on cycles-form corresponding to the permutation group isomorphic to \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
CyclesMap::usage = " <|Int->Cycles[]|>  |  Mapping of elements of \!\(\*SubscriptBox[\(Z\), \(n\)]\) to their corresponding permutation in the isomorphic permutation group.";


(* ::Subsection:: *)
(*Properties of single elements in the group*)


(* ::Subsection:: *)
(*Properties of whole group*)


(* ::Subsection:: *)
(*Subsets and their properties*)


SubgroupOrder::usage = " Int --> Int  |  The order of a subgroup. Subgroups are indexed after their sizes.";
GeneratorSpan::usage = " Int --> {Int}  |  Returns the set of elements that the generator generates, spans.";

ContainmentPaths::usage = " {{{Int}}}  |  Returns all possible non-trivial smaller containments one can form from a containment."
ContainmentPathsIndexes::usage = " {{Int}}  | Returns all possible non-trivial smaller containments one can form from a containment indexed by subgroup's size.";
ContainmentPathsContainmentIndexes::usage = " {{Int}}  | Returns all possible non-trivial smaller containments one can form from a containment indexed by the containment's subgroups.";

NonEntanglingPaths::usage = "  {{{Int}}}  |   Returns the containment paths (see ContainmentPaths) where the inclusion operator works every step along the path.";
NonEntanglingPathsIndexed::usage = "  {{Int}}  |  Like NonEntanglingPaths but with subgroup indexes. See NonEntanglingPaths.";

ContainmentMaximum::usage = " {{Int}} |  The containment paths of the conatinment target group that has most paths to it.";

PartitioningSubgroups::usage = "EXIST NOT";
ZeroMeetFactorisation::usage = "  <|{Int,Int}->Int|>  |  Association of indexes of zero-meeting subgroups to their product. Se ZeroMeetings in basic package.";

SubgroupProducts::usage = " {{Int}}  |  Returns the matrix i.e the multiplication table of all products of two subgroups.";
SubgroupProductBoundary::usage = "{{Int}} |  Returns the differing set between the union of two subgroups and the product of them,";
SubgroupProductBoundaries::usage = " {{{Int}}}  |  Returns all subgroup product boundaries between pairs of subgroups in \!\(\*SubscriptBox[\(Z\), \(n\)]\)";
SubgroupIntersections::usage = " {{Int}}  |  Returns the matrix i.e the multiplication table of all intersections of two subgroups.";
SubgroupSymmetricDiffs::usage = " {{Int}}  |  Returns the matrix i.e the mutiplication table of all symmetric differences of two subgroups.";
													 
ContainmentAndComplement::usage = " Int --> {{Int},{Int}}  |  A subgroups containment together with it\.b4s complement.";
ContainmentAndComplementNonTrivial::usage = " Int --> {{Int},{Int}} |  Like ContainmentAndComplement but with trivial subgroups removed.";

ContainmentIndexesFromGeneratorNonTrivial::usage = "  {{Int}}  |  Gives the containment corresponding to a generator but with the trivial subgroups removed.";


(* ::Subsection:: *)
(*Derived sets and corresponding properites*)


(* ::Subsection:: *)
(*Isomorphy*)


AdditivePermutationGroup::usage = " Int --> PermutationGroup  |  The permutation group isomorphic to \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
AdditivePermutationSubgroups::usage = " Int --> {PermutationGroup}  |  The subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) as isomorphic permutation groups.";


(* ::Subsection:: *)
(*Structure and graphical overview*)


WholeGroupProducts::usage = " {{1|0}}  |  Returns a matrix that shows which pair of subgroups have the product \!\(\*SubscriptBox[\(Z\), \(n\)]\). ";
SubgroupProducing::usage = " {{1|0}}  |  Returns a matrix that shows which pair of subgroups have a certain subgroup as the product. ";
AggregatedSubgroupProducing::usage = " {{{1|0}}}  |  Returns the consecutive matrix sums of subgroup producing pairs, producing subgroup number "<>
													 "one up to the one with index given as argument.";													 

DivergerMatrix::usage = " {{1|0}}  |  Triangular matrix with ones for pairs of subgroups where the smaller subgroup is not fully contained "<>
									  "in the larger and zeroes everywhere else.";
DivergerMatrixBool::usage = " {{True|False}}  |  Like DivergerMatrix but with bool value instead of 1|0.";
DivergerPairs::usage = " {{Int,Int}}  |  The pair of subgroups derived from DivergerMatrix.";

HasseWeb::usage = " {{1|0}}  |  The Hassediagram adjancy matrix (see HasseDiagram in basics package) but also with all other subset paths "<>
								"between subgroups represented as a one. Very many for larger groups and difficult to visualise.";

HasseGraph::usage = " Int --> Graph  |  The the concrete graph visualisation of a Hasse-diagram. See HasseDiagram. ";
HasseGraphEdges::usage = " Int --> Graph  |  Like HasseGraph but with the extension of cyclic group on the edges.";
HasseGraphEdges::usage = " Int --> Graph  |  Like HasseGraphEdges[n] but subgroups with order larger than the second argument is discarded. "<>
											 "The beginning of the diagram for groups that are not feasible to compute.";

HasseWebGraphEdges::usage = " Graph  |  Like HasseGraph but based on a HasseWeb instead of a HasseDiagram. See HasseDiagram in basics package "<>
										"and HasseWeb int this. ";
										
HasseBoth::usage = " Graph  |  The HasseWeb graph with contained HasseGraph higlighted.";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Operators, constants  and primitive sets and mappings*)


(* ::Subsection::Closed:: *)
(*Undependent theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


Permutation[i_]:= With[{G=Zn},Module[{mappings={},j},
									If[G[[i]]==0, 
										Return[(#+1)& /@ G]
										,
										Table[AppendTo[mappings,(G[[j]]\[CirclePlus]G[[i]])+1],{j,N0}];		
										Return[mappings];
									];
						]]
Remove[i,j,G,mappings]

InnerAut[x_,g_]:= g\[CirclePlus]x\[CirclePlus]SuperMinus[g]
Remove[g,x]

CyclesMap[]:= Module[{G=Zn,cycles=<||>,i},
							Table[AppendTo[cycles,G[[i]]->PermutationCycles[Permutation[i]]],{i,N0}];
							Return[cycles];						
						]
Remove[G,cycles,i]

CyclesList[]:= Module[{cycles={},i},
							Table[AppendTo[cycles,PermutationCycles[Permutation[i]]],{i,N0}];
							Return[cycles];						
						]
Remove[cycles,i]


(* ::Subsubsection:: *)
(*Subsets*)


(* ::Subsubsection:: *)
(*Independent instances*)


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


(* ::Subsubsection:: *)
(*Subsets*)


SubgroupOrder[k_]:= SubgroupOrders[N0][[-k]]
Remove[k]

GeneratorSpan[g_]:=AdditiveGroup`Private`TotalSpan[{g}]
Remove[g]

ContainmentPathsIndexes[k_]:=With[{all=Containment[k]},
								Module[{n=Length[all],containmentmatrix,itrees,ipaths,nonzeros},
									containmentmatrix=Table[
															If[ContainsAll[all[[i]],all[[j]]]\[Or]ContainsAll[all[[j]],all[[i]]],i,0],
														{i,1,n},{j,1,n}];
									itrees=Select[containmentmatrix\[Transpose],MemberQ[0]];
									ipaths=MaximalBy[itrees,Length[Select[#,(#==0)&]]&];
									nonzeros=(Select[#,(#!=0)&])& /@ ipaths;
									If[ipaths!={},
										Return[Flatten /@ ((all[[#]]& /@ nonzeros)/.PositionIndex[Sns])],
										Return[{Flatten[all/.PositionIndex[Sns]]}]];
								]];
Remove[k,all,n,containmentmatrix,itrees,ipaths,nonzeros]

ContainmentPathsContainmentIndexes[k_]:=With[{all=Containment[k]},
								Module[{n=Length[all],containmentmatrix,itrees,ipaths,nonzeros},
									containmentmatrix=Table[
															If[ContainsAll[all[[i]],all[[j]]]\[Or]ContainsAll[all[[j]],all[[i]]],i,0],
														{i,1,n},{j,1,n}];
									itrees=Select[containmentmatrix\[Transpose],MemberQ[0]];
									ipaths=MaximalBy[itrees,Length[Select[#,(#==0)&]]&];
									nonzeros=(Select[#,(#!=0)&])& /@ ipaths;
									If[ipaths!={},
										Return[nonzeros],
										Return[{Flatten[all/.PositionIndex[all]]}]];
								]];
Remove[k,all,n,containmentmatrix,itrees,ipaths,nonzeros]

ContainmentPaths[k_]:=With[{all=Containment[k],ipaths=ContainmentPathsContainmentIndexes[k]},
						Module[{nonzeros=(Select[#,(#!=0)&])& /@ ipaths},
								Return[all[[#]]& /@ nonzeros];]]
Remove[k,all,n,ipaths,nonzeros]

NonEntanglingPaths[k_]:= Join[#,{Sns[[k]]}]& /@ Select[((ContainmentPaths[k]\[Transpose])[[1;;-2]])\[Transpose],((Fold[Intersection[#2,#1]&,Sns[[1]],#])!={0})&]

NonEntanglingPathsIndexed[k_]:= Flatten /@ (NonEntanglingPaths[k]/.PositionIndex[Sns]);
Remove[k]

ContainmentMaximum[]:= MaximalBy[NonEntanglingPathsIndexed[#]& /@ Range[N1],Length];

ZeroMeetFactorisation[]:= Module[{zeromeets=ZeroMeetingSubgroups[],partials}, 
							partials=Table[zeromeets[[i,1]]\[CirclePlus]zeromeets[[i,2]],{i,1,Length[zeromeets]}];
							Return[AssociationThread[Flatten[(#/.PositionIndex[Sns])&/@zeromeets,{1,3}],First /@ (partials/.PositionIndex[Sns])]]]
Remove[zeromeets,partials,i]

SubgroupProducts[]:=Table[Sns[[i]]\[CirclePlus]Sns[[j]],{i,1,N1},{j,1,N1}]
Remove[i,j]

SubgroupProductBoundary[G1_,G2_]:=Complement[(G1\[CirclePlus]G2),G1\[Union]G2]
Remove[G1,G2]

SubgroupProductBoundaries[]:=Table[SubgroupProductBoundary[Sns[[i]],Sns[[j]]],{i,1,N1},{j,1,N1}]
Remove[i,j]

SubgroupIntersections[]:=Table[Sns[[i]]\[Intersection]Sns[[j]],{i,1,N1},{j,1,N1}]
Remove[i,j]

SubgroupSymmetricDiffs[]:=Table[SymmetricDifference[Sns[[i]],Sns[[j]]],{i,1,N1},{j,1,N1}]
Remove[i,j]

ContainmentAndComplement[k_]:=With[{in=Containment[k],all=Sns},Return[{in,Complement[all,in]}]];
Remove[k,in,all]
ContainmentAndComplementNonTrivial[k_]:=With[{in=ContainmentNonTrivial[k],all=Sns[[2;;-2]]},Return[{in,Complement[all,in]}]]
Remove[k,in,all]

ContainmentIndexesFromGeneratorNonTrivial[g_]:= ContainmentIndexes[Flatten[Position[Sns,GeneratorSpan[g]]][[1]]]
Remove[g]


(* ::Subsection:: *)
(*Isomorphy*)


AdditivePermutationGroup[]:= PermutationGroup[CyclesList[]]
AdditivePermutationSubgroups[]:= PermutationGroup /@ (Sns/.CyclesMap[])


(* ::Subsection:: *)
(*Structure and graphical overview*)


(* ::Code::Initialization::"Tags"-><|"UppercaseParameter" -> <||>|>:: *)
WholeGroupProducts[]:= If[(#==Sns[[1]]),1,0]& \[Congruent] SubgroupProducts[] 
SubgroupProducing[k_]:= If[(#==Sns[[k]]),1,0]& \[Congruent] SubgroupProducts[] 
AggregatedSubgroupProducing[]:=Table[Fold[Plus,0,Table[SubgroupProducing[i],{i,1,N1-j}]],{j,N1-1,0,-1}]
Remove[i,j,k]

DivergerMatrix[]:=Table[If[ContainsAll[Sns[[i]],Sns[[j]]]\[Or]i>j,0,1],{i,1,N1},{j,1,N1}];
Remove[i,j]
DivergerMatrixBool[]:=Table[ContainsAll[Sns[[i]],Sns[[j]]]\[Or]i>j,{i,1,N1},{j,1,N1}];
Remove[i,j]
DivergerPairs[]:=Flatten[Table[If[ContainsAll[Sns[[i]],Sns[[j]]]\[Or]i>j,Nothing,{Sns[[i]],Sns[[j]]}],{i,1,N1},{j,1,N1}],1]
Remove[i,j]

HasseWeb[] := HasseWeb[N0]
HasseWeb[m_] := With[{HS=Reverse[Subgroups[m]]}, 
							Module[{M,K,i,j,connected},
									l=Length[HS];
									M=Table[0,l,l];
									(*M[[1]]=Table[1,{l}];*)
									M[[1]][[1]]=0;
									M=PadLeft[#1,l]&/@Table[AdditiveGroupBasics`Private`Weight[HS[[i]],HS[[j]],1&],{i,1,l},{j,i,l}]+M;
									
									Table[If[Total[M[[All,i]]]==0,M[[1,i]]=1,None],{i,2,l}];
									Return[{M,HS}];]]
Remove[m,HS,M,K,i,j,connected,l]
							
HasseGraph[]:= HasseGraph[N0]
HasseGraph[m_]:=With[{MHS=HasseDiagram[m],orders=SubgroupOrders[m]}, 
								Module[{M,HS,vtxpnl,i},
									
									vtxpnl[v_]:= Panel[v,FrameMargins->0,Background->Lighter[Green,0.9]];
		
									M=MHS[[1]];HS=MHS[[2]];
									
									Return[AdjacencyGraph[M,
															VertexLabels->
																Table[i->
																Placed["\!\( o \_ "<>ToString[i]<>" \)="<>ToString[orders[[i]]],
																	Center,vtxpnl], 
																	{i,Length[orders]}],
															VertexShapeFunction->"Rectangle"
															(** ,EdgeShapeFunction->{{"Arrow","ArrowSize"->.03}}**)
														]
												];]]
Remove[m,MHS,orders,M,HS,vtxpnl,i,v]

HasseGraphEdges[]:= HasseGraphEdges[N0]
HasseGraphEdges[m_]:= With[{MHS=HasseDiagram[m],orders=SubgroupOrders[m]},
											Module[{
													 HS,expansions,
													 M,A,
													 X,VS,ES,
													 edge,to,from,expansion,
													 vtxpnl,vtxlbltbl,edglbltbl,title},
													
													vtxpnl[v_]:= Panel[v,FrameMargins->0,Background->Lighter[Green,0.9]];
													
													M=MHS[[1]];
													HS=MHS[[2]];
													
													A=AdjacencyGraph[M,DirectedEdges->True];
													ES=EdgeList[A];
													VS=VertexList[A];
													
													
													expansions=Table[
																		edge=ES[[i]];to=edge[[2]];from=edge[[1]];
																		expansion=Complement[HS[[to]],HS[[from]]];
																		edge->ToString[expansion],
																		{i,Length[ES]}];
																			
													vtxlbltbl=Table[i->Placed[
																			ToString[i]<>": "<>"\!\( o \_ " <>ToString[VS[[i]]]<>" \)="<>ToString[orders[[i]]],
																			Center,vtxpnl],
																	{i,Length[orders]}];	
																											
													edglbltbl=Table[
																	expansions[[i]][[1]]->"\[Union]<"<>ToString[ToExpression[expansions[[i]][[2]]][[1]]]<>">",
																	{i,Length[ES]}];
																						
													title=StringJoin[("\!\("<>ToString[#[[1]] ]<>"\^"<>ToString[#[[2]]]<>"\)")&/@FactorInteger[N0]];
													
													X=Graph[ES,
																EdgeWeight->expansions,
																VertexLabels->vtxlbltbl,VertexShapeFunction->"Rectangle",
																EdgeLabels->edglbltbl,
																PlotLabel->title];
													Return[X];]]
Remove[m,MHS,orders,HS,expansions,expansion,M,A,X,VS,ES,edge,to,from,vtxpnl,vtxlbltbl,edglbltbl,title,i]

HasseWebGraphEdges[]:= HasseWebGraphEdges[N0]
HasseWebGraphEdges[m_]:= With[{MHS=HasseWeb[m],orders=SubgroupOrders[m]},
											Module[{
													 HS,expansions,
													 M,A,
													 X,VS,ES,
													 edge,to,from,expansion,
													 vtxpnl,vtxlbltbl,edglbltbl,title},
													
													vtxpnl[v_]:= Panel[v,FrameMargins->0,Background->Lighter[Green,0.9]];
													
													M=MHS[[1]];
													HS=MHS[[2]];
													
													A=AdjacencyGraph[M,DirectedEdges->True];
													ES=EdgeList[A];
													VS=VertexList[A];
													
													expansions=Table[
																		edge=ES[[i]];to=edge[[2]];from=edge[[1]];
																		expansion=Complement[HS[[to]],HS[[from]]];
																		edge->ToString[expansion],
																		{i,Length[ES]}];
																			
													vtxlbltbl=Table[i->Placed[
																			"\!\( o \_ " <>ToString[VS[[i]]]<>" \)="<>ToString[orders[[i]]],
																			Center,vtxpnl],
																	{i,Length[orders]}];	
																											
													edglbltbl=Table[
																	expansions[[i]][[1]]->"\[Union]<"<>ToString[ToExpression[expansions[[i]][[2]]][[1]]]<>">",
																	{i,Length[ES]}];
																						
													title=StringJoin[("\!\("<>ToString[#[[1]] ]<>"\^"<>ToString[#[[2]]]<>"\)")&/@FactorInteger[N0]];
													
													X=Graph[ES,
																EdgeWeight->expansions,
																VertexLabels->vtxlbltbl,VertexShapeFunction->"Rectangle",
																EdgeLabels->edglbltbl,
																PlotLabel->title];
													
													Return[X];]]
Remove[m,MHS,orders,HS,expansions,expansion,M,A,X,VS,ES,edge,to,from,vtxpnl,v,vtxlbltbl,edglbltbl,title,i]

HasseBoth[]:= Module[{nonplanar = HasseWebGraphEdges[], planar = HasseGraphEdges[]},
				Return[HighlightGraph[nonplanar,{planar,GraphDifference[nonplanar,planar]}]];]



Begin["`Private`"];

TotalSpan[list_]:= Module[{span,length,newlength,absid,changing=True},
								span=list;
								length=Length[span];
								While[changing,
									span=span\[Union]Flatten[Outer[#1\[CirclePlus]#2&,span,span]];
									newlength = Length[span];
									changing = length != newlength;
									length = newlength;
								];
								Return[span];
							]						
Remove[list,k,span,length,changing];

End[];


(* ::Section:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
