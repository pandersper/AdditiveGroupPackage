(* ::Package:: *)

(* ::Section:: *)
(*The Additive Group  Theorems Package*)


PrependTo[$ContextPath,"Commons`"];
PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
PrependTo[$ContextPath,"AdditiveGroupBasics`"];
PrependTo[$ContextPath,"AdditiveGroup`"];
PrependTo[$ContextPath,"AdditiveGroupQuotients`"];

BeginPackage["AdditiveGroupTheorems`"];
<< AdditiveGroupQuotients`
Needs["Commons`"]


AdditiveGroupTheoremsPackage::usage = "This is is the fifth module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suite, the AdditiveGroup package suite. It implements the \
										correspondence theorem and the three isomorphism theorems with some adjoining functionality.";


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Correspondence theorem*)


(* * * group side of correspondence theorem * * *)
ContainmentCorrespondenceGroups::usage = " {{Int}} | Groups corresponding to the subgroups in the containment of the given group. See Containment.";
ContainmentCorrespondenceOrder::usage = " Int | The order of the largest group in a containment correspondance. See ContainmentGroups.";
AllContainmentCorrespondencies::usage = " {{Int}} | Sets of sets of corresponding groups (containment correspondancies) for all containments in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
AllContainmentsCorrespondenceOrders::usage = "  {Int} | All containmnent correspondancies's orders.";

(* * * quotient group side of correspondence theorem * * *)
QuotientCorrespondenceGroups::usage = " {{Int}} | Groups corresponding to the quotient groups derived from dividing all larger groups with the one given.";
QuotientCorrespondenceOrder::usage = " Int | The order of the largest group in quotient group correspondance. Se QuotientGroups. ";
AllQuotientCorrespondencies::usage = " {{Int}} | Sets of sets of corresponding groups (quotient group correspondancies) for all quotientgroups.";
AllQuotientCorrespondenceOrders::usage = " {Int} | All quotient group correspondancies's orders.";

CorrespondenceMap::usage = " <|{Int} -> {Int}|> |   The correspondance map for a given subgroup. Correspondance between subgroups containing the subgroup "<>
													"and quotient groups derived by letting larger groups taking quotients with it.";			
CorrespondenceMapIndexed::usage = " <| Int -> Int |>  |  Like CorrespondenceMap but the group in the domain and range are given by indexes. The indexes are "<>
															"given by their size ordering in \!\(\*SubscriptBox[\(Z\), \(n\(.\)\)]\)";
															
UniqueCorrespondenceMap::usage = " <| Int -> Int |>  |  Like CorrespondenceMap but all non-unique occurencies are removed.";
UniqueCorrespondenceMapIndexed::usage = " <| Int -> Int |>  |  Like CorrespondanceMapIndexed but all non-unique occurencies are removed.";

CorrespondenceMethodDictionary::usage = " <| Int -> {Int}  | The same result as CorrespondenceMapSubgroupsIndexed except it is constructed on the fly "<>
															 "by looking up subgroups in \!\(\*SubscriptBox[\(Z\), \(n\)]\). (Algorithm 1 to compute correspondence instead of finding it.)";
CorrespondenceMethodGenerators::usage = " <| Int -> {Int}  | The same result as CorrespondenceMapSubgroupsIndexed except it is constructed on the fly "<>
															 "by using the generators of \!\(\*SubscriptBox[\(Z\), \(n\)]\). (Algorithm 2 to compute correspondence instead of finding it.)";

Collisions::usage = " <| Int -> {{Int,Int}} |>  |  Parts of correspondence map that are non-unique. Where indexed subgroups could be mapped to other  "<>
												   "quotient groups. The latter are given as indices in the QuotientGroups matrix.";
CorrespondenceSizeFactors::usage = " <| Int -> Int |>  |  For a correspondence, how many times smaller the group in the range is than the group in the domain.";

CorrespondenceReduction::usage = " {{{Int},{Int}}}  |  The size reductions in the correspondence map given as a mapping between sizes. See CorrespondanceMap";
CorrespondenceMaxReductions::usage = " {Int}  |  The maximal size reduction for each subgroup in \!\(\*SubscriptBox[\(Z\), \(n\)]\) under the correspondence.;";

ThirdNonIsomorphiesN::usage = "";


(* ::Subsubsection:: *)
(*Isomorphism theorems*)


IsomorphicGroup::usage = " The group isomorphic to the group under a homomorphism.";
IsomorphicCyclicGroup::usage = " The cyclic group isomorphic to the cyclic group given. The homomorphism is implicit. ";

RelabelElements::usage = " Relabel elements of a group with as smal numbering as possible starting from zero as the smallest.";
ExtensionRelabeling::usage = " {{<|Int -> Int|>}}  |  Matrix consisting of maps for relabeling of elements to reveal the isomorphy of the left and right sides of the second isomorphy theorem."<>
													 " Where there are no relabeling ";
ExtensionTransformativeTerms::usage = " {<|Int->Int|>  |  List of the the non-zero-size indices in the extension relabeling matrix.";
ExtensionIsomorphy::usage = " {{Int}}  |  The cyclic isomorphic groups that concludes the isomorphies of the second isomorphism theorem.";

IntersectionsToProducts::usage = " {[Int}} -> {{Int}}  |  Proof-of-principle-methods of the reshufflability of the outputs of SubgroupProducts and SubgroupIntersection to each other. "<>
														 "So in \!\(\*SubscriptBox[\(Z\), \(n\)]\) groupwise multiplication and intersections is interchangeable, that is does the same thing in some sense."<>
														 " See AdditiveGroup`SubgroupProducts and AdditiveGroup`SubgroupIntersections.";
ProductsToIntersections::usage = IntersectionsToProducts::usage;
														 
QuotientSimilars::usage = " Int -> {{Int,Int}}  |  Which set of subggroups (indexes) gives (as denominators) the same quotient group for a given (nominator) subgroup index. ";

SecondIsomorphyDeterminingSubgroups::usage = "";

ThirdIsomorphyTrippelsClasses::usage = "";
ThirdIsomorphyTrippels::usage = "";
ThirdIsomorphyTrippelsIso::usage = "";
ThirdIsomorphyMappingIso::usage = "";
ThirdIsomorphyMappingIndexToIso::usage = "";																																																																																																																																																																																																					
ThirdIsomorphyCheck::usage = "";
ThirdNonIsomorphies::usage = "";
ThirdNonIsomorphiesFull::usage = "";
ThirdIsomorphies::usage = "";
ThirdNonIsomorphicK::usage = "";
ThirdNonIsomorphicM::usage = "";
ThirdNonIsomorphicKIndex::usage = "";
ThirdNonIsomorphicMIndex::usage = "";
ThirdIsomorphyFromGenerators::usage = "";
LargestThirdIsomorphyPair::usage = "";

FirstKMs::usage = "";
FirstKMindexed::usage = "";
FirstFalseKMindexed::usage = "";
FirstTrueKMindexed::usage = "";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Correspondence theorem*)


(* * * group side of correspondence theorem * * *)
ContainmentCorrespondenceGroups[k_]:=Take[#,Length[#]/Length[Sns[[k]]]]&/@Containment[k]
ContainmentCorrespondenceOrder[k_]:=Length[First[ContainmentCorrespondenceGroups[k]]]
AllContainmentCorrespondencies[]:= Table[ContainmentCorrespondenceGroups[i],{i,N1,1,-1}]
AllContainmentCorrespondenceOrders[]:= Table[Length[ContainmentCorrespondenceGroups[i]],{i,N1}]

(* * * quotient group side of correspondence theorem * * *)
QuotientCorrespondenceGroups[k_]:= Return[DeleteDuplicates[QuotientGroup[#,k]& /@ Range[k]]]
QuotientCorrespondenceOrder[k_]:= Return[Length[QuotientCorrespondenceGroups[k]]]
AllQuotientCorrespondencies[]:= Table[QuotientCorrespondenceGroups[i],{i,N1,1,-1}]
AllQuotientCorrespondenceOrders[]:= Table[Length[QuotientCorrespondenceGroups[i]],{i,N1}]
Remove[i,k];

CorrespondenceMap[k_]:= With[{HS=Sns},
									Module[{CS=Containment[k],subgroupindexes,HNsubset},
										subgroupindexes=Flatten[CS/.PositionIndex[HS]];
										HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
										Return[AssociationThread[HS[[subgroupindexes]],HNsubset]];
									]]

UniqueCorrespondenceMap[k_]:= With[{HS=Sns},
											Module[{Cs=Containment[k],HN=QuotientGroups[],subgroupindexes,quotientgroupindexes,HNsubset},
												subgroupindexes=Flatten[Cs/.PositionIndex[HS]];
												HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
												HN=If[MemberQ[Cs,#[[-1]]],#,{{}}]&/@ HN;
												quotientgroupindexes=Position[HN,#]& /@ HNsubset;
												quotientgroupindexes[[-1]]={Last[quotientgroupindexes[[-1]]]};
												Return[AssociationThread[subgroupindexes,HNsubset]];
												]]

CorrespondenceMapIndexed[k_]:= With[{HS=Sns,HN=QuotientGroups[]},
									Module[{Cs=Containment[k],subgroupindexes,quotientgroupindexes,HNsubset},
										subgroupindexes=Flatten[Cs/.PositionIndex[HS]];
										HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
										quotientgroupindexes=Position[HN,#]& /@ HNsubset;
										Return[AssociationThread[subgroupindexes,quotientgroupindexes]];
										]]

UniqueCorrespondenceMapIndexed[k_]:= With[{HS=Sns},
											Module[{Cs=Containment[k],HN=QuotientGroups[],subgroupindexes,quotientgroupindexes,HNsubset},
												subgroupindexes=Flatten[Cs/.PositionIndex[HS]];
												HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
												HN=If[MemberQ[Cs,#[[-1]]],#,{{}}]&/@ HN;
												quotientgroupindexes=Position[HN,#]& /@ HNsubset;
												quotientgroupindexes[[-1]]={Last[quotientgroupindexes[[-1]]]};
												Return[AssociationThread[subgroupindexes,quotientgroupindexes]];
												]]
Remove[HS,CS,k,subgroupindexes,quotientgroupindexes,HN,HNsubset];

Begin["`Extras`"];		
	AdditiveGroupTheorems`Extras`CorrespondenceMapSubgroupsIndexed::usage = " <| Int -> {Int}  |  Like CorrespondenceMap except the domain of subgroups are all-subgroups-sizes-indexed.";						
	AdditiveGroupTheorems`Extras`CorrespondenceMapSubgroupsIndexed[k_]:= With[{HS=Sns},
												Module[{Cs=Containment[k],subgroupindexes,HNsubset},
													subgroupindexes=Flatten[Cs/.PositionIndex[HS]];
													HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
													Return[AssociationThread[subgroupindexes,HNsubset]];
												]]
	AdditiveGroupTheorems`Extras`CorrespondenceMapCosetsIndexed::usage = " <| {Int} -> Int |>  |  Like CorrespondenceMap except the range of quotient groups are all-quotientgroup-sizes-indexed. ";
	AdditiveGroupTheorems`Extras`CorrespondenceMapCosetsIndexed[k_]:= With[{HS=Sns,HN=QuotientGroups[]},
											Module[{Cs=Containment[k],subgroupindexes,quotientgroupindexes,HNsubset},
												subgroupindexes=Flatten[Cs/.PositionIndex[HS]];
												HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
												quotientgroupindexes=Position[HN,#]& /@ HNsubset;
												Return[AssociationThread[Sns[[subgroupindexes]],quotientgroupindexes]];
												]]
End[];

CorrespondenceMethodDictionary[k_]:=With[{map=AdditiveGroupTheorems`Extras`CorrespondenceMapSubgroupsIndexed[k]},
										Module[{lengths = Length /@ Sns[[N1+1-Reverse[Keys[map]]]],valueranges},(* how to compute these lengths without generating subgroup? *)
												valueranges=Table[Take[Sns[[Keys[map]]][[i]],lengths[[i]]],{i,1,Length[lengths]}];
												Return[AssociationThread[Keys[map],valueranges]];
											]];

CorrespondenceMethodGenerators[k_]:=With[{map=AdditiveGroupTheorems`Extras`CorrespondenceMapSubgroupsIndexed[k]},
										Module[{lengths = Length /@ Sns[[N1+1-Reverse[Keys[map]]]],valueranges,generators},
												valueranges=Range[0,#-1]&/@lengths;
												generators=Keys[map];(* indexes is actually generators here - do it properly *)
												Return[AssociationThread[Keys[map],Inner[Times,generators,valueranges,List]]];
											]];
Remove[map,k,lengths,valueranges,generators];

Collisions[k_]:=With[{map=CorrespondenceMapIndexed[k]},Return[Select[map,(Length[#]>1)&]]];

CorrespondenceSizeFactors[k_]:=With[{map=CorrespondenceMap[k]},
						Return[AssociationThread[Keys[map],MapThread[Quotient[Length[#1],Length[#2]]&,{Keys[map],Values[map]}]]];
					]
Remove[k,map];

Begin["`Extras`"];
	CorrespondenceReport::usage = " Obscure investigation of the correspondences. To be removed.";
	CorrespondenceReport[k_]:= Module[{Ns,Cs,GN,Is,n1,HNs,
										dims,sizes,contis,contimap,corrimap,map,imap,
										images,partitioning,
										results,maps,report},
										
										Ns=Sns[[k]];
										Cs=Containment[k];
										GN=AdditiveGroupQuotients`Private`QuotientGroupFull[k];
										n1=Length[GN];
										Is=Subgroups[n1];
	
										dims=<|"|G|"->N0,"|N|"->Length[Ns],"|G\\N|"->n1,"N"->Ns|>;
										sizes=<|"{|C|}"->(Length/@Cs),"{|H|},{|\[Phi]N|}"->{Length/@Sns,Length/@Is}|>;
	
										contis=Flatten[Cs/.PositionIndex[Sns]];
										contimap=AssociationThread[contis,Sns[[contis]]];
										
										HNs=AdditiveGroupQuotients`Private`QuotientGroupFull[#,k]&/@contis;
										corrimap=AssociationThread[contis,HNs];
										map=AssociationThread[Sns[[contis]],HNs];
										
										images=({}\[Union]#)& /@ Map[Canonical[k,#]&,Cs,{2}];
										
										partitioning=<|"|\[Phi](\!\(\*SubscriptBox[\(H\), \(i\)]\))|"->Length/@images,
													   "{|\[Phi](\!\(\*SubscriptBox[\(h\), \(i\)]\))|}"->MatrixForm[Map[Length,images,{2}]],
													   "\[Phi](\!\(\*SubscriptBox[\(h\), \(i\)]\))"->MatrixForm[images]|>;
													   
										report={dims,sizes,partitioning};
										results=Values[report];
										maps={contimap,corrimap,map};
										Return[{results,MatrixForm[maps],report}];	 
									]
	Remove[Ns,Cs,GN,Is,n1,HNs,dims,sizes,contis,contimap,corrimap,map,imap,images,partitioning,results,maps,report];	

	CorrespondenceReductionList::usage = " {{{Int},{Int}}}  |  The size reductions in the correspondence map. See CorrespondanceMap";
	CorrespondenceReductionList[]:= Module[{subsets, quotientgroups,list={}},
											For[i=1,i<=N1,i++,
												subsets=SubsettingPaths[i];
												quotientgroups=subsets/.CorrespondenceMap[i];
												AppendTo[list,{Length \[Congruent] subsets,Length \[Congruent] quotientgroups}];
											];
											Return[list];
									]
	Remove[subsets,quotientgroups,list];
End[];

CorrespondenceReduction[]:= (MatrixForm[#[[1]]]->MatrixForm[#[[2]]])& /@ AdditiveGroupTheorems`Extras`CorrespondenceReductionList[]

CorrespondenceMaxReductions[]:= With[{all=AdditiveGroupTheorems`Extras`CorrespondenceReductionList[]}, Return[Last /@ (First /@ (First /@ all))];]
Remove[all];

ThirdNonIsomorphiesN[l_]:=Module[{nonisos,isoindexes,generators,orders,isos={}},
								For[i=10,i<100,i++,
										MakeGroup[i,True];
										nonisos=ThirdNonIsomorphicMIndex[#,N1+1-#]& /@ Range[1,N1];
										If[Length[Flatten[nonisos]]!=0,
											isoindexes = Select[nonisos,#!={}&];
											generators=SubgroupGenerators[][[#]]& /@ isoindexes;
											orders = GeneratorsAndOrder[][#]&\[Congruent] generators;
											AppendTo[isos,i-> {isoindexes,generators,orders}];
											,
											Nothing];
									];
									Return[isos]];
Remove[n,m,nonisos,isoindexes,generators,orders,isos];


(* ::Subsubsection:: *)
(*Isomorphism theorems*)


IsomorphicGroup[G_,op_]:=If[Length[G]==1,Return[G],Return[RelabelElements[G,op][[1]]]] (* first row if common modular addition *)

IsomorphicCyclicGroup[G_]:=If[Length[G]==1,Return[G],IsomorphicGroup[G,Mod[#,G[[-1]]+(G[[-1]]-G[[-2]])]& @* Plus]]
Remove[G,op]

RelabelElements[G_,op_]:=With[{n=Length[G],relabelingmap=AssociationThread[G,Range[0,Length[G]-1]]},
							Module[{multiplicationtable},
									multiplicationtable=Table[op[G[[i]],G[[j]]],{i,1,n},{j,1,n}];
									Return[multiplicationtable/.relabelingmap];
							]]
ExtensionRelabeling[]:=Table[
							long=(Sns[[n]]\[CirclePlus]Sns[[m]])\[Backslash]Sns[[m]]; 
							short=Sns[[n]]\[Backslash](Sns[[n]]\[Intersection]Sns[[m]]); 
							If[First /@ long == First /@ short,
								<||>,
								AssociationThread[First /@ short,First /@ long]],
									{n,1,N1},{m,1,N1}]
Remove[G,op,relabelingmap,multiplicationtable]

ExtensionTransformativeTerms[]:= Flatten[If[#!=<||>,#,Nothing]& \[Congruent] ExtensionRelabeling[]]

ExtensionIsomorphy[]:=Module[{long,short,longgroup,shortgroup},
						Table[
							long=(Sns[[n]]\[CirclePlus]Sns[[m]])\[Backslash]Sns[[m]]; 
							short=Sns[[n]]\[Backslash](Sns[[n]]\[Intersection]Sns[[m]]);
							longgroup=IsomorphicCyclicGroup[Sort[First /@ long]];
							shortgroup=IsomorphicCyclicGroup[Sort[First /@ short]]; 
							If[longgroup == shortgroup, shortgroup, False]
						,{n,1,N1},{m,1,N1}]]

IntersectionsToProducts[gs_]:= With[{indexes=(gs/.PositionIndex[Sns])},
								Module[{reverse = AssociationThread[Range[N1],Reverse[Range[N1]]],mapped}, 
									mapped=Reverse[Reverse /@ indexes]/.reverse;
									Return[First \[Congruent] (mapped/.AssociationThread[Range[N1],Sns])];			
									]]

ProductsToIntersections[gs_]:= With[{indexes = Flatten \[Congruent] (gs/.PositionIndex[Sns])},
								Module[{reverse = AssociationThread[Range[N1],Reverse[Range[N1]]],mapped}, 
									mapped=Reverse[Reverse /@ (indexes/.reverse)];									
									Return[First \[Congruent] (mapped/.AssociationThread[Range[N1],Sns])];			
									]]
Remove[indexes,gs,reverse,mapped,long,short,longgroup,shortgroup];

QuotientSimilars[gi_]:=Table[First /@ Position[Table[ContainsExactly[Sns[[j]]\[Backslash]Sns[[gi]],Sns[[i]]\[Backslash]Sns[[gi]]],{i,1,N1}],True],{j,1,N1}]
Remove[gi];

SecondIsomorphyDeterminingSubgroups[k_]:=If[Length[#]==Length[Sns[[1]]\[Backslash]Sns[[k]]],1,0]& \[Congruent](SubgroupIntersections[])
Remove[k];
										
ThirdIsomorphyTrippelsClasses[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsFlattened[n,m,k]},
											Module[{reps=AdditiveGroupTheorems`Package`RemoveNonDisjointSuccessors[pqQ[[3]][[1]]]},
													Return[{SortBy[pqQ[[1]],Min],SortBy[pqQ[[2]],Min],SortBy[reps,Min]}]]]

ThirdIsomorphyTrippels[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{pqQ[[1]],pqQ[[2]],First /@ pqQ[[3]]}]]

ThirdIsomorphyTrippelsIso[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{IsomorphicCyclicGroup[First /@ pqQ[[1]]],
															IsomorphicCyclicGroup[First /@ pqQ[[2]]],
															IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]

ThirdIsomorphyMappingIso[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{Sns[[n]] -> IsomorphicCyclicGroup[First /@ pqQ[[1]]],
															Sns[[m]] -> IsomorphicCyclicGroup[First /@ pqQ[[2]]],
															First /@ (Sns[[n]]\[Backslash]Sns[[m]])->IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]

ThirdIsomorphyMappingIndexToIso[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{n -> IsomorphicCyclicGroup[Sort[First /@ pqQ[[1]]]],
															m -> IsomorphicCyclicGroup[Sort[First /@ pqQ[[2]]]],
															First /@ (Sns[[n]]\[Backslash]Sns[[m]])->IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]
																																																																																																																																																																																																															
ThirdIsomorphyCheck[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{IsomorphicCyclicGroup[Sort[First /@ (Sns[[n]]\[Backslash]Sns[[m]])]],IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]

ThirdNonIsomorphies[n_]:= Module[{check}, Select[Table[check=ThirdIsomorphyCheck[n,i,j];If[check[[1]]!=check[[2]],{i,j},Nothing],{i,n+1,N1},{j,i,N1}],#!={}&]];
ThirdNonIsomorphiesFull[n_]:= Module[{check}, Select[Table[check=ThirdIsomorphyCheck[n,i,j];If[check[[1]]!=check[[2]],{{i,j},check},Nothing],{i,n+1,N1},{j,i,N1}],#!={}&]];
ThirdIsomorphies[n_]:= Module[{check}, Select[Table[check=ThirdIsomorphyCheck[n,i,j];If[check[[1]]==check[[2]],{i,j},Nothing],{i,n+1,N1},{j,i+1,N1}],#!={}&]];

ThirdNonIsomorphicK[n_,m_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,m,i];If[checks[[1]]!=checks[[2]],Sns[[i]],Nothing],{i,N1,m+1,-1}]]
ThirdNonIsomorphicM[n_,k_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,i,k];If[checks[[1]]!=checks[[2]],Sns[[i]],Nothing],{i,k,n+1,-1}]]

ThirdNonIsomorphicKIndex[n_,m_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,m,i];If[checks[[1]]!=checks[[2]],i,Nothing],{i,N1,m+1,-1}]]
ThirdNonIsomorphicMIndex[n_,k_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,i,k];If[checks[[1]]!=checks[[2]],i,Nothing],{i,k,n+1,-1}]]
Remove[check,checks,pqQ,reps]

ThirdIsomorphyFromGenerators[m_,k_]:= Module[{K,M,N},
													K=AdditiveGroup`Private`TotalSpan[{k},N0];
													M=K\[CirclePlus]AdditiveGroup`Private`TotalSpan[{m},N0];
													N=Zn;
													Return[ThirdIsomorphyMappingIso[Commons`Pos[N,Sns],
																					Commons`Pos[M,Sns],
																					Commons`Pos[K,Sns]]];
												]
Remove[k,m,n,K,M]
												
LargestThirdIsomorphyPair[]:= With[{gens=NonUnitMinimalGenerators[]},
							Module[{isos=Table[Values[ThirdIsomorphyFromGenerators[gens[[i]],gens[[j]]][[3]]],{i,1,Length[gens]},{j,1,Length[gens]}],
									lengths,max,generatorpairs},
									lengths=Length \[Congruent] isos;
									lengths=lengths-DiagonalMatrix[Diagonal[lengths]];
									max=Commons`NthLargest[2,lengths];(* squaring *);
									If[max==0,lengths=Length \[Congruent] isos; max=Commons`NthLargest[2,lengths],None];
									generatorpairs={gens[[#[[1]]]],gens[[#[[2]]]]}& /@ Position[lengths,max];
									Return[DeleteDuplicatesBy[generatorpairs,Sort]];
							]];
Remove[gens,isos,lengths,max,generatorpairs];

FirstKMs[]:=With[{pairs=LargestThirdIsomorphyPair[]},
					Return[{AdditiveGroup`Private`TotalSpan[{#[[1]]},N0],AdditiveGroup`Private`TotalSpan[{#[[2]]},N0]}& /@ pairs];
					]
Remove[pairs];
					
FirstKMindexed[]:=Flatten[FirstKMs[]/.PositionIndex[Sns],{1,3}];

FirstFalseKMindexed[]:=With[{nis=ThirdNonIsomorphies[1]},
						Module[{kms,false={},pair},
								kms=FirstKMindexed[];
								For[i=1,i<=Length[kms],i++,
									pair=kms[[i]];
									If[ContainsAny[Flatten[nis,1],{pair}],AppendTo[false,pair],None];
									];
								Return[false];]];

FirstTrueKMindexed[]:=With[{is=ThirdIsomorphies[1]},
						Module[{kms,true={},pair},
								kms=FirstKMindexed[];
								For[i=1,i<=Length[kms],i++,
									pair=kms[[i]];
									If[ContainsAny[Flatten[is,1],{pair}],AppendTo[true,pair],None];
									];
								Return[true];]];
Remove[i,j,k,l,nis,is,kms,true,false,pair];


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];

AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsDirect[n_,m_,k_]:= With[{p=Sns[[n]]\[Backslash]Sns[[k]], q=Sns[[m]]\[Backslash]Sns[[k]]},Return[{p,q,p\[Backslash]q}]]

AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsNoDups[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsDirect[n,m,k]},
												Return[{pqQ[[1]],pqQ[[2]],DeleteDuplicates[Map[DeleteDuplicates,pqQ[[3]],{2}]]}];]

AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsFlattened[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsNoDups[n,m,k]},
												Return[{pqQ[[1]],pqQ[[2]],Flatten[pqQ[[3]],{3}]}]]
Remove[n,m,k,p,q,pqQ];

AdditiveGroupTheorems`Package`RemoveNonDisjointSuccessors[ls_]:= Module[{keep,kept={},exhaust=ls},
										keep=First[exhaust];
										exhaust=Rest[exhaust];
										AppendTo[kept,keep];
										While[exhaust!={},
											keep=First[exhaust];
											exhaust=Rest[exhaust];
											If[Fold[Or,False,ContainsAny[keep,#]& /@ kept],
												None,
												AppendTo[kept,keep]];
										];
										Return[kept];
									] 
Remove[keep,kept,exhaust,ls];

End[];


(* ::Subsection:: *)
(*Seldom used*)


ExtensionRelabelingIndexed[]:=Table[
							long=(Sns[[n]]\[CirclePlus]Sns[[m]])\[Backslash]Sns[[m]]; 
							short=Sns[[n]]\[Backslash](Sns[[n]]\[Intersection]Sns[[m]]); 
							If[First /@ long == First /@ short,
								First /@ long,
								{{n,m},AssociationThread[First /@ short,First /@ long]}],
									{n,1,N1},{m,1,N1}]
Remove[n,m,long,short]


(* ::Section:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
