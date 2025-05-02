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


AdditiveGroupTheoremsPackage::usage = "This is is the fifth module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suite, the AdditiveGroup package suite. It implements the "<>
									"correspondence theorem and the three isomorphism theorems with some adjoining functionality.";

Print["AdditivegroupTheorems`: See Docs[\"Theorems\"] for documentation."]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Correspondence theorem*)


CorrespondenceMap::usage = " <|{Int} -> {Int}|> |   The correspondence map for a given subgroup. Correspondence between subgroups containing the subgroup "<>
													"and quotient groups derived by letting larger groups taking quotients with it.";			
															
CorrespondenceMapIndexed::usage = " <| Int -> Int |>  |  Like CorrespondenceMap but all non-unique occurrences is displayed since indexes indexing same groups are used.";

CorrespondenceImageSizes::usage = " {Int}  |  The sizes )not orders( of all subgroups on the quotient group side of the correspondence theorem, the 'other' side.";

CorrespondenceFullGroupImages::usage = " {Int}  |  The images of the full group, that is on the quotient group side of the correspondence theorem, the 'other' side.";

Collisions::usage = " <| Int -> {{Int,Int}} |>  |  Parts of correspondence map that are non-unique. Where subgroups are mapped to by many  "<>
												   " quotient groups. The latter are given as indices in the QuotientGroups matrix.";


(* ::Subsubsection:: *)
(*Isomorphism theorems*)


ExtensionRelabeling::usage = " {{<|Int -> Int|>}}  |  Matrix consisting of maps for relabeling of elements to reveal the isomorphy of the left and right sides of the second isomorphy theorem."<>
													 " Where there are no relabeling ";

ExtensionTransformativeTerms::usage = " {<|Int->Int|>  |  List of the the non-zero-size indices in the extension relabeling matrix.";

ExtensionIsomorphy::usage = " {{Int}}  |  The cyclic isomorphic groups that concludes the isomorphies of the second isomorphism theorem.";

IntersectionsToProducts::usage = " {[Int}} -> {{Int}}  |  Proof-of-principle-methods of the reshufflability of the outputs of SubgroupProducts and SubgroupIntersection to each other. "<>
														 "So in \!\(\*SubscriptBox[\(Z\), \(n\)]\) groupwise multiplication and intersections is interchangeable, that is does the same thing in some sense."<>
														 " See AdditiveGroup`SubgroupProducts and AdditiveGroup`SubgroupIntersections.";
														 
ProductsToIntersections::usage = IntersectionsToProducts::usage;
														 
SecondIsomorphyLiftFactors::usage = " {{Int}}  |  The scalar quotients giving the expansion of groups going from lower side (intersections) to upper side (products) in the second isomorphy theorem "<>
												" or the upper side of the diamond graph.";
																											
ThirdIsomorphyTrippels::usage = " Int,Int,Int -> {{Int},{Int},{int}}  |  The representatives set of the nominator, denominator and the quotient on the other "<>
														" side of \!\(\*SubscriptBox[\(Z\), \(n\)]\)\\!\(\*SubscriptBox[\(Z\), \(m\)]\) in the equation of the third isomorphy theorem.";

ThirdIsomorphyMappingIndexToIso::usage = " Int,Int,Int -> {Int->{Int},Int->{Int},int->{Int}}  |  Mapping between the respective nominators, denominators and the full quotient "<>
																				" the two sides of the third isomorphy theorem. The domain groups are indexed ny size.";

ThirdIsomorphyCheck::usage = " Int,Int,Int -> Bool  |  Checks if the third isomorphy theorem is true, that is if i\.08t's assumptions is true, that is if the trippels are "<>
													"(normal) subgrouping sequence.";

ThirdNonIsomorphies::usage = " Int -> {{Int,Int}}  |  For a subgroup index it returns all pairs of group indexes for which the third isomophy theorem's assumptions are fullfilled."; 

ThirdIsomorphies::usage = " Int -> {{Int,Int}}  |  The negation of ThirdNonIsomorphies, that is where the theorem's assumptions are not fullfilled. " ;

ThirdNonIsomorphicKIndex::usage = " Int,Int -> Int |  Returns one out of three groups index in (\!\(\*SubscriptBox[\(Z\), \(n\)]\)\\K)\(M\\K) for which the third isomorphy theorem does not hold "<>
													" and where the group with index m,k gives this method postfix M,K.";

ThirdNonIsomorphicMIndex::usage = " Int,Int -> Int |  Returns one out of three groups index in (\!\(\*SubscriptBox[\(Z\), \(n\)]\)\\K)\(M\\K) for which the third isomorphy theorem does not hold "<>
													" and where the group with index m,k gives this method postfix M,K.";

ThirdIsomorphyFromGenerators::usage = " Int,Int -> {{Int}->{Int},{Int}->{Int},{int}->{Int}}  |  Like ThirdIsomorphyMappingIso but takes two generators, "<>
																								" generating two subgroups, as arguments.";

LargestThirdIsomorphyGeneratorPairs::usage = " {Int,Int}  |  Returns pairs of generators of different one-generator-subgroups that generates the largest and equal size of the "<>
															"left and right side of the third isomorphy theorem.";

FirstLargestThirdIsomorphyGroupsIndexed::usage = " {{Int,Int}}  |  Groups for which left side isomorphic group equals right side isomophic group. Groups are indexed asu usual.";

FirstFalseThirdIsomorphyGroupsIndexed::usage = " {{Int,Int}}  |  Like FirstLargestThirdIsomorphyGroupsIndexed but the pair(s) are not isomporphic.";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Theory*)


(* ::Subsubsection:: *)
(*Correspondence theorem*)


CorrespondenceMap[k_]:= With[{HS=Sns},
									Module[{CS=Containment[k],subgroupindexes,HNsubset},
										subgroupindexes=Flatten[CS/.PositionIndex[HS]];
										HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
										Return[AssociationThread[HS[[subgroupindexes]],HNsubset]];
									]]
Remove[HS,CS,subgroupindexesS,HNsubset];
									
CorrespondenceMapIndexed[k_]:= With[{HS=Sns,HN=QuotientGroups[]},
									Module[{CS=Containment[k],subgroupindexes,quotientgroupindexes,HNsubset},
										subgroupindexes=Flatten[CS/.PositionIndex[HS]];
										HNsubset=QuotientGroup[#,k]&/@subgroupindexes;
										quotientgroupindexes=Position[HN,#]& /@ HNsubset;								
										Return[AssociationThread[subgroupindexes,quotientgroupindexes]];
										]]
Remove[HS,CS,subgroupindexes,quotientgroupindexes,HN,HNsubset];

CorrespondenceImageSizes[]:= Length /@ (CorrespondenceMap[#]& /@ Range[N1]);

CorrespondenceFullGroupImages[]:= (First /@ (CorrespondenceMap[#]& /@ Range[N1]));

Collisions[k_]:=With[{map=CorrespondenceMapIndexed[k]},Return[Select[map,(Length[#]>1)&]]];
Remove[map];


(* ::Subsubsection:: *)
(*Isomorphism theorems*)


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
																																
SecondIsomorphyLiftFactors[]:= { Table[(Length\[Congruent]SubgroupProducts[])[[i]]/(Length /@ Sns),{i,1,N1}],
								 Table[Length[Sns[[i]]]/(Length\[Congruent]SubgroupIntersections[])[[i]],{i,1,N1}]}
						
ThirdIsomorphyTrippels[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{First /@ pqQ[[1]], First /@ pqQ[[2]],First /@ pqQ[[3]]}]]

ThirdIsomorphyMappingIndexToIso[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{n -> IsomorphicCyclicGroup[Sort[First /@ pqQ[[1]]]],
															m -> IsomorphicCyclicGroup[Sort[First /@ pqQ[[2]]]],
															First /@ (Sns[[n]]\[Backslash]Sns[[m]])->IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]

ThirdIsomorphyCheck[n_,m_,k_]:= With[{pqQ=AdditiveGroupTheorems`Private`ThirdIsomorphyTrippelsClasses[n,m,k]},
													Return[{IsomorphicCyclicGroup[Sort[First /@ (Sns[[n]]\[Backslash]Sns[[m]])]],IsomorphicCyclicGroup[First /@ pqQ[[3]]]}]]

ThirdIsomorphies[n_]:= Module[{check}, Select[Table[check=ThirdIsomorphyCheck[n,i,j];If[check[[1]]==check[[2]],{i,j},Nothing],{i,n+1,N1},{j,i+1,N1}],#!={}&]];
ThirdNonIsomorphies[n_]:= Module[{check}, Select[Table[check=ThirdIsomorphyCheck[n,i,j];If[check[[1]]!=check[[2]],{i,j},Nothing],{i,n+1,N1},{j,i,N1}],#!={}&]];

ThirdNonIsomorphicKIndex[n_,m_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,m,i];If[checks[[1]]!=checks[[2]],i,Nothing],{i,N1,m+1,-1}]]
ThirdNonIsomorphicMIndex[n_,k_]:= Module[{checks},Table[checks=ThirdIsomorphyCheck[n,i,k];If[checks[[1]]!=checks[[2]],i,Nothing],{i,k,n+1,-1}]]
Remove[check,checks,pqQ,reps]

ThirdIsomorphyFromGenerators[m_,k_]:= Module[{K,M,N,mapping},
													K=AdditiveGroup`Private`TotalSpan[{k}];
													M=K\[CirclePlus]AdditiveGroup`Private`TotalSpan[{m}];
													N=Zn;
													mapping = ThirdIsomorphyMappingIndexToIso[Commons`Pos[N,Sns],Commons`Pos[M,Sns],Commons`Pos[K,Sns]];
													Return[(If[ListQ[Keys[#]],Keys[#]->Values[#],Sns[[Keys[#]]]->Values[#]])& /@ mapping] (* not just indexes but groups as keys *)
												]
Remove[m,n,K,M,mapping]
												
LargestThirdIsomorphyGeneratorPairs[]:= With[{gens=Drop[SubgroupGenerators[],2]},
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

FirstLargestThirdIsomorphyGroupsIndexed[]:=With[{pairs=LargestThirdIsomorphyGeneratorPairs[]},
												Module[{groups},
														groups={AdditiveGroup`Private`TotalSpan[{#[[1]]}],AdditiveGroup`Private`TotalSpan[{#[[2]]}]}& /@ pairs;
														Return[Flatten[groups/.PositionIndex[Sns],{1,3}]];
												]];
Remove[pairs,groups];
					
FirstFalseThirdIsomorphyGroupsIndexed[]:=With[{nis=ThirdNonIsomorphies[1]},
						Module[{kms,false={},pair},
								kms=FirstLargestThirdIsomorphyGroupsIndexed[];
								For[i=1,i<=Length[kms],i++,
									pair=kms[[i]];
									If[ContainsAny[Flatten[nis,1],{pair}],AppendTo[false,pair],None];
									];
								Return[false];]];
Remove[i,j,k,l,nis,is,kms,true,false,pair];


(* ::Subsection::Closed:: *)
(*Helpers*)


Begin["`Private`"];

ThirdIsomorphyTrippelsClasses::usage = " Int,Int,Int -> {{{Int}},{{Int}},{{int}}}  |  The full coset set of the nominator, denominator and quotient on the other "<>
																		"side of \!\(\*SubscriptBox[\(Z\), \(n\)]\)\\!\(\*SubscriptBox[\(Z\), \(m\)]\) in the equation of the third isomorphy theorem.";

ThirdIsomorphyTrippelsClasses[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsFlattened[n,m,k]},
											Module[{reps=RemoveNonDisjointSuccessors[pqQ[[3]][[1]]]},
													Return[{SortBy[pqQ[[1]],Min],SortBy[pqQ[[2]],Min],SortBy[reps,Min]}]]]

ThirdIsomorphyTrippelsDirect[n_,m_,k_]:= With[{p=Sns[[n]]\[Backslash]Sns[[k]], q=Sns[[m]]\[Backslash]Sns[[k]]},Return[{p,q,p\[Backslash]q}]]

ThirdIsomorphyTrippelsNoDups[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsDirect[n,m,k]},
												Return[{pqQ[[1]],pqQ[[2]],DeleteDuplicates[Map[DeleteDuplicates,pqQ[[3]],{2}]]}];]

ThirdIsomorphyTrippelsFlattened[n_,m_,k_]:= With[{pqQ=ThirdIsomorphyTrippelsNoDups[n,m,k]},
												Return[{pqQ[[1]],pqQ[[2]],Flatten[pqQ[[3]],{3}]}]]
Remove[n,m,k,p,q,pqQ];

RemoveNonDisjointSuccessors[ls_]:= Module[{keep,kept={},exhaust=ls},
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


(* ::Section::Closed:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
