(* ::Package:: *)

(* ::Section:: *)
(*The Additive Group Quotientgroup Package*)


PrependTo[$ContextPath,"Commons`"];
PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
PrependTo[$ContextPath,"AdditiveGroupBasics`"];
PrependTo[$ContextPath,"AdditiveGroup`"];

BeginPackage["AdditiveGroupQuotients`"];
<< AdditiveGroup`


AdditiveGroupQuotientsPackage::usage = "This is is the fourth module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit, the AdditiveGroup package suite. It adds quotient group functionality which "<>
										" takes theory to higher grounds and prepares for the classical theorems.";
										
Print["AdditivegroupQuotients`: See Docs[\"Quotients\"] for documentation."]
Off[Image::shdw]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection:: *)
(*Operators, constants  and primitive set and mappings	*)


Css::usage = " {{{Int}}}  |  Cosets of all quotient groups constructible by the subgroups of Zn.";

CircleDot::usage = "  {Int},{Int} --> {Int} |  The operator giving the product set (subgroup) of two sets (subgroups) of sets(cosets). "<>
											  " All members (cosets) are added to each other respectively.";

SmallCircle::usage = "  {Int},{Int} --> {Int} |  Like CircleDot but  uses \[Diamond] instead of \[CirclePlus].";

Backslash::usage = "  {Int},{Int} --> {{Int}} |  The quotient group of two groups. Most often the second argument is a subgroup of the first but most important is that they both have a common supergroup.";

Canonical::usage = " Int --> {Int}  |  The canonical isomophism which takes an element to the coset it is a representative for.";

Kernel::usage = "  {Int}  |  The elements in the domain that a homomorphism maps to the zero element in the image.";	

QuotientKernels::usage = " {{Int}}  |  The kernels of all quotient groups. ";

QuotientKernelSpacings::usage = " {{Int}}  |  The distance (difference) between elements of all quotient group's kernels.";

Extend::usage = " Brevity operator for investigating left hand side of second isomorphy theorem. ";

Extending::usage = " Brevity operator for investigating coset multiplication. ";


(* ::Subsection:: *)
(*Whole group*)


MakeGroup::usage = " Int,Bool --> Void  |  Constructs \!\(\*SubscriptBox[\(Z\), \(n\)]\) it's subgroups and if second argument is True also all the group's quotient groups are prepared.";

IsNormal::usage = " Since \!\(\*SubscriptBox[\(Z\), \(n\)]\) is cyclic all groups concerned are Abelian. ";


(* ::Subsection:: *)
(*Instances*)


MakeFullGroupInstance::usage = " Int --> {{Int},Int,Int,Int,{Int},{Int}}  |  Computes an additive group instance containing it's subgroups and also it's quotient groups. This function temporarily "<>
																			" alters but restores the current context. If you only need it's subgroup use MakeMinimalGroupInstance instead. ";
																		
InstanceQuotientgroups::usage = " Int -->  {{Int}}  |  Returns the quotient groups of the additive group of a given order. It temporarily alters but restores the current context group.";

QuotientPermutationGroup::usage  =  " Int|Int,Int --> PermutationGroup  |  The mathematica permutation group isomorphic with this package's quotient group of two subgroups.";


(* ::Subsection:: *)
(*Derived sets*)


Coset::usage = " Int,Int,Int --> {Int}  |  The coset of \!\(\*SubscriptBox[\(Z\), \(n\)]\) with representant g.";

Cosets::usage = " Int,Int[,Int] --> {{Int}}  |  The cosets of two subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) if left and right cosets are same, otherwise -1. If no last argument the whole group is used.";

CosetsPowerSets::usage = " The set of all possible cosets in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

QuotientGroup::usage = " Int,Int[,Int] --> {Int}  |  The quotient group of two subgroups in \!\(\*SubscriptBox[\(Z\), \(n\)]\). If no last argument the whole group is used.";

QuotientGroups::usage = " {{Int}}  |  All the quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

Indexes::usage = " {Int}  |  All indexes in relation to the full group \!\(\*SubscriptBox[\(Z\), \(n\)]\). See Index.";

QuotientGroupGenerators::usage = " {Int}  |  The generators of the quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

QuotientGroupsRanks::usage = " {Int} | NOT IMPLEMENTED YET. WAITING FOR ELEGANT WAY TO COMPUTE IT. NOT JUST GENERATE IT BY BRUTE FORCE.";

RelabelElements::usage = " Relabel elements of a group with as smal numbering as possible starting from zero as the smallest.";

IsomorphicCyclicGroup::usage = " The cyclic group isomorphic to the cyclic group given. The homomorphism is implicit. ";


(* ::Subsection::Closed:: *)
(*Structure and graphical overview*)


CosetSizeExtremes::usage = " {-1|0|1}  |  Indications of local maximum and minimum values of the coset lengths for all quotient groups.";


(* ::Section:: *)
(*Code*)


(* ::Subsection::Closed:: *)
(*Operators, constants  and primitive sets and mappings*)


CircleDot[gs1_,gs2_]:= DeleteDuplicates@Sort@Flatten[Outer[CirclePlus,gs1,gs2,1],1]
Remove[gs1,gs2];

SmallCircle[gs1_,gs2_]:= Flatten[Map[DeleteDuplicates@*Flatten,Outer[Diamond,gs1,gs2,1],{2}],{1,2}]
Remove[gs1,gs2];

Backslash[H1_,H2_]:=DeleteDuplicates[Sort/@Table[(H1[[i]]\[CirclePlus]#)&/@H2,{i,1,Length[H1]}]]

Canonical[k_,g_]:= With[{Ns=Sns[[k]]},Return[Sort[(g\[CirclePlus]#)&/@Ns]]]
Remove[k,g,Ns];

Kernel[phi_,k_]:= Module[{ker={},zero=Sns[k],i},
						For[i=1,i<=N0,i++,
							If[MemberQ[phi[k,Zn[[i]]],0],AppendTo[ker,Zn[[i]]],None];];
						Return[ker];]
						
ClearAttributes[Image,{Protected,ReadProtected}];
Image[phi_]:= Module[{},Return[Table[phi[k,#]&/@Zn,{k,1,N1}]];]
Image::usage = " {Int}  |  The set of all values that a homomorphism takes over it's domain.";
Remove[phi,k,ker,zero];

QuotientKernels[]:= Table[Select[Sns[[i]],(First[Canonical[j,#]]== 0)&],{i,1,N1},{j,1,N1}];

QuotientKernelSpacings[]:= ((If[#=={},0,First[#]])& \[Congruent]( Differences \[Congruent] QuotientKernels[]));

Extend[H1_,H2_]:=(H1\[CircleDot]H2)\[Backslash]H2;

Extending[H1_,g_]:=(g\[CirclePlus]#)& /@ H1;
Remove[H1,H2,Ns,g];


(* ::Subsection:: *)
(*Whole group*)


MakeGroup[n_,quotientgroups_]:= Module[{t=0},
											If[quotientgroups,
												t = Commons`EstimatedTimeCosets[n];
												If[t>60, 
													If[t>3600, 
														Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/3600]<>" hours."];
														,
														Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/60]<>" minutes."];
													],
													Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t]<>" seconds."]];	
												,
												t = Commons`EstimatedTime[n];
												If[t>60, 
													If[t>3600, 
														Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/3600]<>" hours."];
														,
														Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t/60]<>" minutes."];
													],
													Commons`PrintIf[Commons`TIMING,"Estimated time is: "<>ToString[t]<>" seconds."]
												];	
											];
											t=TimeUsed[]; 
											N0=n; 
											Zn=Range[0,N0-1]; 
											Sns=Subgroups[];
											N1=Length[Sns];
											Css={};
											If[quotientgroups,
												Css=CosetsPowerSets[];
												N2=Length[Css];
												,
												None];
											t=TimeUsed[]-t;
											If[t > 60, 
												If[t > 3600, Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t/3600]<>" hours."];,
															 Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t/60]<>" minutes."];],
															 Commons`PrintIf[Commons`TIMING,"Computation time: "<>ToString[t]<>" seconds."]];	
											runtime={N[N0],N[t]};
											Save[$HomeDirectory<>"\\makegrouplog.ma",runtime];
										];
					
Remove[t,quotientgroups,n,runtime];

IsNormal = True; (* all cyclic groups are normal *)


(* ::Subsection:: *)
(*Instances*)


MakeFullGroupInstance[n_,timing_]:= Module[{C0,C1},
											C0 = {Zn,N0,N1,N2,Sns,Css};
											MakeGroup[n,timing];
											C1 = {Zn,N0,N1,N2,Sns,Css};
											{Zn,N0,N1,N2,Sns,Css}=C0;
											Return[C1];
										]
								
InstanceQuotientgroups[n_,timing_]:= Map[First,MakeFullGroupInstance[n,timing][[6]],{3}];
Remove[C0,C1,n,x1,x2,timing];

QuotientPermutationGroup[n_,m_]:= PermutationGroup[CyclesMap[] /@ QuotientGroup[n,m]]
QuotientPermutationGroup[m_]:= PermutationGroup[CyclesMap[] /@ QuotientGroup[N1,m]]
Remove[n,m];																																										


(* ::Subsection:: *)
(*Derived sets*)


Cosets[k_,l_]:= Return[AdditiveGroupQuotients`Private`CosetsLeft[k,l]]
Cosets[k_]:= Cosets[1,k]
Remove[k,l];

Coset[k_,r_]:= SelectFirst[Cosets[k],MemberQ[#,r]&]
Remove[k,r];

CosetsPowerSets[]:= If[Css=={},Table[Cosets[i,j],{i,1,N1},{j,i,N1}],Return[Css]]
											
QuotientGroup[k_,l_]:= First /@ AdditiveGroupQuotients`Private`QuotientGroupFull[k,l]
 
QuotientGroups[]:= Map[First,AdditiveGroupQuotients`Private`QuotientGroupsFull[],{3}]

Indexes[]:= AdditiveGroupQuotients`Private`Index[1,#]&/@ Range[N1];
																																				
QuotientGroupGenerators[]:= (#[[2]])& \[Congruent](Rest /@ QuotientGroups[][[1;;-2]])

QuotientGroupsRanks[] := Return["Not implemented"];

RelabelElements[G_,op_]:=With[{n=Length[G],relabelingmap=AssociationThread[G,Range[0,Length[G]-1]]},
							Module[{multiplicationtable},
									multiplicationtable=Table[op[G[[i]],G[[j]]],{i,1,n},{j,1,n}];
									Return[multiplicationtable/.relabelingmap];
							]]

IsomorphicCyclicGroup[G_]:=If[Length[G]==1,Return[G],AdditiveGroupQuotients`Private`IsomorphicGroup[G,Mod[#,G[[-1]]+(G[[-1]]-G[[-2]])]& @* Plus]]


(* ::Subsection::Closed:: *)
(*Structure and graphical overview*)


CosetSizeExtremes[]:= AdditiveGroupQuotients`Private`CosetSizeMaxima[] + AdditiveGroupQuotients`Private`CosetSizeMinima[]
Remove[L];


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];
	IsomorphicGroup::usage = " The group isomorphic to the group under a homomorphism.";
	IsomorphicGroup[G_,op_]:=If[Length[G]==1,Return[G],Return[RelabelElements[G,op][[1]]]] (* first row if common modular addition *)

	CosetsLeft::usage = " Int,Int[,Int] --> {{Int}}  |  The left cosets of two subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\). If no last argument the whole group is used.";
	CosetsLeft[k_,l_]:= With[{HS=Sns,G=Sns[[k]],H=Sns[[l]]},
										If[k<=N1\[Or]l>=k,None,Print["Indexing error"]; Return[-1]];
										Module[{g,coset={},cosets={},i},
											For[i=1,i<=Length[G],i++,
												g=G[[i]];
												coset=Sort[(g\[CirclePlus]#)&/@H];
												cosets=cosets\[Union]{coset};
												coset={};
											];
											Return[cosets];]]
	CosetsLeft[k_]:= CosetsLeft[1,k]
	Remove[HS,G,H,k,l,g,coset,cosets];
	
	CosetSizeMinima[]:= With[{L=CosetSizes[]}, 
									Table[
										If[AdditiveGroupQuotients`Private`dipp[#],-1,0,0]& /@ Partition[L[[i]],3,1],
											{i,1,N1-2}]]
	CosetSizeMaxima[]:= With[{L=CosetSizes[]}, 
									Table[
										If[AdditiveGroupQuotients`Private`peak[#],1,0,0]& /@ Partition[L[[i]],3,1],
											{i,1,N1-2}]]

	Index[k_,l_]:= If[IsNormal,Return[Length[CosetsLeft[k,l]]],Print["not commutative group"];Return[-1]]
	Index[k_]:= Index[1,k];
	Remove[k,l];

	dipp::usage = " Int,Int,Int  --> True|False  |  Returns true if the middle value is strictly the smallest of the three given.";
	peak::usage = " Int,Int,Int  --> True|False  |  Returns true if the middle value is strictly the largest of the three given.";
	
	dipp[xs_]:=(xs[[2]]<xs[[1]]&&xs[[2]]<xs[[3]])
	peak[xs_]:=(xs[[2]]>xs[[1]]&&xs[[2]]>xs[[3]])
	Remove[xs];
	
	QuotientGroupFull[k_,l_]:=Module[{Q},
										If[Css=={}, 
											Q=Cosets[k,l]; N2=Length[Q]; Return[Q]
											, 
											Return[Css[[k,l-k+1]]]]]
	
	QuotientGroupsFull[]:= If[Css=={}, Return[Table[QuotientGroupFull[i,j],{i,1,N1},{j,i,N1}]],Return[Css]];
	
	QuotientGroupFull[l_]:= QuotientGroupFull[1,l]
	Remove[k,l,Q];
	
	QuotientGroupFullTable[n_,m_] := With[{H=QuotientGroupFull[n,m]},Module[{mtable},
										mtable=Table[H[[i]]\[SmallCircle]H[[j]],{i,1,Length[H]},{j,1,Length[H]}];
										Return[mtable];
									]]
	Remove[n,m,H,mtable];

	
	Begin["`SeldomUsed`"];
		CosetsLeftInvertable[k_,l_]:= Module[{HS=Sns,G,H,g,coset={},cosets={},i},
											If[k<=N1,G=HS[[k]],Print["Indexing error"];Return[-1];];
											If[l<=N1,H=HS[[l]],Print["Indexing error"];Return[-1];];
											For[i=1,i<=Length[G],i++,
												g=G[[i]];
												coset=Sort[(g\[CirclePlus]#)&/@H];
												cosets=cosets\[Union]{coset};
												coset={};
											];
											Return[cosets];]
		Remove[HS,G,H,g,coset,cosets,k,l];	
																									
		CosetsPowerSetInvertable[]:= Table[CosetsLeftInvertable[i,j],{i,1,N1},{j,1,N1}]
		
		QuotientGroupsPartitionFull[]:= Table[
										{(AdditiveGroupQuotients`Private`QuotientGroupFull[#,i]& /@ Range[1,i]),
										None,
										(AdditiveGroupQuotients`Private`QuotientGroupFull[i,#]& /@ Range[i,N1])},
										{i,1,N1}]
		
		Remove[i,j]
	End[];

End[];


(* ::Section::Closed:: *)
(*Author: Anders Persson (persssonandersper@gmail.com)*)


EndPackage[];
SetAttributes[System`Image,{Protected,ReadProtected}]
On[Image::shdw]
