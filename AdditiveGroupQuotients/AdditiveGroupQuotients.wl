(* ::Package:: *)

(* ::Section::Closed:: *)
(*The Additive Group Quotientgroup Package*)


PrependTo[$ContextPath,"Commons`"];
PrependTo[$ContextPath,"AdditiveGroupMinimal`"];
PrependTo[$ContextPath,"AdditiveGroupBasics`"];
PrependTo[$ContextPath,"AdditiveGroup`"];

BeginPackage["AdditiveGroupQuotients`"];
<< AdditiveGroup`


AdditiveGroupQuotientsPackage::usage = "This is is the fourth module of the \!\(\*SubscriptBox[\(Z\), \(n\)]\) package-suit, the AdditiveGroup package suite. It adds quotient group functionality which "<>
										" takes theory to higher grounds and prepares for the classical theorems.";


(* ::Section::Closed:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*Operators, constants  and primitive mappings	*)


Unprotect[Image];

N2::usage = " Int  |  The in-use quotient group's size and order."; 
Css::usage = " {{{Int}}}  |  Cosets of all quotient groups constructible by the subgroups of Zn.";

CircleDot::usage = "  {Int},{Int} --> {Int} |  The operator giving the product set (subgroup) of two sets (subgroups) of sets(cosets). "<>
											  " All members (cosets) are added to each other respectively.";

SmallCircle::usage = "  {Int},{Int} --> {Int} |  Like CircleDot but  uses \[Diamond] instead of \[CirclePlus].";

Backslash::usage = "  {Int},{Int} --> {{Int}} |  The quotient group of two groups. Most often the second argument is a subgroup of the first but most important is that they both have a common supergroup.";

Canonical::usage = " Int --> {Int}  |  The canonical isomophism which takes an element to the coset it is a representative for.";
Kernel::usage = "  {Int}  |  The elements in the domain that a homomorphism maps to the zero element in the image.";	
Image::usage = " {Int}  |  The set of all values that a homomorphism takes over it's domain.";

QuotientKernels::usage = " {{Int}}  |  The kernels of all quotient groups. ";
NonTrivialQuotientKernels::usage = " {{Int}}  |  The kernels of all quotient groups with trivial elements (group) removed.";

Extend::usage = " Brevity operator for investigating left hand side of second isomorphy theorem. ";
Extending::usage = " Brevity operator for investigating coset multiplication. ";


(* ::Subsection:: *)
(*Properties of single elements in the group*)


(* ::Subsection::Closed:: *)
(*Properties of whole group*)


MakeGroup::usage = " Int,Bool --> Void  |  Constructs \!\(\*SubscriptBox[\(Z\), \(n\)]\) it's subgroups and if second argument is True also all the group's quotient groups are prepared.";

IsNormal::usage = " Since \!\(\*SubscriptBox[\(Z\), \(n\)]\) is cyclic all groups concerned are Abelian. ";


(* ::Subsection:: *)
(*Subsets and their properties*)


(* ::Subsection::Closed:: *)
(*Derived sets and corresponding properites*)


CosetsLeft::usage = " Int,Int[,Int] --> {{Int}}  |  The left cosets of two subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\). If no last argument the whole group is used.";

Coset::usage = " Int,Int,Int --> {Int}  |  The coset of \!\(\*SubscriptBox[\(Z\), \(n\)]\) with representant g.";
Cosets::usage = " Int,Int[,Int] --> {{Int}}  |  The cosets of two subgroups of \!\(\*SubscriptBox[\(Z\), \(n\)]\) if left and right cosets are same, otherwise -1. If no last argument the whole group is used.";

CosetsPowerSet::usage = " The set of all possible cosets in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

CosetByInverse::usage = " Int,Int --> {Int}  |  The coset of a subgroup of \!\(\*SubscriptBox[\(Z\), \(n\)]\) with representant \!\(\*SuperscriptBox[\(g\), \(-\)]\). Subgroup is given by size indexation.";

QuotientGroup::usage = " Int,Int[,Int] --> {Int}  |  The quotient group of two subgroups in \!\(\*SubscriptBox[\(Z\), \(n\)]\). If no last argument the whole group is used.";
QuotientGroups::usage = " {{Int}}  |  All the quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
QuotientGroupGenerators::usage = " {Int}  |  The generators of the quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
QuotientGroupLasts::usage = " {Int}  |  The largest element in each quotient group (n-1) which always generates the whole cyclic group and always is of order n. ";
QuotientGroupSizes::usage = " {Int}  |  Size of the quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\).";
QuotientGroupsRanks::usage = " {Int} | NOT IMPLEMENTED YET. WAITING FOR ELEGANT WAY TO COMPUTE IT. NOT JUST GENERATE IT BY BRUTE FORCE.";

Index::usage = " Int,Int[,Int] --> Int  |  The number of cosets of two subgroup in \!\(\*SubscriptBox[\(Z\), \(n\)]\). If no last argument the whole group is used."
Indexes::usage = " {Int}  |  All indexes in relation to the full group \!\(\*SubscriptBox[\(Z\), \(n\)]\). See Index.";


(* ::Subsection:: *)
(*Isomorphy*)


(* ::Subsection:: *)
(*Structure and graphical overview*)


CosetSizes::usage = " Int --> {{Int}}  |  The size of the cosets in all possible quotient groups of subgrups of Zn.";
CosetNumbers::usage = " Int --> {{Int}}  |  The number of cosets in all possible quotient groups of subgrups of Zn.";

CosetSizeExtremes::usage = " {-1|0|1}  |  Indications of local maximum and minimum values of the coset lengths for all quotient groups.";
CosetSizeMaxima::usage = "  {1|0}  |  Indication of local maximum values of the cosets lengths for all quotient groups.";
CosetSizeMinima::usage = "  {-1|0}  |  Indication of local minimum values of the cosets lengths for all quotient groups.";

QuotientGroupsPartition::usage = " {{Int}}  |  For each index up to N0 a row of the table consists of quotient groups that higher index subgroups takes with it followed " <>
								 " by the quotient groups it takes with lower index subgroups. None is placed on the diagonal.";
QuotientGroupTable::usage = "  {{Int}}  |  The multiplication table of a quotient group.";
 
QuotientPermutationGroups::usage  =  " PermutationGroup  |  The mathematica permutation group isomorphic with this package's representation of \!\(\*SubscriptBox[\(Z\), \(n\)]\).";

QuotientPGroupTables::usage = " Int --> {{GroupMultiplicationTable}}  |  The multiplication tables of all possible quotient groups in \!\(\*SubscriptBox[\(Z\), \(n\)]\)";
QuotientPGroupSizes::usage = " Int -->  {{Int}}  |  The sizes of all possible quotient groups in Zn. ";
QuotientPGroupGenerators::usage = " Int --> {{Cycles}}  |  The set of sets of generators for all quotient groups in Zn.";
QuotientPGroupRanks::usage = " Int --> {{Int}}  |  Matrix of ranks of all the quotient groups in Zn.";
QuotientPGroupRedundancies::usage = " {Real}  |  The quotient between the order (size) of the group and its rank (number of generators). ";
QuotientPGroupReduction::usage = " {Real}  |  The association of rules that map all subgroup orders (sizes) to respective rank (number of generators). ";


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Operators, constants  and primitive sets and mappings*)


CircleDot[gs1_,gs2_]:= DeleteDuplicates@Sort@Flatten[Outer[CirclePlus,gs1,gs2,1],1]
Remove[gs1,gs2];

SmallCircle[gs1_,gs2_]:= Flatten[Map[DeleteDuplicates@*Flatten,Outer[Diamond,gs1,gs2,1],{2}],{1,2}]
Remove[gs1,gs2];
(*
Backslash[H1_,H2_]:= If[Dimensions[H1[[1]]]=={}\[And]Dimensions[H2[[1]]]=={},
							Return[DeleteDuplicates@Sort@Table[Sort[Mod[H1[[i]]+#,N0]& /@ H2],{i,1,Length[H1]}]],
							Return[-1]];
							*)
Backslash[H1_,H2_]:=DeleteDuplicates[Sort/@Table[(H1[[i]]\[CirclePlus]#)&/@H2,{i,1,Length[H1]}]]


Canonical[k_,g_]:= With[{Ns=Sns[[k]]},Return[Sort[(g\[CirclePlus]#)&/@Ns]]]
Remove[k,g,Ns];

Kernel[phi_,k_]:= Module[{ker={},zero=Sns[k],i},
						For[i=1,i<=N0,i++,
							If[MemberQ[phi[k,Zn[[i]]],0],AppendTo[ker,Zn[[i]]],None];];
						Return[ker];]
						
Unprotect[Image];Image[phi_]:= Module[{},Return[Table[phi[k,#]&/@Zn,{k,1,N1}]];]
Remove[phi,k,ker,zero];

NonTrivialQuotientKernels[]:= Table[Select[Sns[[i]],(First[Canonical[j,#]]== 0)&],{i,1,N1-1},{j,2,N1-1}];
QuotientKernels[]:= Table[Select[Sns[[i]],(First[Canonical[j,#]]== 0)&],{i,1,N1},{j,1,N1}];

Extend[H1_,H2_]:=(H1\[CircleDot]H2)\[Backslash]H2;
Extending[H1_,g_]:=(g\[CirclePlus]#)& /@ H1;
Remove[H1,H2,Ns,g];


(* ::Subsection:: *)
(*Group criterias*)


(* ::Subsection:: *)
(*Undependent theory*)


(* ::Subsubsection:: *)
(*Elementwise*)


MakeGroup[n_,quotientgroups_]:= Module[{t,runtime},
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
												Css=CosetsPowerSet[];
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
					
Remove[t,runtime,quotientgroups,n];

Coset[k_,r_]:= SelectFirst[Cosets[k],MemberQ[#,r]&]
Remove[k,r];

CosetByInverse[k_,g_]:= Sort[(SuperMinus[g]\[CirclePlus]#)& /@ Sns[[k]]]
Remove[k,g];



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


(* ::Subsubsection:: *)
(*Derived sets*)


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

IsNormal = True; (* all cyclic groups are normal *)

Cosets[k_,l_]:= If[IsNormal,Return[CosetsLeft[k,l]],Print["not commutative group"];Return[-1]]
Cosets[k_]:= Cosets[1,k]
Remove[k,l];

CosetsPowerSet[]:= If[Css=={},Table[Cosets[i,j],{i,1,N1},{j,i,N1}],Return[Css]]

CosetSizes[]:= Length \[Congruent] CosetsPowerSet[]
CosetNumbers[]:= (N0/Length[#])& \[Congruent] CosetsPowerSet[]

CosetSizeMinima[]:= With[{L=CosetSizes[]}, 
									Table[
										If[AdditiveGroupQuotients`Private`dipp[#],-1,0,0]& /@ Partition[L[[i]],3,1],
											{i,1,N1-2}]]
CosetSizeMaxima[]:= With[{L=CosetSizes[]}, 
									Table[
										If[AdditiveGroupQuotients`Private`peak[#],1,0,0]& /@ Partition[L[[i]],3,1],
											{i,1,N1-2}]]
CosetSizeExtremes[]:= CosetSizeMaxima[] + CosetSizeMinima[]
Remove[L];
											
Index[k_,l_]:= If[IsNormal,Return[Length[CosetsLeft[k,l]]],Print["not commutative group"];Return[-1]]
Index[k_]:= Index[1,k];
Remove[k,l];

Indexes[]:= Index[1,#]&/@Range[N1];

QuotientGroup[k_,l_]:= First /@ AdditiveGroupQuotients`Private`QuotientGroupFull[k,l]
 
QuotientGroups[]:= Map[First,AdditiveGroupQuotients`Private`QuotientGroupsFull[],{3}]
																		
QuotientGroupTable[n_,m_] := With[{H=QuotientGroup[n,m]},
								Module[{r=N0,mtable},
										N0=H[[-1]]+(H[[-1]]-H[[-2]]);
										mtable=Table[H[[i]]\[CirclePlus]H[[j]],{i,1,Length[H]},{j,1,Length[H]}];
										N0=r;
										Return[mtable];
									]];
Remove[n,m,H,r,mtable];

QuotientGroupGenerators[]:= (#[[2]])& \[Congruent](Rest /@ QuotientGroups[][[1;;-2]])							
QuotientGroupLasts[]:= (#[[-1]])& \[Congruent](QuotientGroups[][[1;;-2]])		
QuotientGroupSizes[]:= Length \[Congruent] QuotientGroups[]
QuotientGroupsRanks[] := Print["not implemented yet"]					
																											
QuotientPermutationGroups[]:= Table[PermutationGroup[CyclesMap[] /@ (QuotientGroup[i,j])],{i,1,N1},{j,i,N1}]
QuotientGroupsPartition[]:= Table[ {(QuotientGroup[#,i]& /@ Range[1,i]),
									 None,
									(QuotientGroup[i,#]& /@ Range[i,N1])},{i,1,N1}]

QuotientPGroupGenerators[]:= GroupGenerators \[Congruent] QuotientPermutationGroups[]
QuotientPGroupSizes[]:= GroupOrder \[Congruent] QuotientPermutationGroups[]
QuotientPGroupRanks[]:= Length \[Congruent] QuotientPGroupGenerators[]
QuotientPGroupTables[]:= (MatrixForm@*GroupMultiplicationTable) \[Congruent] QuotientPermutationGroups[]
QuotientPGroupReduction[]:={GroupOrder[#]->Length[GroupGenerators[#]]}& \[Congruent] QuotientPermutationGroups[]
QuotientPGroupRedundancies[]:=(GroupOrder[#]/Length[GroupGenerators[#]])& \[Congruent] QuotientPermutationGroups[]


(* ::Subsection:: *)
(*Helpers*)


Begin["`Private`"];

	dipp::usage = " Int,Int,Int  --> True|False  |  Returns true if the middle value is strictly the smallest of the three given.";
	peak::usage = " Int,Int,Int  --> True|False  |  Returns true if the middle value is strictly the largest of the three given.";
	
	dipp[xs_]:=(xs[[2]]<xs[[1]]&&xs[[2]]<xs[[3]])
	peak[xs_]:=(xs[[2]]>xs[[1]]&&xs[[2]]>xs[[3]])
	
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
	Remove[n,m,mtable];

	
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
