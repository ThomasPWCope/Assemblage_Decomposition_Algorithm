%This is Example 3 in the main paper.

ExtremalAssemTest5=MUBAssemprime(3);
%Note however, that the above assemblage is in "Skypryk" form.
ExtremalAssemTest5=permute(ExtremalAssemTest5,[3,4,1,2]);
ExtremalAssemTest6=0.4*ExtremalAssemTest5+0.6*identity_assemblage_creator(4,3,3); 
ExtremalReduced=ExtremalAssemTest6(:,1:2,:,:);
size(ExtremalReduced)
OutCell6h=zero_recursion_H([],ExtremalReduced,1,2,10^-6);


ReducedExtremal6=congregate_extremals(OutCell6h);

EqCheck=ReducedExtremal6;
EqCheck{3728,1}=1;
EqCheck{3728,2}=ExtremalAssemTest5(:,1:2,:,:);
congregate_extremals(EqCheck,10^-3)
