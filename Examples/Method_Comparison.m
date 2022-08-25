% In this notebook I will compare three different versions of the algorithm
% all with varying levels of protection from small perturbations problems. 

%The algorithms will be compared on the following test assemblages:

%2 trivial assemblages with 1 measurement.
%Mutually Unbiased Basis for d=2,3 (and maybe 4).
% Noisy mutually unbiased basis (i.e. mixed with identity matrix)


%They will be compared for:
%1) whether the decomposition gives again the original input.
%2) Whether each decomposition element is a proper assemblage.
%3) Correct answer: (where applicable)

ResultsCell=cell(8,3);

%Test 1: a trivial assemblage with two extremals.
ExtremalAssemTest1=zeros(2,1,2,2);
ExtremalAssemTest1(1,1,:,:)=[[1,0];[0,0]]/2;
ExtremalAssemTest1(2,1,:,:)=[[0,0];[0,1]]/2;
OutCell1o=zero_recursion_original([],ExtremalAssemTest1,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest1,OutCell1o);
OutCell1h=zero_recursion_H([],ExtremalAssemTest1,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest1,OutCell1h);
OutCell1r=zero_recursion_rank([],ExtremalAssemTest1,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest1,OutCell1r);

ResultsCell{1,1}=OutCell1o;
ResultsCell{1,2}=OutCell1h;
ResultsCell{1,3}=OutCell1r;
save('DecompAnalysis','ResultsCell')

%So far, the Original and Hermitian perform best.
ExtremalAssemTest2=zeros(2,1,2,2);
ExtremalAssemTest2(1,1,:,:)=eye(2)/8 +[[1/2,1/2];[1/2,1/2]]/4;
ExtremalAssemTest2(2,1,:,:)=eye(2)/8 +[[1/2,-1/2];[-1/2,1/2]]/4;
OutCell2o=zero_recursion_original([],ExtremalAssemTest2,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest2,OutCell2o);
OutCell2h=zero_recursion_H([],ExtremalAssemTest2,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest2,OutCell2h);
OutCell2r=zero_recursion_rank([],ExtremalAssemTest2,1,1,10^-6);
assemblage_analysis(ExtremalAssemTest2,OutCell2r);

ResultsCell{2,1}=OutCell2o;
ResultsCell{2,2}=OutCell2h;
ResultsCell{2,3}=OutCell2r;

%Now we start to see trade offs- but the Hermitian seems to do best so far.
ExtremalAssemTest3=MUBAssemprime(2);
%Note however, that the above assemblage is in "Skypryk" form.
ExtremalAssemTest3=permute(ExtremalAssemTest3,[3,4,1,2]);
OutCell3o=zero_recursion_original([],ExtremalAssemTest3,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest3,OutCell3o);
OutCell3h=zero_recursion_H([],ExtremalAssemTest3,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest3,OutCell3h);
OutCell3r=zero_recursion_rank([],ExtremalAssemTest3,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest3,OutCell3r);

ResultsCell{3,1}=OutCell3o;
ResultsCell{3,2}=OutCell3h;
ResultsCell{3,3}=OutCell3r;

IdAssem2=0*ExtremalAssemTest3;
ExtremalAssemTest4=0.8*ExtremalAssemTest3+0.1*identity_assemblage_creator(3,2,2);
OutCell4o=zero_recursion_original([],ExtremalAssemTest4,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest4,OutCell4o);
% Now we start to see problems during the original function: However, the
% final answer tallies with what we expect.
OutCell4h=zero_recursion_H([],ExtremalAssemTest4,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest4,OutCell4h);
OutCell4r=zero_recursion_rank([],ExtremalAssemTest4,1,3,10^-6);
assemblage_analysis(ExtremalAssemTest4,OutCell4r);


ResultsCell{4,1}=OutCell4o;
ResultsCell{4,2}=OutCell4h;
ResultsCell{4,3}=OutCell4r;

save('DecompAnalysis','ResultsCell')

ExtremalAssemTest5=MUBAssemprime(3);
%Note however, that the above assemblage is in "Skypryk" form.
ExtremalAssemTest5=permute(ExtremalAssemTest5,[3,4,1,2]);
tic
OutCell5o=zero_recursion_original([],ExtremalAssemTest5,1,4,10^-6);
toc
assemblage_analysis(ExtremalAssemTest5,OutCell5o);
tic
OutCell5h=zero_recursion_H([],ExtremalAssemTest5,1,4,10^-6);
toc
assemblage_analysis(ExtremalAssemTest5,OutCell5h);
tic
OutCell5r=zero_recursion_rank([],ExtremalAssemTest5,1,4,10^-6);
toc
assemblage_analysis(ExtremalAssemTest5,OutCell5r);

ResultsCell{5,1}=OutCell5o;
ResultsCell{5,2}=OutCell5h;
ResultsCell{5,3}=OutCell5r;

save('DecompAnalysis','ResultsCell')

ExtremalAssemTest6=0.4*ExtremalAssemTest5+0.6*identity_assemblage_creator(4,3,3);
tic
OutCell6o=zero_recursion_original([],ExtremalAssemTest6,1,4,10^-6);
toc
%The above function fails: it cannot do a perturbation.

assemblage_analysis(ExtremalAssemTest6,OutCell6o);
tic
OutCell6h=zero_recursion_H([],ExtremalAssemTest6,1,4,10^-6);
toc
assemblage_analysis(ExtremalAssemTest6,OutCell6h);

tic
OutCell6r=zero_recursion_rank([],ExtremalAssemTest6,1,4,10^-6);
toc
assemblage_analysis(ExtremalAssemTest6,OutCell6r);

ResultsCell{6,1}=OutCell6o;
ResultsCell{6,2}=OutCell6h;
ResultsCell{6,3}=OutCell6r;

save('DecompAnalysis','ResultsCell')

ExtremalAssemTest7=ExtremalAssemTest3;
ExtremalAssemTest7(1,3,:,:)=eye(2)/2;
ExtremalAssemTest7(2,3,:,:)=0;
tic
OutCell7o=zero_recursion_original([],ExtremalAssemTest7,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest7,OutCell7o);
tic
OutCell7h=zero_recursion_H([],ExtremalAssemTest7,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest7,OutCell7h);
tic
OutCell7r=zero_recursion_rank([],ExtremalAssemTest7,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest7,OutCell7r);


ResultsCell{7,1}=OutCell7o;
ResultsCell{7,2}=OutCell7h;
ResultsCell{7,3}=OutCell7r;

save('DecompAnalysis','ResultsCell')

ExtremalTestAssem8=0.91*ExtremalAssemTest7+0.09*identity_assemblage_creator(3,2,2);
tic
OutCell8o=zero_recursion_original([],ExtremalAssemTest8,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest8,OutCell8o);
tic
OutCell8h=zero_recursion_H([],ExtremalAssemTest8,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest8,OutCell8h);
tic
OutCell8r=zero_recursion_rank([],ExtremalAssemTest8,1,3,10^-6);
toc
assemblage_analysis(ExtremalAssemTest7,OutCell8r);

ResultsCell{8,1}=OutCell8o;
ResultsCell{8,2}=OutCell8h;
ResultsCell{8,3}=OutCell8r;

save('DecompAnalysis','ResultsCell')








