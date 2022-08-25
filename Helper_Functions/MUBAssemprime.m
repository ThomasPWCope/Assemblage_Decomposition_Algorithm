function[assemblage]=MUBAssemprime(p)
state=MaxEntangled(p)*transpose(MaxEntangled(p));
assemblages=zeros(p,p,p,p+1);
barr=MUBPrime(p);
for measurement=1:p+1
for outcome=1:p
assemblage(:,:,outcome,measurement)=PartialTrace(kron(conj(squeeze(barr(measurement,outcome,:)))*transpose(squeeze(barr(measurement,outcome,:))),eye(p))*state,1,p);
end
end



