%QR Method
function [Vals,Vecs,k]=QR(A)

kmax=250
copyA=A
[Q,R]=qr(copyA)
k=1
while k<kmax
    if istriu(round(Q'*copyA*Q*10000))==1
        valmatrix=Q'*copyA*Q
        for i=1:1:length(valmatrix)
            Vals(i)=valmatrix(i,i)
        end
        break
    elseif k==kmax-1
        valmatrix=Q'*copyA*Q
        for i=1:1:length(valmatrix)
            Vals(i)=valmatrix(i,i)
        end
        break
    end
    copyA=Q'*copyA*Q
    [Q,R]=qr(copyA)
    k=k+1
end
