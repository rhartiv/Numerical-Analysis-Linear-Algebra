function x=DiagonallyDominant(A)
n=length(A)
l=1
for i=1:n
    k=0
    for j=1:n
        if j~=i
            k=k+abs(A(i,j))
        end
        if k>abs(A(i,i))
            l=false
        end
    end
end
if l==1
    x='This matrix is diagonally dominant'
else
    x='This matrix is not diagonally dominant'
end
end
        