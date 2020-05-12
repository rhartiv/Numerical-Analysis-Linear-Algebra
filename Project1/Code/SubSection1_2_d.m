%SubSection1_2_d
% 1.2 Script

% Matrix construction
tic
a=-1;
k=1:
n=10.^k;
for i=1:length(n)
    Mnew{i}=eye(n(i))+diag(a*ones(1,n(i)-1),1);
    bnew{i}=zeros(n(i),1);
    bnew{i}(1)=1;
end

%Matrix properties
    for i=1:length(n)
        % Running the functions and plotting values of omega for matrices

        [xJnew{i},kJnew{i}]=Jacobi(Mnew{i},bnew{i});
        [xGSnew{i},kGSnew{i}]=GaussSeidel(Mnew{i},bnew{i});
        [xSDnew{i},kSDnew{i}]=SteepestDescent(Mnew{i},bnew{i});
        [xCGnew{i},kCGnew{i}]=ConjugateGradient(Mnew{i},bnew{i});
        
        lnew{i}=length(Mnew{i});
        
         % Finding, numerically, the optimal omega
        for ii=1:200
            [txnew{i},tknew{i}]=SOR(Mnew{i},bnew{i},ii/100);
            knew{i}(ii+1)=tknew{i};
        end
        knew{i}(1)=[]
        optomegaiterationsnew{i}=min(knew{i});
        oolocnew{i}=find(knew{i}==optomegaiterationsnew{i});
        inew{i}=[0.02:0.01:2.01];

        figure(i)
        plot(inew{i},knew{i})
        hold on
        plot(inew{i}(oolocnew{i}),knew{i}(oolocnew{i}), '-p')
        title(['Iteration vs. Omega for matrix of size ' num2str(lnew{i}, '%.0f')])
        if i==1
            xlabel('Omega value')
            ylabel('Number of iterations')
        else
            set(gca, 'YTick', [])
        end
        
        % Running the SOR function with the numerically calculated optimum
        % omega
        [xSORnew{i},kSORnew{i}]=SOR(Mnew{i},bnew{i},mean(oolocnew{i}/100))
    end
    for i=1:length(n)
        lnew{i}=length(Mnew{i});
        eMnew{i}=eig(Mnew{i});
        % Condition Number
        CNMnew{i}=['The (2) condition number of the matrix of size ', num2str(lnew{i}, '%.0f'), ' is ', num2str(cond(Mnew{i}),'%.0f')];
        disp(CNMnew{i})

        % Spectral Radius
        Dnew{i}=diag(Mnew{i}).*eye(length(Mnew{i}));
        Unew{i}=-triu(Mnew{i})+Dnew{i};
        Lnew{i}=-tril(Mnew{i})+Dnew{i};
        Tjnew{i}=inv(Dnew{i})*(Lnew{i}+Unew{i});
        Tgnew{i}=inv(Dnew{i}-Lnew{i})*Unew{i};
        Tsornew{i}=inv(Dnew{i}-mean(oolocnew{i}/100)*Lnew{i})*(mean(oolocnew{i}/100)*Unew{i}+(1-mean(oolocnew{i}/100))*Dnew{i});
        Jradnew{i}=max(abs(eig(Tjnew{i})));
        GSradnew{i}=max(abs(eig(Tgnew{i})));
        SORradnew{i}=max(abs(eig(Tsornew{i})));

        SRMnew{i}=['The spectral radii of the matrix of size', num2str(lnew{i}), ' are ', num2str(Jradnew{i},'%.3f'), ' for jacobi, ', num2str(GSradnew{i},'%.3f'), ' for Gauss-Seidel, and ', num2str(SORradnew{i},'%.3f'), ' for SOR.'];
        disp(SRMnew{i})

        % Test if symmetric
        if issymmetric(Mnew{i})==true
            symm{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is Symmetric'];
            disp(symm{i})
        else
            symm{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is not Symmetric'];
            disp(symm{i})
        end

        % Test if positive definite
        pd(i)=1; 
        for ii=1:lnew{i}
            if eMnew{i}(ii)<=0
                pd(i)=0;
            end
        end
        if pd(i)==1
            posdef{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is Positive definite'];
            disp(posdef{i})
        else
            posdef{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is not Positive definite'];
            disp(posdef{i})
        end

        % Test if tridiagonal
        if isbanded(Mnew{i},1,1)==true
            tridiag{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is tridiagonal, with diagonal elements greater than zero but the rest less than zero'];
            disp(tridiag{i})
        else
            tridiag{i}=['The matrix of size ', num2str(lnew{i}, '%.0f'), ' is not tridiagonal']
            disp(tridiag{i})
        end

        % Displaying optimal omega (from formula and from numerical generation)
        allOneStringnew{i} = sprintf('%.2f,' , inew{i}(oolocnew{i}));
        if Jradnew{i}<1
            optomegaformulanew{i}=2/(1+sqrt(1-(max(abs(Jradnew{i})))^2));
            optomeganew{i}=['The optimal numerically calculated value(s) of omega is/are ', allOneStringnew{i}, ' and the formula for the optimal omega yields ', num2str(optomegaformulanew{i},'%.3f')];
        else
            optomeganew{i}=['The optimal numerically calculated value(s) of omega is/are ', allOneStringnew{i}, ' and because the spectral radius of the Jacobi T matrix is greater than or equal to 1, the formula for the optimal omega is imaginary'];
        end
        disp(optomeganew{i})

        % Test if Conjugate Gradient converges in n iterations
        if kCGnew{i}<=lnew{i}
            CG{i}=['Conjugate gradient of the matrix of size ', num2str(lnew{i}, '%.0f'), ' converges in n iterations'];
            disp(CG{i})
        else
            CG{i}=['Conjugate gradient of the matrix of size ', num2str(lnew{i}, '%.0f'), ' does not converge in n iterations'];
            disp(CG{i})
        end
    end
toc