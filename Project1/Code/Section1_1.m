% 1.1 c script
tic
M4=readmatrix('matrix4.txt');
M14=readmatrix('matrix14.txt');

    l4=length(M4);
    b4=randi(4,l4,1);
    
    % Running the functions and plotting values of omega for matrix 4
    
    [xJ4,kJ4]=Jacobi(M4,b4);
    [xGS4,kGS4]=GaussSeidel(M4,b4);
    [xSD4,kSD4]=SteepestDescent(M4,b4);
    [xCG4,kCG4]=ConjugateGradient(M4,b4);
    
     % Finding, numerically, the optimal omega for matrix 4
    for i=1:200
        [tx4,tk4]=SOR(M4,b4,i/100);
        k4(i+1)=tk4;
    end
    k4(1)=[];
    optomegaiterations4=min(k4);
    ooloc4=find(k4==optomegaiterations4);
    i4=[0.01:0.01:2];
    
    figure(1)
    plot(i4,k4)
    hold on
    title(['Number of iterations for convergence of matrix 4 using varied values of omega in SOR'])
    plot(i4(ooloc4),k4(ooloc4), '-p')
    xlabel('Omega value')
    ylabel('Number of iterations')
    
    [xSOR4,kSOR4]=SOR(M4,b4,mean(ooloc4/100));
   
    l14=length(M14);
    b14=randi(14,l14,1);
    
    % Running the functions and plotting values of omega for matrix 14
    
    [xJ14,kJ14]=Jacobi(M14,b14);
    [xGS14,kGS14]=GaussSeidel(M14,b14);
    [xSD14,kSD14]=SteepestDescent(M14,b14);
    [xCG14,kCG14]=ConjugateGradient(M14,b14);
    
     % Finding, numerically, the optimal omega for matrix 14
    for i=1:200
        [tx14,tk14]=SOR(M14,b14,i/100);
        k14(i+1)=tk14;
    end
    k14(1)=[];
    optomegaiterations14=min(k14);
    ooloc14=find(k14==optomegaiterations14);
    i14=[0.01:0.01:2];
    
    figure(2)
    plot(i14,k14)
    hold on
    plot(i14(ooloc14),k14(ooloc14), '-p')
    title(['Number of iterations for convergence of matrix 14 using varied values of omega in SOR'])
    xlabel('Omega value')
    ylabel('Number of iterations')

    [xSOR14,kSOR14]=SOR(M14,b14,mean(ooloc14/100));
    
% Matrix 4 properties
    
    % Condition Number
    CNM4=['The (2) condition number of matrix 4 is ', num2str(cond(M4),'%.0f')];
    disp(CNM4)
    
    % Spectral Radii
    D4=diag(M4).*eye(length(M4));
    U4=-triu(M4)+D4;
    L4=-tril(M4)+D4;
    Tj4=inv(D4)*(L4+U4);
    Tg4=inv(D4-L4)*U4;
    Tsor4=inv(D4-mean(ooloc4/100)*L4)*(mean(ooloc4/100)*U4+(1-mean(ooloc4/100))*D4);
    Jrad4=max(abs(eig(Tj4)));
    GSrad4=max(abs(eig(Tg4)));
    SORrad4=max(abs(eig(Tsor4)));
    
    SRM4=['The spectral radii of matrix 4 are ', num2str(Jrad4,'%.3f'), ' for jacobi, ', num2str(GSrad4,'%.3f'), ' for Gauss-Seidel, and ', num2str(SORrad4,'%.3f'), ' for SOR.'];
    disp(SRM4)
    
    % Test if symmetric
    if issymmetric(M4)==true
        disp('Matrix 4 is Symmetric')
    else
        disp('Matrix 4 is not Symmetric')
    end

    % Test if positive definite
    eM4=eig(M4);
    pd=1; 
    for i=1:l4
        if eM4(i)<=0
            pd=0;
        end
    end
    if pd==1
        disp('Matrix 4 is positive definite')
    else
        disp('Matrix 4 is not positive definite')
    end
    
    % Test if tridiagonal
    if isbanded(M4,1,1)==true
        disp('Matrix 4 is tridiagonal, with diagonal elements greater than zero but the rest less than zero')
    else
        disp('Matrix 4 is not tridiagonal')
    end
    
    % Displaying optimal omega (from formula and from numerical generation)
    allOneString4 = sprintf('%.2f,' , i4(ooloc4));
    if Jrad4<1
        optomegaformula4=2/(1+sqrt(1-(max(abs(Jrad4)))^2));
        optomega4=['The optimal numerically calculated value(s) of omega is/are ', allOneString4, ' and the formula for the optimal omega yields ', optomegaformula4];
    else
        optomega4=['The optimal numerically calculated value(s) of omega is/are ', allOneString4, ' and because the spectral radius of the Jacobi T matrix is greater than or equal to 1, the formula for the optimal omega is imaginary'];
    end
    disp(optomega4)
    
    % Test if Conjugate Gradient converges in n iterations
    if kCG4<=l4
        disp('Conjugate gradient of matrix 4 converges in n iterations')
    else
        disp('Conjugate gradient of matrix 4 does not converge in n iterations')
    end 
    
% Matrix 14 properties
    
    % Condition Number
    CNM14=['The (2) condition number of matrix 14 is ', num2str(cond(M14),'%.0f')];
    disp(CNM14)
    
    % Spectral Radii
    D14=diag(M14).*eye(length(M14));
    U14=-triu(M14)+D14;
    L14=-tril(M14)+D14;
    Tj14=inv(D14)*(L14+U14);
    Tg14=inv(D14-L14)*U14;
    Tsor14=inv(D14-mean(ooloc14/100)*L14)*(mean(ooloc14/100)*U14+(1-mean(ooloc14/100))*D14);
    Jrad14=max(abs(eig(Tj14)));
    GSrad14=max(abs(eig(Tg14)));
    SORrad14=max(abs(eig(Tsor14)));
    
    SRM14=['The spectral radii of matrix 14 are ', num2str(Jrad14,'%.3f'), ' for jacobi, ', num2str(GSrad14,'%.3f'), ' for Gauss-Seidel, and ', num2str(SORrad14,'%.3f'), ' for SOR.'];
    disp(SRM14)
    
    % Test if symmetric
    if issymmetric(M14)==true
        disp('Matrix 14 is Symmetric')
    else
        disp('Matrix 14 is not Symmetric')
    end

    % Test if positive definite
    eM14=eig(M14);
    pd=1; 
    for i=1:l14
        if eM14(i)<=0
            pd=0;
        end
    end
    if pd==1
        disp('Matrix 14 is positive definite')
    else
        disp('Matrix 14 is not positive definite')
    end
    
    % Test if tridiagonal
    if isbanded(M14,1,1)==true
        disp('Matrix 14 is tridiagonal, with diagonal elements greater than zero but the rest less than zero')
    else
        disp('Matrix 14 is not tridiagonal')
    end
    
    % Displaying optimal omega (from formula and from numerical generation)
    allOneString14 = sprintf('%.2f,' , i14(ooloc14));
    if Jrad14<1
        optomegaformula14=2/(1+sqrt(1-Jrad14^2));
        optomega14=['The optimal numerically calculated value(s) of omega is/are ', allOneString14, ' and the formula for the optimal omega yields ', num2str(optomegaformula14,'%.3f')];
    else
        optomega14=['The optimal numerically calculated value(s) of omega is/are ', allOneString14, ' and because the spectral radius of the Jacobi T matrix is greater than or equal to 1, the formula for the optimal omega is imaginary'];
    end
    disp(optomega14)
    
    % Test if Conjugate Gradient converges in n iterations
    if kCG14<=l14
        disp('Conjugate gradient of matrix 14 converges in n iterations')
    else
        disp('Conjugate gradient of matrix 14 does not converge in n iterations')
    end 
toc    

