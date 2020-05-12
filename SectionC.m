% 1.1 c script

M4=readmatrix('matrix4.txt');
M14=readmatrix('matrix14.txt');

    l4=length(M4);
    b4=randi(4,l4,1);
    eM4=eig(M4);
    
    % Running the function and plotting values of omega for matrix 4
    
    [xJ4,kJ4]=Jacobi(M4,b4);
    [xGS4,kGS4]=GaussSeidel(M4,b4);
    [xSD4,kSD4]=SteepestDescent(M4,b4);
    [xCG4,kCG4]=ConjugateGradient(M4,b4);
    
     % Finding the optimal omega
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
    plot(i4(ooloc4),k4(ooloc4), '-p')
    title('Number of iterations for convergence of matrix 4 using varied values of omega in SOR')
    xlabel('Omega value')
    ylabel('Number of iterations')
    
    
    l14=length(M14);
    b14=randi(14,l14,1);
    eM14=eig(M14);
    
    % Running the function and plotting values of omega for matrix 14
    
    [xJ14,kJ14]=Jacobi(M14,b14);
    [xGS14,kGS14]=GaussSeidel(M14,b14);
    [xSD14,kSD14]=SteepestDescent(M14,b14);
    [xCG14,kCG14]=ConjugateGradient(M14,b14);
    
     % Finding the optimal omega
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
    plot(i4(ooloc14),k4(ooloc14), '-p')
    title('Number of iterations for convergence of matrix 14 using varied values of omega in SOR')
    xlabel('Omega value')
    ylabel('Number of iterations')

disp('\linebreak')

% Matrix 4
    
    % Condition Number
    CNM4=['The (2) condition number of matrix 4 is ', num2str(cond(M4),'%.0f')];
    disp(CNM4)
    
    % Spectral Radius
    SRM4=['The spectral radius of matrix 4 is ', num2str((max(abs(eM4))),'%.3f')];
    disp(SRM4)
    
    % Test if symmetric
    if issymmetric(M4)==true
        disp('Matrix 4 is Symmetric')
    else
        disp('Matrix 4 is not Symmetric')
    end

    % Test if positive definite
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
        disp('Matrix 4 is tridiagonal')
        for i=1:l4
            for j=1:l4~=i
                if M4(i,i)>0
                    if M4(i,j)<0
                        disp(', with diagonal elements greater than zero but the rest less than zero')
                    end
                end
            end
        end
    else
        disp('Matrix 4 is not tridiagonal')
    end
    
    % Displaying optimal omega
    allOneString4 = sprintf('%.2f,' , i4(ooloc4));
    optomega4=['The optimal value(s) of omega is/are ', allOneString4];
    disp(optomega4)
    
    % Test if Conjugate Gradient converges in n iterations
    if kCG4<=l4
        disp('Conjugate gradient of matrix 4 converges in n iterations')
    else
        disp('Conjugate gradient of matrix 4 does not converge in n iterations')
    end 
    
% Matrix 14
    
    % Condition Number
    CNM14=['The (2) condition number of matrix 14 is ', num2str(cond(M14),'%.0f')];
    disp(CNM14)
    
    % Spectral Radius
    SRM14=['The spectral radius of matrix 14 is ', num2str((max(abs(eM14))),'%.3f')];
    disp(SRM14)
    
    % Test if symmetric
    if issymmetric(M14)==true
        disp('Matrix 14 is Symmetric')
    else
        disp('Matrix 14 is not Symmetric')
    end

    % Test if positive definite
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
        disp('Matrix 14 is tridiagonal')
        for i=1:l14
            for j=1:l14~=i
                if M14(i,i)>0
                    if M14(i,j)<0
                        disp(', with diagonal elements greater than zero but the rest less than zero')
                    end
                end
            end
        end
    else
        disp('Matrix 14 is not tridiagonal')
    end
    
    % Displaying optimal omega
    allOneString14 = sprintf('%.2f,' , i14(ooloc14));
    optomega14=['The optimal value(s) of omega is/are ', allOneString14];
    disp(optomega14)
    
    % Test if Conjugate Gradient converges in n iterations
    if kCG14<=l14
        disp('Conjugate gradient of matrix 14 converges in n iterations')
    else
        disp('Conjugate gradient of matrix 14 does not converge in n iterations')
    end 
    

