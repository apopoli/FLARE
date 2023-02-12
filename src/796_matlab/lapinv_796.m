function [tarray,fzinv,ifzeval,ifail] = lapinv_796(fz,tarray,tol,ssbar,nmax,sigma0)
%lapinv_796 FOURIER BASED METHOD FOR THE NUMERICAL INVERSION OF LAPLACE TRANSFORMS 
% fz: handle to function in complex domain, to be antitrasformed to f(t)
% tarray=linspace(0,2*pi,250);
% tol = 10^-5;
% ssbar = 0;
% nmax = 550;
% sigma0 = 0.;
%
% EXAMPLE
% fz = @(z) 1/(1^2+z^2);
% tarray = linspace(0,4*pi,250);
% tol = 10^-5; ssbar = 0; nmax = 550; sigma0 = 0.;
% [tarray,fzinv,~,~] = lapinv_796(fz,tarray,tol,ssbar,nmax,sigma0);

%   FOURIER BASED METHOD FOR THE NUMERICAL INVERSION

%                  OF LAPLACE TRANSFORMS


%   AUTHORS: D'AMORE LUISA, GIULIANO LACCETTI, ALMERICO MURLI

%   DRIVER  PROGRAM TO TEST  THE ROUTINE INVLTF
%   FOR THE INVERSION OF A LAPLACE TRANSFORM FUNCTION.
%   THIS VERSION USES BOTH REAL AND COMPLEX DOUBLE PRECISION OPERATIONS.
%   THE MAIN PROGRAM ALLOWS THE INVERSION OF A SET OF LAPLACE TRANSFORM
%   FUNCTIONS WHICH CAN BE ADDRESSED  BY A NATURAL NUMBER  BETWEEN 1 AND 34.
%   FOR THE COMPLETE LIST SEE THE TABLE IN THE COMPANION PAPER.

%   THIS IS A SELF-CONTAINED DRIVER FOR THE INVLTF ROUTINE
%   WITH SUBPROGRAMS QDACC AND THE LAPLACE TRANSFORM PAIRS FZ AND FEX

%   THIS DRIVER ALLOWS TO OBTAIN UP TO 44 VALUES OF THE
%   INVERSE LAPLACE FUNCTION

% DEFAULT INPUT
%nmax = ...               % ON ENTRY, NMAX SPECIFIES THE MAXIMUM 
                           % NUMBER OF EVALUATIONS OF FZ ALLOWED.
%sigma0 = ...           % ON ENTRY, SIGMA0 CONTAINS THE VALUE
                           % OF THE ABSCISSA OF CONVERGENCE OF
%CASE (18)                     % THE LAPLACE TRANSFORM FUNCTION TO BE
%sigma0=3.;                    % INVERTED OR AN UPPER BOUND TO THIS.
                           % IT IS RECOMMENDED THAT A CORRECT VALUE TO
%CASE (23)                     % SIGMA0 OR A CORRECT UPPER BOUND TO
%sigma0=0.25;                  % SIGMA0 IS PROVIDED. IF AN INCORRECT VALUE
                           % IS USED THE ROUTINE APPEARS TO WORK WELL
%CASE (29)                     % BUT CONVERGES TO COMPLETELY WRONG RESULTS.
%sigma0=2;                     % THERE IS NO WAY IN WHICH THE ROUTINE CAN
                           % DETECT THIS.
                           % NOTE: CHANGE THE VALUE OF SIGMA FOR CASES
                           % 18, 23 AND 29.
                           
                                               
                           
nt=size(tarray,2);       % IS THE NUMBER OF T-VALUES USED HERE FOR THE TEST  FUNCTIONS.
                         % THE MAXIMUM ALLOWED IS THE DIMENSION OF THE ARRAY TVAL

%    *****************************************************************
%    INITIALIZATION
%    *****************************************************************

difabs=zeros(nt,1);
difrel=zeros(nt,1);
fzinv=zeros(nt,1);
error=zeros(3,nt);
ERR=zeros(3,1);
part0_b=0;
part1_b=0;
part2_b=0;

%     *****************************************************************
%     IN THIS TEST PROGRAM THE INVERSE FUNCTION IS REQUESTED IN THE
%     FOLLOWING VALUES OF T:
%     T=1,20, STEP=0.5 AND T=20,65 STEP=5
%     REMARK:
%     ======
%     FOR T SMALL, THAT IS FOR  T=1,20, STEP=0.5, ALL THE RESULTS ARE
%     QUITE ACCURATE, WHILE FOR T LARGE, THAT IS FOR T=20,65 STEP=5,
%     IN SOME CASES THE RESULTS COULDN'T BE ACCURATE. IN SUCH CASES
%     THE AUTHORS SUGGEST TO USE AN ASYMPTOTIC INVERSION METHOD.
%     *****************************************************************

% ************************************************************
%                    CALL OF THE ROUTINE INVLTF
% ************************************************************

fprintf('Numerical inversion \n')
for  i=1:nt
    %chiamata a antitrasformata analitica per calcolare accuratezza soluzione
    % exf(i) = fex(tarray(i)); 
    
    [fzinv(i), ERR,ifzeval(i),ifail(i)] = invltf(tol, tarray(i), fz, sigma0, ssbar, nmax);

    error(1,i) = ERR(1);
    error(2,i) = ERR(2);
    error(3,i) = ERR(3);

    fprintf('it %i/%i; t=%g s; Ip=%g A; ifzeval=%i; ifail=%i \n',i,nt,tarray(i),fzinv(i),ifzeval(i),ifail(i));
end

% CALCOLO DELLA ACCURATEZZA
% TOGLIERE COMMENTI SE NOTA exf, SOLUZIONE NEL TEMPO

% ****************************************************************************

%     ESTIMATION OF THE ERRORS

% ****************************************************************************

% for  i = 1:nt
% 
% if(ifail(i) == 0 ) 
% 
%   difabs(i) = abs(exf(i) - fzinv(i));
%   if (exf(i) == 0)
%        difrel(i) = difabs(i);
%     
%   else
%    difrel(i) = difabs(i)/abs(exf(i));
%   end
% end
%  if (ifail(i) == -2) 
%    error(:,i)=0;
%    difabs(i)=0;
% end
% end
% 
% %TABLE CONSTRUCTION
% 
% T = tarray';
% FACAL = fzinv;
% FEX = exf';
% ESTREL = error(1,:)';
% RELERR=difrel;
% TRUNERR = error(3,:)';
% ESTABS=error(2,:)';
% ABSERR=difabs;
% N=ifzeval';
% IFAIL=ifail';
% 
% C = table(T,FACAL,FEX,ESTREL,RELERR,TRUNERR,ESTABS,ABSERR,N,IFAIL);

end