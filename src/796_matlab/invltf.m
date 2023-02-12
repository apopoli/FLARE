function [fzinv, error, ifzeval, ifail] = invltf(tol, valt, fz, sigma0, ssbar, nmax)

%    ALGORITHM 796, COLLECTED ALGORITHMS FROM ACM.
%    THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
%    VOL. 25,NO. 3, September, 1999, P.  306--315.


% *********************************************************************

%  PURPOSE:
%  =======
%  THIS ROUTINE COMPUTES AN APPROXIMATE VALUE OF AN INVERSE LAPLACE TRANSFORM,
%       BY MEANS OF AN EXPANSION IN SINE AND COSINE FOURIER SERIES.
%       THE METHOD WHICH THIS SOFTWARE IS BASED ON UTILIZES THE Q-D ALGORITHM
%       TO ACCELERATE THE CONVERGENCE OF THE SERIES.   THE SUMMATION OF
%       SUCCESSIVE TERMS IS TRUNCATED WHEN THE ESTIMATED TRUNCATION ERROR
%       IS LESS OR EQUAL THAN AN INPUT PROVIDED TOLERANCE. THE DISCRETIZATION
%       ERROR IS MADE LESS OR EQUAL THAN THIS TOLERANCE BY SUITABLY CHOOSING
%       SOME METHOD'S PARAMETERS.

% INPUT PARAMETERS
% ================


% tol
%                       ON ENTRY TOL CONTAINS THE RELATIVE REQUIRED ACCURACY.

% 
% valt                  ON ENTRY, VALT CONTAINS A POSITIVE VALUE OF T
%                       FOR WHICH THE INVERSE LAPLACE TRANSFORM IS REQUIRED.
%                       VALT  HAS TO BE GREATER THAN ZERO.

% fz
%                       NAME OF THE FUNCTION SUBPROGRAMS FOR THE COMPLEX
%                       VALUED LAPLACE TRANSFORM TO BE INVERTED.
%                      
%                       FZ MUST BE DECLARED AS EXTERNAL IN THE PROGRAM
%                       FROM WHICH INVLTF IS CALLED.

% sigma0
%                       ON ENTRY, SIGMA0 CONTAINS THE VALUE OF THE ABSCISSA
%                       OF CONVERGENCE OF THE LAPLACE TRANSFORM FUNCTION TO
%                       BE INVERTED OR AN UPPER BOUND OF THIS.
%                       IT IS RECOMMENDED THAT A CORRECT VALUE TO SIGMA0 OR A
%                       CORRECT UPPER BOUND TO SIGMA0 IS PROVIDED.
%                       IF AN INCORRECT VALUE IS USED THE ROUTINE APPEARS TO
%                       WORK WELL BUT CONVERGES TO COMPLETELY WRONG RESULTS.
%                       THERE IS NO WAY IN WHICH THE ROUNTINE CAN DETECT THIS.

% ssbar
%                       ON ENTRY, IT SPECIFIES THE VALUE OF THE PARAMETER SS
%                       (> 2) TO BE USED IN CALCULATING THE PARAMETER D*T.
%                       TO OBTAIN DEFAULT OPTION (SS=4.1d0) ONE
%                       MAY SET   SSBAR = 0.

%
% nmax                  ON ENTRY, NMAX SPECIFIES THE MAXIMUM
%                       NUMBER OF FZ EVALUATIONS ALLOWED.


% OUTPUT PARAMETERS
% =================

% fzinv
%                      ON EXIT, FZINV CONTAINS  THE APPROXIMATION
%                      OF THE INVERSE LAPLACE TRANSFORM  AT THE POINT VALT.

% error                REAL  ARRAY OF DIMENSION 3.
%                      ON EXIT, ERROR(1) CONTAINS AN ESTIMATE OF THE RELATIVE ERROR
%                      WHICH SHOULD BE AN  UPPER BOUND FOR

%                      ABS(TRUEVALUE AT VALT - FZINV)/ABS(TRUEVALUE AT VALT)

%                      ON EXIT, ERROR(2) CONTAINS AN ESTIMATE OF THE ABSOLUTE ERROR
%                      WHICH SHOULD BE AN UPPER BOUND FOR

%                      ABS(TRUEVALUE AT VALT - FZINV)

%                      ON EXIT, ERROR(3) CONTAINS AN ESTIMATE OF THE TRUNCATION ERROR
%                      MADE IN  CALCULATING FZINV.

% ifzeval
%                      ON EXIT, IFZEVAL CONTAINS THE
%                      NUMBER OF EVALUATIONS OF FZ USED TO OBTAIN FZINV.

% ifail
%                      ON EXIT, A DIAGNOSTIC.
%                      = 0   ALL INPUTS WITHIN LIMITS.
%                            ALGORITHM APPERENTLY WORKED PROPERLY.
%                      =  1  TOL IS EQUAL OR GREATER THAN 1.
%                      =  2  VALT NOT POSITIVE.
%                      = -1  ACCURACY NOT REACHED AFTER NMAX FUNCTION EVALUATIONS.
%                      = -2  ALGORITHM DID NOT APPEAR TO BE CONVERGING, POSSIBLY
%                            DUE TO AN UNSUITABLE VALUE OF PARAMETER SSBAR.


%    *****************************************************************
%    INITIALIZATION
%    *****************************************************************

error=zeros(3,1);

ss=ssbar+4.1;

% ****************************************************************************

%  INPUT DATA CHECK

% ****************************************************************************


ifail=0;


if tol >= 1 
  ifail = 1;
 return 
end 
    
if valt < 0
  ifail = 2;
  return 
end 

%  SETTING THE PARAMETERS : ALPHA, TOL1, RELP
%  THE PRODUCT D*T IS COMPUTED IN SUCH A WAY THAT
%  THE ERROR IS LESS THAN TOL



%     RELP = X02AJF()
%     RELP  =0.222044d-15
%     FOR AN IBM 6000 RISC STATION MOD. 32H


if tol < eps 
    tol = eps;
end

tol1= tol*1E-1; % ORIGINALE: tol*10^-2;
alpha = 1.4;
relp=tol1;

dt = -log(tol1*2*atan(1))/(ss-log(10));



dt=dt*alpha;
gamma = dt/valt+ sigma0;
t = ss/2*valt;
pi_t =pi/t;

ma = 2; 
mp = ma-1;



%    CALL TO THE Q-D ALGORITHM



[fzinv,ifzeval,ifail,err1,err2,error]=qdacc(valt,fz,gamma,nmax,dt,pi_t,ifail,t,sigma0,tol1);


if ifail~=0 
    return
end


%            COMPUTATION OF THE TRUNCATION ERROR
% ERR1 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION AND THE PREVIOUS ONE,
% WHILE ERR2 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION
% AND THE LAST-THIRD


if (ifzeval>nmax)

if (err1 ~= 0 && err2 ~= 0) 
    if (err1 > relp && err2 > relp) 
      if (err1 > tol1 || err2 > tol1)
        ifail=-1;
      end
      
    end

    
elseif (err1 == 0 && err2 ~= 0) 
    if (err2 > relp || err2 > tol1) 
      ifail =-1;
      
    elseif (err2 == 0 && err1 ~= 0) 
      if (err1 > relp || err1 > tol1) 
        ifail =-1;
      end 
    end
end

end
%ESTIMATION OF THE ABSOLUTE ERROR, TRUNCATION ERROR AND THE SUM OF THEM


end

