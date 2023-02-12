function [fzinv,ifv,ifail,err1,err2,error]=qdacc(tv,fz,gamma,nmax,dt,pi_t,ifail,t,sigma0,tol1)


% ****************************************************************
%    PURPOSE:
%    ======
%    THIS SUBROUTINE IMPLEMENTS A COLUMN-VERSION OF THE Q-D ALGORITHM.
%    THE QD SCHEME IS BUILT UP PROCEEDING FROM LEFT TO RIGHT ALONG THE
%    DIAGONAL DIRECTION. THE ELEMENTS OF EACH DIAGONALS ARE COMPUTED
%    PROCEEDING BOTTOM-UP IN THE TABLE. THE LAST ELEMENT COMPUTED IN
%    EACH DIAGONAL IS THE COEFFICIENT OF THE CONTINUED FRACTION.

%    DESCRIPTION OF PARAMETERS:
%    =========================

%    INPUT PARAMETERS:

%   tv              ARRAY WITH TIME VALUE.

%   nmax            ON ENTRY, NMAX SPECIFIES THE MAXIMUM
%                   NUMBER OF FZ EVALUATIONS ALLOWED. 

%   dt               TIME STEP.
%

%   sigma0           ON ENTRY, SIGMA0 CONTAINS THE VALUE OF THE ABSCISSA
%                    OF CONVERGENCE OF THE LAPLACE TRANSFORM FUNCTION TO
%                    BE INVERTED OR AN UPPER BOUND OF THIS.
%                    IT IS RECOMMENDED THAT A CORRECT VALUE TO SIGMA0 OR A
%                    CORRECT UPPER BOUND TO SIGMA0 IS PROVIDED.
%                    IF AN INCORRECT VALUE IS USED THE ROUTINE APPEARS TO
%                    WORK WELL BUT CONVERGES TO COMPLETELY WRONG RESULTS.
%                    THERE IS NO WAY IN WHICH THE ROUNTINE CAN DETECT THIS.




%   fz               COMPLEX FUNCTION, SUPPLIED BY THE USER.
%                    FZ MUST RETURN THE VALUE OF THE LAPLACE TRANSFORM FUNCTION TO
%                    BE INVERTED, AT A GIVEN POINT. ITS SPECIFICATION IS:

%                    FZ MUST BE DECLARED AS EXTERNAL IN THE PROGRAM FROM WHICH INVLTF IS
%                    CALLED.


% gamma               IS OBTAINED AS dt/valt+ sigma0.

% ifail
%                      ON EXIT, A DIAGNOSTIC.
%                      = 0   ALL INPUTS WITHIN LIMITS.
%                            ALGORITHM APPERENTLY WORKED PROPERLY.
%                      =  1  TOL IS EQUAL OR GREATER THAN 1.
%                      =  2  VALT NOT POSITIVE.
%                      = -1  ACCURACY NOT REACHED AFTER NMAX FUNCTION EVALUATIONS.
%                      = -2  ALGORITHM DID NOT APPEAR TO BE CONVERGING, POSSIBLY
%                            DUE TO AN UNSUITABLE VALUE OF PARAMETER SSBAR.



%tol1                IS tol*10^-2 



     
%    OUTPUT PARAMETERS:

% ****************************************************************

%   fzinv            FZINV IS THE APPROXIMATED VALUE OF THE INVERSE LAPLACE
%                    TRANSFORM.

%   ifv              INTEGER. ON EXIT, IFV CONTAINS THE
%                    NUMBER OF EVALUATIONS OF FZ USED TO OBTAIN PART2.


%   ifail            INTEGER. ON EXIT, IFAIL CONTAINS POSSIBLE ERRORS DETECTED

%   err1             THE DIFFERENCE BETWEEN THE LAST APPROXIMATION AND THE
%                    PREVIOUS ONE.
%                   
%
%   err2             THE DIFFERENCE BETWEEN THE LAST APPROXIMATION.

%   error            STIMATION OF THE ABSOLUTE ERROR, TRUNCATION ERROR AND
%                    THE SUM OF THEM.





% **********************************************************
%   INITIALIZATION STEP                               
% **********************************************************



ma=2;
mp=ma-1;

%part0=0;           %PART0 CONTAINS THE APPROXIMATION
                    %OF THE LAST-THIRD ACCELERATED TERM.

part1=0;            %PART1 CONTAINS THE APPROXIMATION
                    %OF THE LAST-SECOND ACCELERATED TERM.
part2=0;            %PART2 CONTAINS THE APPROXIMATION
                    %OF THE LAST ACCELERATED.
                    
work=zeros(2,2*nmax+1); %COMPLEX  ARRAY OF DIMENSION(2,0:2*NMAX).  WORKSPACE AREAS.

first=1;
ifv=0;
err1=100;
err2=100;
fa3 = exp(dt)/t;

% STOPPING CRITERION

%     THE ALGORITHM STOPS WHEN BOTH ERR1 AND ERR2 ARE LESS THAN RELP
%                         OR
%     THE NUMBER OF FUNCTION EVALUATIONS, OF COURSE, IS >= NMAX
error=zeros(3,1);


while  (ifv < nmax && (( err1 > eps) && ( err2 > eps) && (err1 > tol1 || err2 > tol1)) || (((err1 == 0 && err2 ~= 0) && (err2 > eps || err2 > tol1)) || ((err2 == 0 && err1 ~= 0)&& (err1 >eps || err1 > tol1)))) 



if ifail  ~=0 
return
end

ibegin = 2*mp;
mm = 2*ma;


if(first==1)
ifv = ifv+1;



h = fz(gamma)/2;

mx = h;
work(2,1) = h;

% **************************************************************
%   INITIALIZE THE WORKSPACE AREA: THE ARRAYS Q and D OF
%   THE COLUMN VERSION OF THE QD ALGORITHM
%   THE VECTOR Q HAS BEEN STORED IN THE FIRST ROW OF THE ARRAY WORK
%   THE VECTOR D HAS BEEN STORED IN THE SECOND ROW OF THE ARRAY WORK
% ***************************************************************
for  k = 0:1
work(1,k+1) = mx;
mx = fz(gamma+1i*pi_t*(k+1));
ifv =ifv+1;
work(1,k+1) = mx/work(1,k+1);
work(2,k+2) =  work(1,k+1);
end

aux1 = work(1,2);
work(1,2) =  work(1,2) -work(1,1);
work(2,3) = - work(1,2);
work(1,1) = aux1;
work(2,2) = -work(2,2);
work(1,3) = 0;

end
u1 = (mm-1);
ifail=0;
%     DO 15 K = 3, UL
%       WORK(1,K) = 0.0D0
%15    CONTINUE

% *****************************************************************
% COMPUTATION OF THE COEFFICIENTS NEEDED IN THE QD ALGORITHM
% *****************************************************************
i=ibegin;

while(i <= u1 && ifail == 0 )
riv = true;
aux1 = work(1,1);
work(1,1) = mx;
mx = fz(gamma+1i*pi_t*(i+1));

work(1,1) = mx/work(1,1);
aux2 = work(1,2);
work(1,2) = work(1,1)-aux1;
ifv = ifv + 1;

j = 2;
while(j <= i && ifail == 0)
 if(j == i)
  work(1,j+1) = 0;
 end
aux3 = work(1,j+1);
if(riv)
if(aux2 ~=0)
work(1,j+1) = aux1*work(1,j)/aux2;
else
ifail = -2;
end
else
work(1,j+1) = work(1,j)-aux2+aux1;
end
aux1 = aux2;
aux2 = aux3;
riv = ~riv;
j =j+1;

end
 

work(2,i+2) = -work(1,i+1);
i =i+1;

end


if (first==0)

ibegin= (ibegin + 1);
end
if(first==1)
part0 = 0;
part1 = 0;

else
part0 = part1;
part1 = part2;
end


%*****************************************************************************
%  COMPUTATION OF THE  CONTINUED FRACTION BY THE  BACKWARD FORMULA
% *****************************************************************************

mm = 2*ma;

z = exp(1i*pi_t*tv);


h2m = ((work(2,mm)-work(2,mm+1))*z+1)/2;
r2m = sqrt((work(2,mm+1)*z/(h2m*h2m)) + 1);
r2m = 1-r2m;
r2m =-h2m*r2m;
cc= r2m*z;

for i=mm-1:-1:1
  
cc= (work(2,i+1)*z)/(1+cc);
 
end


cc = work(2,1)/(1+cc);
part2 = real(cc);
mp = ma;
ma = ma+1;
first=0;




er=[err1;err2];


err1 = fa3*abs(part2 - part1);
err2 = fa3*abs(part2 - part0);


tol1=tol1*abs(fa3*part2) + eps*abs(fa3*part2);

if tol1 <eps
    tol1 = eps;
elseif tol1>0
    tol1=1E-8; % ORIGINAL 10^-8
end  



end





%            COMPUTATION OF THE TRUNCATION ERROR
% ERR1 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION AND THE PREVIOUS ONE,
% WHILE ERR2 IS THE DIFFERENCE BETWEEN THE LAST APPROXIMATION
% AND THE LAST-THIRD




error(3) = err1+err2;
error(2) = tol1*abs(fa3*part2);
error(2) = error(2)+(err1+err2);
error(1) = tol1 + error(3);



%**********************************************************************

%     FZINV IS THE APPROXIMATED VALUE OF THE INVERSE LAPLACE TRANSFORM
 
% *********************************************************************


fzinv = exp(sigma0*tv)*fa3*part2;

   
end




