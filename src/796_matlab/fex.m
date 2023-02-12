function [fn_val] = fex(x)
%   CASE (1)
 fn_val = 1;

%  CASE (2)
%    fn_val = (1 -exp(-x))/(x*sqrt(pi*x));

%  CASE (3)
%    fn_val = 1/sqrt(pi*x);

%  CASE (4)
%   fn_val = x*cos(x);

%  CASE (5)
%    fn_val = x*exp(-x);

%  CASE (6)
%    fn_val = x;

%  CASE (7)
%    fn_val = sin(x);

%  CASE (8)
%    fn_val = exp(-0.5*x);


%  CASE (10)
%    fn_val = cos(2*sqrt(x))/sqrt(pi*x);


%  CASE (11)
%    fn_val = 2*exp(-4/x)/(x*sqrt(pi*x));

%  CASE (12)
%    fn_val = sin(x)/x;

%  CASE (13)
%    fn_val = exp(-0.2*x)*sin(x);

%  CASE (14)
%    fn_val = 0.5*x^2;

%  CASE (15)
%    if (x > 2) 
%     fn_val = 1;

%   elseif(x < 2) 
%     fn_val = 0;

%   else
%    fn_val = 0.5;
%   end



%  CASE (17)
%    fn_val = 2/sqrt(3)*exp(-x/2)*sin(x*sqrt(3)/2);

%  CASE (18)
%   fn_val = sinh(3*x);

%  CASE (19)
%    fn_val = x^5;

%  CASE (20)
%    fn_val = x/2*sin(x);

%  CASE (21)
%    fn_val = exp(-x) - exp(-1000*x);

%  CASE (22)
%    fn_val = cos(x);

%  CASE (23)
%    fn_val = x*exp(x/4);

%  CASE (24)
%    fn_val = 2*sqrt(x/pi);

%  CASE (25)
%    fn_val =exp(-x)/sqrt(pi*x);

%  CASE (26)
%    fn_val = (1+ 4*x)/sqrt(pi*x);

%  CASE (27)
%    fn_val = (sin(x) - x*cos(x))/2;

%  CASE (28)
%    fn_val = 1 - exp(-x)*(1 + x);

%  CASE (29)
%    fn_val = exp(-x)/12*(exp(3*x) - cos(sqrt(3)*x) -sqrt(3)*sin(sqrt(3)*x));



%  CASE (30)
%    fn_val = 2* (cos(2*x) - cos(x))/x;

%  CASE (31)
%    fn_val = (1- exp(-x))/x;



%  CASE(32)
%eulero = 0.5772156649015329;

%fn_val = -eulero - log(x);

%  CASE (33)
%    if (x >= 0 && x <= 1) 
%      fn_val = x;
%    else
%      fn_val = 1;
%    end







end


