%  *****************************************************
%          TRIANGULAR membership function
%  *****************************************************
function uyelik=ucgen(a,b,c,x)
%
%*******************************************************
%                          *
%                        * . *
%                      *   .   *
%                    *     .     *
%                  *       .       *
%            ********************************
%                  a       b        c       ==> x
%
%                   a < b < c
%
%    inputs: the corners of the triangle: a, b, c, 
%    inputs: the variable x
%     AAA=[a b c x]
         mu1=(x-a)/(b-a);
         mu2=(c-x)/(c-b);
         mu3=0;
         uyelik=max(min(mu1,mu2),mu3);


