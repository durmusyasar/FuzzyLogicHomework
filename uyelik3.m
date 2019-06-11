function FS=uyelik3(N,S,P,X)
% Copyright(c): Doç. Dr. Ýsmail H. ALTAÞ, 2001  - uyelik3.m
% -------------- sub program for  fuzzification  
%                using TRIANGLE type membership functions
%
%          N         S        P                   S  
%         . \      /. \      /.                  /.\
%         .  \   /  .   \   /  .               /  . \
%         .   \/    .    \/    .             /    .  \
%         .  / \    .    / \   .           /      .   \
%     ------------------------------------SL-----ST---SR--------
%      N(1) S(1) N(3) S(2) P(1) S(3) P(3)          
%      N(2)                     P(2)
% -------------- Data to be supplied :
% X,  N=[NL NT NR],  S=[SL ST SR],   P=[PL PT PR]
% ---- resultant membership values are :  FS(n); for X
%      where n=1,2,3, ....,5
% ---- determine the membership values for X
	        if X <= N(2)
        	    FS(1) = 1;
                FS(2) = 0.0;
   		        FS(3) = 0.0;
	        elseif X <= S(1)
         	    FS(1) = abs((X-N(3)))/abs(N(2)-N(3));
                FS(2) = 0.0;
   		        FS(3) = 0.0;
             elseif X <= N(3)
                FS(1) = abs(X-N(3))/abs(N(2)-N(3));
                FS(2) = abs(X-S(1))/abs(S(2)-S(1));
   		        FS(3) = 0.0;
	        elseif X <= S(2)
                FS(1) = 0;
                FS(2) = abs(X-S(1))/abs(S(2)-S(1));
   		        FS(3) = 0.0;
	        elseif X <= P(1)
                FS(1) = 0;
                FS(2) = abs(X-S(3))/abs(S(2)-S(3));
   		        FS(3) = 0.0;
	        elseif X <= S(3)
        	    FS(1) = 0.0;
		        FS(2) = abs(X-S(3))/abs(S(2)-S(3));
	            FS(3) = abs(X-P(1))/abs(P(2)-P(1));
 	        elseif X <= P(2)
        	    FS(1) = 0.0;
		        FS(2) = 0.0;
	            FS(3) = abs(X-P(1))/abs(P(2)-P(1));
      	  else
        	    FS(1) = 0.0;
        	    FS(2) = 0.0;
        	    FS(3) = 1.0;
   		  end

	