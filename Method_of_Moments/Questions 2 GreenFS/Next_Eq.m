function Next = Next_Eq(E_q,q,x)
%  x      is the argument of the Exponential integral
%  E_q    value of the previous Exponential integral
%  q      Order of the new Exponential integral

    Next = (exp(-x)-x.*E_q)/q;

end