function [z]=Uvsolve_row(M,v,z,n)

%
% Inner product version of in-place solving Uz=v
% with upper triangular U 
% stored in a 2d array M and vector v stored in a 1d array. M need not match dimensions with U 
% if done as a real code then M would be declared with its
% dimensions Mrdim and Mcdim (or it would come with them in an object)
% here they are extracted for reference and debugging.
% The same is true for the 1d array storing v and the 1d output array 
% storing the vector z.
% n is the dimension of U, v, and z 
% with U stored in the leading block of M in the usual fashion.
%
% This version would be preferred in a language that
% stores two dimensional arrays in a row major ordering.
% In a language that uses column major ordering for 
% two dimensional arrays, a column oriented middle product
% form of algorithm would be preferred.
%
% The output vector z  would be declared in a real code.
% It is passed here as an input variable to indicate
% how many languages would provide a data structure 
% to be updated with the computed output, e.g.,
% Here its dimensions are not needed but they would be passed in
% and used in that declaration.
% Here they are extracted for reference and debugging.
%
%
%
% The following extract the dimensions
% of the M, v and z arrays from the information
% Matlab maintains for them in their particular
% objects. They would be used in the declarations
% of M, v and z mentioned above in a code with 
% a compiled and typed language.
% 

  Mrdim=size(M,1);
  Mcdim=size(M,2); 
  zdim=size(z,1);
  vdim=size(v,1);  
%
% note the code relies on Matlab suppressing operations on 
% vectors of length 0, i.e., when index range r:s is such that r > s
%


%
% the i loop goes over the elements of the output solution vector z
%
% Row i of U has a nonzero in the (i,i) and all of its possibly nonzero
% elements are in positions (i,i+1:n).
% The last row is U(n,n)*e_n^T.


  z(n)=v(n)/M(n,n);
  for i = n-1:-1:1
    z(i) = (v(i)-M(i,i+1:n)*z(i+1:n))/M(i,i);
  end % end i loop
return