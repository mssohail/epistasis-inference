% Written: 31-Ma5, 2020
% Author: M Saqib Sohail


% function to convert upper trianglular matrix to a vector. Takes elements
% row-wise and places them in a long row vector that is n(n-1)/2 in length

function outVec = matUpperTriuToVec(mtx)

mtxDim = size(mtx,1);
vecLen = mtxDim*(mtxDim-1)/2;
outVec = zeros(1, vecLen);

counter = 0;
for i = 1:(mtxDim-1)
    counter = counter + 1;
    outVec(counter:(counter+mtxDim-(i + 1))) = mtx(i,i+1:end); 
    counter = counter + (mtxDim) - (i+1);
end
