function boolStrongly = isStrongly(adjacencyMat)

%Is it Strongly Connected? Can't find any zero entries after multiplying
%the matrix by itself n times. n is the dimension of the square matrix
n = size(adjacencyMat,1);

resultaMat = adjacencyMat^n;

zerosSearch = find(resultaMat == 0);

boolStrongly = isempty(zerosSearch);

end