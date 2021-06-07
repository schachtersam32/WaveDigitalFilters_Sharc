function [B,A] = placeNullor(B,A,NullNode)
B(NullNode(1)) = plus(B(NullNode(1)),1);
B(NullNode(2)) = plus(B(NullNode(2)),-1);

A(NullNode(3)) = plus(A(NullNode(3)),1);
A(NullNode(4)) = plus(A(NullNode(4)),-1);

end