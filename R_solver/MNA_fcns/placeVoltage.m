function A = placeVoltage(A,Vnum,nodePos,nodeNeg)
A(nodePos,Vnum) = plus(A(nodePos,Vnum),1);
A(nodeNeg,Vnum) = plus(A(nodeNeg,Vnum),-1);
end