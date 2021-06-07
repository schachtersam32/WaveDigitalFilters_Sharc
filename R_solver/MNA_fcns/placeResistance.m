function Y = placeResistance(Y,G,node1,node2)

Y(node1,node1) = plus(Y(node1,node1),G);
Y(node1,node2) = plus(Y(node1,node2),-G);
Y(node2,node1) = plus(Y(node2,node1),-G);
Y(node2,node2) = plus(Y(node2,node2),G);
end
