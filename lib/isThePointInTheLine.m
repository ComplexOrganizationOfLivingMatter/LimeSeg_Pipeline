function [difference] = isThePointInTheLine(firstPointOfLine, secondPointOfLine, queryPoint)
%ISTHEPOINTINTHELINE Summary of this function goes here
%   Detailed explanation goes here

        n1 = (queryPoint(1) - firstPointOfLine(1))/(secondPointOfLine(1) - firstPointOfLine(1));
        n2 = (queryPoint(2) - firstPointOfLine(2))/(secondPointOfLine(2) - firstPointOfLine(2));
        n3 = (queryPoint(3) - firstPointOfLine(3))/(secondPointOfLine(3) - firstPointOfLine(3));
        difference = abs(n1 - n2) + abs(n1 - n3);
end

