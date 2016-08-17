function sInd = getSourceInd(source)


switch source
    %Bottom edge (left to right)
    case 1;sInd=34;  case 2;sInd=38;  case 3;sInd=41;  case 4;sInd=45;
    %Right edge  (bottom to top)
    case 5;sInd=207;  case 6;sInd=311;  case 7;sInd=389;  case 8;sInd=493;
    %Top edge    (left to right)
    case 9;sInd=632;  case 10;sInd=636; case 11;sInd=639; case 12;sInd=643;
    %Left edge   (bottom to top)
    case 13;sInd=184; case 14;sInd=288; case 15;sInd=366; case 16;sInd=470;
end