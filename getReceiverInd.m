function rInd = getReceiverInd(receivers)


rInd = [];
for k = receivers
    switch k
        %Bottom edge (left to right)
        case 1;  rInd = [rInd 34];  case 2;  rInd = [rInd 38];
        case 3;  rInd = [rInd 41];  case 4;  rInd = [rInd 45];
        %Right edge  (botom to top)
        case 5;  rInd = [rInd 207];  case 6;  rInd = [rInd 311];
        case 7;  rInd = [rInd 389];  case 8;  rInd = [rInd 493];
        %Top edge    (left to right)
        case 9;  rInd = [rInd 632];  case 10; rInd = [rInd 636];
        case 11; rInd = [rInd 639];  case 12; rInd = [rInd 643];
        %Left edge   (bottom to top)
        case 13; rInd = [rInd 184];  case 14; rInd = [rInd 288];
        case 15; rInd = [rInd 366];  case 16; rInd = [rInd 470];
    end
end