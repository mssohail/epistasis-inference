
function [protListCell, aaList] = mapFromNTInd(testNTSites, protNamesCell, allProtIndMapingMtx)

protListCell = repmat({'     '}, length(testNTSites), 1);
aaList = zeros(length(testNTSites), 1);
for i = 1:length(testNTSites)
    thisTestNTSite = testNTSites(i);
    tempInd = find(allProtIndMapingMtx(1,:) == thisTestNTSite);
    if(isempty(tempInd))
        protListCell{i} = 'NonCoding';
        aaList(i) = -1;
    else
        protListCell{i} = protNamesCell{allProtIndMapingMtx(3,tempInd)};
        aaList(i) = allProtIndMapingMtx(2,tempInd);
    end
end