function writePhenoToDisk(phenotypes,resolution,dir)
%WRITEPHENOTODISK Write phenotype images to disk

for i=1:length(phenotypes)

    pixelCoordinates = ceil(phenotypes{i}.Vertices*(resolution)/2)+resolution/2;
    pixelCoordinates(all(isnan(pixelCoordinates)'),:) = [];
    bw = poly2mask(pixelCoordinates(:,1),pixelCoordinates(:,2),resolution,resolution);
    imwrite(bw,[dir '/' int2str(i) '.png'],'png');
end
end

