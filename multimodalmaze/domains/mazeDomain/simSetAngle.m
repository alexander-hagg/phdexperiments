function bareFileName = simSetAngle(filename,theta,varargin)

fileExt = ''; if nargin > 2; fileExt = char(varargin{1}); end

bareFileName = [filename fileExt];

pathName = 'domains/mazeDomain/simulator/worlds/';
% Write to headless configuration file
xmlInFile = fullfile([pathName filename '.xml']);
xmlOutFile = [pathName filename fileExt '.xml'];

DOMnode = xmlread(xmlInFile);
fastsimElem = DOMnode.getElementsByTagName('fastsim');
robotElem = fastsimElem.item(0).getElementsByTagName('robot');
robotAttribs = robotElem.item(0).getAttributes;
thetaElem = robotAttribs.getNamedItem('theta');
thetaElem.setNodeValue(num2str(theta));
thetaElem = robotAttribs.getNamedItem('theta');
xmlwrite(xmlOutFile,DOMnode);
disp(['Written to ' xmlOutFile]);

% Write to headfull configuration file
xmlInVisFile = fullfile([pathName filename 'Vis.xml']);
xmlOutVisFile = [pathName filename fileExt 'Vis.xml'];

DOMnode = xmlread(xmlInVisFile);
fastsimElem = DOMnode.getElementsByTagName('fastsim');
robotElem = fastsimElem.item(0).getElementsByTagName('robot');
robotAttribs = robotElem.item(0).getAttributes;
thetaElem = robotAttribs.getNamedItem('theta');
thetaElem.setNodeValue(num2str(theta));
thetaElem = robotAttribs.getNamedItem('theta');
xmlwrite(xmlOutVisFile,DOMnode);
disp(['Written to ' xmlOutVisFile]);
end
