function c=Nevronfijimatlab3(b)
% run fiji from inside matlab
% install the Nevron Macro v02 plugin in fiji/imagej
% see http://fiji.sc/Miji for more details

Miji

b=single(b);
MIJ.createImage('test',b,true)
MIJ.run('8-bit');
MIJ.run('Skeletonize (2D/3D)');
c=MIJ.getCurrentImage();
MIJ.run('Close')
MIJ.exit