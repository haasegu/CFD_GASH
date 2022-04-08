
clear all
clc




model = createpde(4);
importGeometry(model,'NS_3D.stl');
pdegplot(model)

%pdegplot(model,'FaceLabels','on')

%mesh_Hmax = generateMesh(model,'Hmax',0.5);
mesh_Hmax = generateMesh(model,'Hmax',0.1,'GeometricOrder','linear');

[p,e,t] = meshToPet(mesh_Hmax);
pdemesh(mesh_Hmax)

ascii_write_mesh_3D( p, t, e, mfilename);