#include "InheritOMesh.h"

class SingleHoleMesh : public SingleHole {
public:
	OMTriMesh holeMesh;
	OpenMesh::VPropHandleT<int> topDis;
	OpenMesh::VPropHandleT<double> meshDis;

	//std::vector<OpenMesh::VertexHandle> VertBoundaryHs;//boundary of origin mesh
	std::vector<OpenMesh::VertexHandle> VertBoundaryHolemeshHs;

	std::map<int, int> VertOri2NewId;//<origin, new> //, when extract maybe need another reverse map
	std::map<int, int> VertNew2OriId;//<new,origin>

	std::map<int, int> FacesOri2NewId;//<origin, new> 
	std::map<int, int> FacesNew2OriId;//<new,origin>

	int InitMesh(OMTriMesh & OriginMesh, int neighborRings, double distanceLimited);// extract
	int Mergeback(OMTriMesh & OriginMesh,bool updateOldVertecis = false);// and clean ?

	int FillHole(int type);
	int FillHoleCenterPoint();//0
	//int angle

	int SmoothMesh(bool changeOrigin);
private:
	int MoveToCenter(OpenMesh::Vec3f center);//move all point to a coordinate origin for calculate  accuracy
};

class HoleMeshHandle : public MeshHandle {
public:
	std::vector<SingleHoleMesh> holes;
	
	int ExtractHolesRegin(int neighborRings,double distanceLimited = 3.0);
	int FillHoles(int method);
	int MergerBack(bool updateOldVertecis = false);
protected:
private:
	int FindHoles();
};